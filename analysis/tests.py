from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

from django.core.management import call_command
from django.core.management.base import CommandError
from django.db import DatabaseError, connection
from django.test import TestCase, TransactionTestCase

from analysis.services.recurrent_axis.project_match import (
    build_workflow_axis_identity,
    enrich_axis_final_response_with_project_matches,
)
from database.models import (
    AxisContextPresence,
    AxisDatasetContext,
    AxisFinalEvidence,
    AxisFinalRecurrentSummary,
    AxisObservation,
    AxisRecurrentSummary,
    AxisResultArtifact,
    AxisSignatureProjectIndex,
    AxisStructureRecurrentSummary,
    CanonicalAxis,
    DatasetAxisFinalOccurrence,
    DatasetAxisFinalProject,
    DatasetMetadata,
    SpongeEvidence,
)


class ClearRecurrentAxisTests(TransactionTestCase):
    imported_counts = {
        AxisDatasetContext: 1,
        AxisResultArtifact: 3,
        CanonicalAxis: 2,
        AxisObservation: 3,
        AxisFinalEvidence: 2,
        SpongeEvidence: 1,
        AxisContextPresence: 1,
        AxisStructureRecurrentSummary: 1,
        AxisFinalRecurrentSummary: 1,
    }

    def setUp(self):
        self.metadata = DatasetMetadata.objects.create(
            dataset="TCGA_ACC_mRNA",
            programme="TCGA",
            obs_type="sample",
            cancer_type="ACC",
            gene_bio_type="mRNA",
        )
        project = DatasetAxisFinalProject.objects.create(
            source="TCGA", module="module2", dataset_name=self.metadata.dataset
        )
        signature = "mRNA-miRNA|hsa-miR-210|ESPL1||"
        occurrence = DatasetAxisFinalOccurrence.objects.create(
            project=project, row_index=0, axis_id="legacy", axis_signature=signature
        )
        self.legacy_rows = (
            project,
            occurrence,
            AxisSignatureProjectIndex.objects.create(
                project=project,
                occurrence=occurrence,
                axis_signature=signature,
                source=project.source,
                module=project.module,
                dataset_name=project.dataset_name,
            ),
            AxisRecurrentSummary.objects.create(axis_signature=signature),
        )

        # Import two versions into the same context to exercise historical rows.
        self.import_module2_result("CDCA3")
        self.import_module2_result("ESPL1")

    @staticmethod
    def import_module2_result(axis_final_mrna):
        with TemporaryDirectory() as root:
            dataset_dir = Path(root) / "TCGA_ACC"
            dataset_dir.mkdir()
            (dataset_dir / "TCGA_ACC_ceRNA_axis_final.csv").write_text(
                f"mRNA,miRNA,axis_regulation\n{axis_final_mrna},hsa-miR-210,up_down\n",
                encoding="utf-8",
            )
            (dataset_dir / "TCGA_ACC_sponge_result.csv").write_text(
                "mRNA,miRNA,cor,pcor,mscor\nESPL1,hsa-miR-210,0.7,0.5,0.2\n",
                encoding="utf-8",
            )
            call_command(
                "import_recurrent_axis", "--module2-root", root, stdout=StringIO()
            )

    def assert_row_counts(self, expected_counts):
        for model, count in expected_counts.items():
            with self.subTest(table=model._meta.db_table):
                self.assertEqual(model.objects.count(), count)

    def test_default_and_explicit_dry_runs_report_counts_without_deleting(self):
        for arguments in ((), ("--dry-run",)):
            with self.subTest(arguments=arguments):
                output = StringIO()
                call_command("clear_recurrent_axis", *arguments, stdout=output)

                self.assertIn("Dry run on database 'default'", output.getvalue())
                for model, count in self.imported_counts.items():
                    self.assertIn(f"{model._meta.db_table}: {count}", output.getvalue())
                self.assert_row_counts(self.imported_counts)

    def test_confirm_clears_all_nine_tables_and_preserves_shared_and_early_data(self):
        self.assert_row_counts(self.imported_counts)
        self.assertEqual(AxisResultArtifact.objects.filter(is_active=False).count(), 1)

        # Repeating the same command after the database is empty also succeeds.
        for _ in range(2):
            output = StringIO()
            call_command("clear_recurrent_axis", "--confirm", stdout=output)
            self.assertIn("Cleared 9 recurrent-axis tables", output.getvalue())
            self.assert_row_counts(dict.fromkeys(self.imported_counts, 0))

        self.metadata.refresh_from_db()
        self.assertEqual(self.metadata.dataset, "TCGA_ACC_mRNA")
        for row in self.legacy_rows:
            row.refresh_from_db()

    def test_database_failure_rolls_back_the_clear(self):
        def fail_after_deleting(execute, sql, params, many, context):
            result = execute(sql, params, many, context)
            if sql.lstrip().upper().startswith(("TRUNCATE ", "DELETE FROM ")):
                raise DatabaseError("simulated failure after deletion")
            return result

        with connection.execute_wrapper(fail_after_deleting):
            with self.assertRaisesMessage(CommandError, "simulated failure after deletion"):
                call_command("clear_recurrent_axis", "--confirm", stdout=StringIO())

        self.assert_row_counts(self.imported_counts)

    def test_clear_allows_import_to_rebuild_evidence_and_summaries(self):
        call_command("clear_recurrent_axis", "--confirm", stdout=StringIO())
        self.import_module2_result("ESPL1")

        self.assertEqual(AxisResultArtifact.objects.filter(is_active=True).count(), 2)
        self.assertEqual(AxisFinalEvidence.objects.count(), 1)
        self.assertEqual(SpongeEvidence.objects.count(), 1)
        summary = AxisStructureRecurrentSummary.objects.get()
        self.assertEqual(summary.dataset_count, 1)
        self.assertEqual(summary.both_result_context_count, 1)
        self.assertEqual(
            AxisFinalRecurrentSummary.objects.get().dominant_axis_regulation,
            "up_down",
        )


class RecurrentAxisProjectMatchTests(TestCase):
    workflow_axis = {"miRNA": "hsa-miR-210", "mRNA": "ESPL1"}
    second_workflow_axis = {"miRNA": "hsa-miR-210", "mRNA": "CDCA3"}

    @classmethod
    def setUpTestData(cls):
        cls.axis = cls.create_axis(cls.workflow_axis)
        cls.second_axis = cls.create_axis(cls.second_workflow_axis)

    @staticmethod
    def create_axis(record):
        identity = build_workflow_axis_identity(record, row_index=0)
        return CanonicalAxis.objects.create(
            **record,
            axis_key=identity["axis_key"],
            axis_signature=identity["axis_signature"],
            axis_type=identity["canonical_axis_type"],
        )

    @staticmethod
    def create_dataset(name, source, obs_type):
        return DatasetMetadata.objects.create(
            dataset=name, programme=source, obs_type=obs_type, gene_bio_type="mRNA"
        )

    def create_context(self, metadata, source, module, **group):
        context = AxisDatasetContext.objects.create(
            dataset_metadata=metadata, dataset_source=source, module=module, **group
        )
        self.add_match(context, self.axis)
        return context

    @staticmethod
    def add_match(context, axis):
        AxisContextPresence.objects.create(
            context=context,
            axis=axis,
            has_axis_final=True,
            observation_count=1,
            axis_final_observation_count=1,
        )

    def assert_group_matches(self, source, obs_type, group_by, group_values):
        metadata = self.create_dataset(f"{source}_test", source, obs_type)
        contexts = [
            self.create_context(
                metadata, source, "module3",
                group_type="other", group_by=group_by, group_value=value,
            )
            for value in group_values
        ]

        response = enrich_axis_final_response_with_project_matches(response_data={
            "results": [self.workflow_axis],
        })
        summary = response["axis_reference_match_summary"]
        self.assertEqual(summary["matched_dataset_count"], 1)
        self.assertEqual(summary["matched_context_count"], 2)
        self.assertEqual(summary["axis_context_match_count"], 2)
        self.assertEqual(summary["matched_project_count"], 2)
        self.assertEqual(summary["project_hits"], summary["context_hits"])
        self.assertEqual(response["axis_project_match_summary"], summary)

        record = response["results"][0]
        self.assertEqual(record["reference_context_match_count"], 2)
        for matches in (summary["context_hits"], record["reference_context_matches"]):
            self.assertEqual(
                {match["context_id"]: match["group_value"] for match in matches},
                {context.id: context.group_value for context in contexts},
            )
            self.assertEqual(
                {match["project_key"] for match in matches},
                {f"{source}|module3|{metadata.dataset}|other|{group_by}"},
            )

    def test_sc_group_values_count_as_distinct_contexts(self):
        self.assert_group_matches("SC", "cell", "cell_type", ("T cell", "B cell"))

    def test_st_group_values_count_as_distinct_contexts(self):
        self.assert_group_matches("ST", "spot", "domain", ("region 1", "region 2"))

    def test_bulk_context_counts_and_distinct_axis_counts(self):
        tcga = self.create_dataset("TCGA_ACC_mRNA", "TCGA", "sample")
        tcga_context = self.create_context(tcga, "TCGA", "module2")
        timedb = self.create_dataset("GSE20194", "TIMEDB", "sample")
        grade = self.create_context(
            timedb, "TIMEDB", "module3", group_type="grade", group_by="Grade"
        )
        stage = self.create_context(
            timedb, "TIMEDB", "module3", group_type="stage", group_by="Stage"
        )
        self.add_match(grade, self.second_axis)

        response = enrich_axis_final_response_with_project_matches(response_data={
            "results": [self.workflow_axis, self.workflow_axis, self.second_workflow_axis],
        })
        summary = response["axis_reference_match_summary"]
        self.assertEqual(summary["matched_dataset_count"], 2)
        self.assertEqual(summary["matched_context_count"], 3)
        self.assertEqual(summary["axis_context_match_count"], 4)
        self.assertEqual(summary["matched_record_count"], 3)
        self.assertEqual(summary["matched_unique_axis_count"], 2)
        self.assertEqual(
            {hit["context_id"]: hit["matched_axis_count"] for hit in summary["context_hits"]},
            {tcga_context.id: 1, grade.id: 2, stage.id: 1},
        )
