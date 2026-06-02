import pandas as pd

from django.core.management.base import BaseCommand
from django.db import transaction

from database.models import RNANode, InteractionDatabase, RNAInteraction


class Command(BaseCommand):
    help = "Import RNA interaction data from CSV"

    def add_arguments(self, parser):
        parser.add_argument("csv_path", type=str)
        parser.add_argument("--chunksize", type=int, default=100000)
        parser.add_argument("--batch-size", type=int, default=5000)

    def handle(self, *args, **options):
        csv_path = options["csv_path"]
        chunksize = options["chunksize"]
        batch_size = options["batch_size"]

        self.stdout.write("Step 1: importing databases...")
        self.import_databases(csv_path, chunksize, batch_size)

        self.stdout.write("Step 2: importing RNA nodes...")
        self.import_nodes(csv_path, chunksize, batch_size)

        self.stdout.write("Step 3: importing interactions and M2M databases...")
        self.import_interactions(csv_path, chunksize, batch_size)

        self.stdout.write(self.style.SUCCESS("Import finished."))

    def read_csv_chunks(self, csv_path, chunksize):
        return pd.read_csv(
            csv_path,
            chunksize=chunksize,
            dtype=str,
            keep_default_na=False,
        )

    def parse_target_type(self, interaction_type):
        """
        miRNA-mRNA -> mRNA
        miRNA-lncRNA -> lncRNA
        """
        if not interaction_type or "-" not in interaction_type:
            return "unknown"

        target_type = interaction_type.split("-", 1)[1]

        valid_types = {"miRNA", "mRNA", "lncRNA", "circRNA"}
        return target_type if target_type in valid_types else "unknown"

    def split_databases(self, database_value):
        if not database_value:
            return []

        return [
            x.strip()
            for x in database_value.split(";")
            if x.strip()
        ]

    def import_databases(self, csv_path, chunksize, batch_size):
        database_names = set()

        for chunk in self.read_csv_chunks(csv_path, chunksize):
            for value in chunk["database"].dropna().unique():
                database_names.update(self.split_databases(value))

        objs = [
            InteractionDatabase(name=name)
            for name in database_names
        ]

        InteractionDatabase.objects.bulk_create(
            objs,
            batch_size=batch_size,
            ignore_conflicts=True,
        )

        self.stdout.write(f"Databases imported: {len(database_names)}")

    def import_nodes(self, csv_path, chunksize, batch_size):
        for i, chunk in enumerate(self.read_csv_chunks(csv_path, chunksize), start=1):
            node_keys = set()

            for row in chunk.itertuples(index=False):
                mirna = row.miRNA
                cerna = row.ceRNA
                species = row.species
                interaction_type = row.type

                if mirna:
                    node_keys.add((mirna, "miRNA", species))

                if cerna:
                    target_type = self.parse_target_type(interaction_type)
                    node_keys.add((cerna, target_type, species))

            objs = [
                RNANode(
                    name=name,
                    rna_type=rna_type,
                    species=species,
                )
                for name, rna_type, species in node_keys
            ]

            RNANode.objects.bulk_create(
                objs,
                batch_size=batch_size,
                ignore_conflicts=True,
            )

            self.stdout.write(f"Node chunk {i} imported, unique nodes in chunk: {len(node_keys)}")

    def import_interactions(self, csv_path, chunksize, batch_size):
        through_model = RNAInteraction.databases.through

        database_map = {
            db.name: db.id
            for db in InteractionDatabase.objects.all()
        }

        for i, chunk in enumerate(self.read_csv_chunks(csv_path, chunksize), start=1):
            chunk = chunk[["miRNA", "ceRNA", "species", "database", "type"]].copy()

            source_node_keys = set()
            target_node_keys = set()

            for row in chunk.itertuples(index=False):
                species = row.species
                target_type = self.parse_target_type(row.type)

                source_node_keys.add((row.miRNA, "miRNA", species))
                target_node_keys.add((row.ceRNA, target_type, species))

            names = {
                key[0]
                for key in source_node_keys | target_node_keys
                if key[0]
            }

            species_set = {
                key[2]
                for key in source_node_keys | target_node_keys
                if key[2]
            }

            nodes = RNANode.objects.filter(
                name__in=names,
                species__in=species_set,
            ).only("id", "name", "rna_type", "species")

            node_map = {
                (node.name, node.rna_type, node.species): node.id
                for node in nodes
            }

            interaction_objs = []
            interaction_keys = set()

            for row in chunk.itertuples(index=False):
                species = row.species
                interaction_type = row.type
                target_type = self.parse_target_type(interaction_type)

                source_id = node_map.get((row.miRNA, "miRNA", species))
                target_id = node_map.get((row.ceRNA, target_type, species))

                if not source_id or not target_id:
                    continue

                key = (source_id, target_id, species, interaction_type)

                if key in interaction_keys:
                    continue

                interaction_keys.add(key)

                interaction_objs.append(
                    RNAInteraction(
                        source_id=source_id,
                        target_id=target_id,
                        species=species,
                        interaction_type=interaction_type,
                    )
                )

            with transaction.atomic():
                RNAInteraction.objects.bulk_create(
                    interaction_objs,
                    batch_size=batch_size,
                    ignore_conflicts=True,
                )

                source_ids = {key[0] for key in interaction_keys}
                target_ids = {key[1] for key in interaction_keys}
                interaction_types = {key[3] for key in interaction_keys}
                species_values = {key[2] for key in interaction_keys}

                interactions = RNAInteraction.objects.filter(
                    source_id__in=source_ids,
                    target_id__in=target_ids,
                    species__in=species_values,
                    interaction_type__in=interaction_types,
                ).only("id", "source_id", "target_id", "species", "interaction_type")

                interaction_map = {
                    (
                        obj.source_id,
                        obj.target_id,
                        obj.species,
                        obj.interaction_type,
                    ): obj.id
                    for obj in interactions
                }

                m2m_objs = []
                m2m_keys = set()

                for row in chunk.itertuples(index=False):
                    species = row.species
                    interaction_type = row.type
                    target_type = self.parse_target_type(interaction_type)

                    source_id = node_map.get((row.miRNA, "miRNA", species))
                    target_id = node_map.get((row.ceRNA, target_type, species))

                    if not source_id or not target_id:
                        continue

                    interaction_id = interaction_map.get(
                        (source_id, target_id, species, interaction_type)
                    )

                    if not interaction_id:
                        continue

                    for db_name in self.split_databases(row.database):
                        db_id = database_map.get(db_name)

                        if not db_id:
                            continue

                        m2m_key = (interaction_id, db_id)

                        if m2m_key in m2m_keys:
                            continue

                        m2m_keys.add(m2m_key)

                        m2m_objs.append(
                            through_model(
                                rnainteraction_id=interaction_id,
                                interactiondatabase_id=db_id,
                            )
                        )

                through_model.objects.bulk_create(
                    m2m_objs,
                    batch_size=batch_size,
                    ignore_conflicts=True,
                )

            self.stdout.write(
                f"Interaction chunk {i} imported, "
                f"interactions: {len(interaction_objs)}, "
                f"m2m rows: {len(m2m_objs)}"
            )
