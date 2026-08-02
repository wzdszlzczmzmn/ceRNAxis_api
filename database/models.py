from django.core.exceptions import ValidationError
from django.db import models
from django.db.models import Q


class RNANode(models.Model):
    RNA_TYPE_CHOICES = [
        ("miRNA", "miRNA"),
        ("mRNA", "mRNA"),
        ("lncRNA", "lncRNA"),
        ("circRNA", "circRNA"),
        ("unknown", "unknown"),
    ]

    name = models.CharField(max_length=255)
    rna_type = models.CharField(
        max_length=20,
        choices=RNA_TYPE_CHOICES,
        default="unknown",
    )
    species = models.CharField(
        max_length=100,
        default="Homo sapiens",
    )

    class Meta:
        db_table = "rna_node"
        constraints = [
            models.UniqueConstraint(
                fields=["name", "rna_type", "species"],
                name="unique_rna_node",
            )
        ]
        indexes = [
            models.Index(fields=["name"]),
            models.Index(fields=["rna_type"]),
            models.Index(fields=["species"]),
            models.Index(fields=["name", "species"]),
        ]

    def __str__(self):
        return f"{self.name} ({self.rna_type}, {self.species})"


class InteractionDatabase(models.Model):
    name = models.CharField(max_length=100, unique=True)

    class Meta:
        db_table = "interaction_database"
        indexes = [
            models.Index(fields=["name"]),
        ]

    def __str__(self):
        return self.name


class RNAInteraction(models.Model):
    source = models.ForeignKey(
        RNANode,
        on_delete=models.CASCADE,
        related_name="outgoing_interactions",
    )
    target = models.ForeignKey(
        RNANode,
        on_delete=models.CASCADE,
        related_name="incoming_interactions",
    )

    species = models.CharField(
        max_length=100,
        default="Homo sapiens",
    )

    interaction_type = models.CharField(max_length=50)

    databases = models.ManyToManyField(
        InteractionDatabase,
        related_name="interactions",
        blank=True,
    )

    class Meta:
        db_table = "rna_interaction"
        constraints = [
            models.UniqueConstraint(
                fields=["source", "target", "species", "interaction_type"],
                name="unique_rna_interaction",
            )
        ]
        indexes = [
            models.Index(fields=["source"]),
            models.Index(fields=["target"]),
            models.Index(fields=["species"]),
            models.Index(fields=["interaction_type"]),
            models.Index(fields=["source", "target"]),
        ]

    def __str__(self):
        return (
            f"{self.source.name} -> {self.target.name} "
            f"({self.interaction_type})"
        )


class FilterField(models.Model):
    table_name = models.CharField(max_length=100)
    field_name = models.CharField(max_length=100)
    field_label = models.CharField(max_length=255)

    sort_order = models.PositiveIntegerField(default=0)
    is_active = models.BooleanField(default=True)

    class Meta:
        db_table = "filter_field"
        constraints = [
            models.UniqueConstraint(
                fields=["table_name", "field_name"],
                name="unique_filter_field",
            )
        ]
        indexes = [
            models.Index(fields=["table_name"]),
            models.Index(fields=["table_name", "is_active"]),
            models.Index(fields=["table_name", "field_name"]),
        ]

    def __str__(self):
        return f"{self.table_name}.{self.field_name}"


class FilterOption(models.Model):
    field = models.ForeignKey(
        FilterField,
        on_delete=models.CASCADE,
        related_name="options",
    )

    value = models.CharField(max_length=255)

    sort_order = models.PositiveIntegerField(default=0)
    is_active = models.BooleanField(default=True)

    class Meta:
        db_table = "filter_option"
        constraints = [
            models.UniqueConstraint(
                fields=["field", "value"],
                name="unique_filter_option",
            )
        ]
        indexes = [
            models.Index(fields=["field"]),
            models.Index(fields=["field", "is_active"]),
            models.Index(fields=["field", "sort_order"]),
        ]

    def __str__(self):
        return f"{self.field}: {self.value}"


class DatasetMetadata(models.Model):
    GENE_BIO_TYPE_CHOICES = [
        ("miRNA", "miRNA"),
        ("mRNA", "mRNA"),
        ("lncRNA", "lncRNA"),
        ("circRNA", "circRNA"),
    ]

    OBS_TYPE_CHOICES = [
        ("sample", "sample"),
        ("cell", "cell"),
        ("spot", "spot"),
    ]

    dataset = models.CharField(max_length=100, unique=True)

    programme = models.CharField(max_length=50)

    obs_type = models.CharField(
        max_length=50,
        choices=OBS_TYPE_CHOICES,
        db_index=True,
    )

    reference = models.CharField(
        max_length=50,
        blank=True,
        default="",
    )

    cancer_type = models.CharField(max_length=50)

    cancer_type_full_name = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    gene_bio_type = models.CharField(
        max_length=20,
        choices=GENE_BIO_TYPE_CHOICES,
        db_index=True,
    )

    workflow = models.CharField(
        max_length=100,
        blank=True,
        default="",
    )

    sample_nums = models.PositiveIntegerField(default=0)
    cell_nums = models.PositiveIntegerField(default=0)
    spot_nums = models.PositiveIntegerField(default=0)

    class Meta:
        db_table = "dataset_metadata"
        indexes = [
            models.Index(fields=["programme"]),
            models.Index(fields=["obs_type"]),
            models.Index(fields=["cancer_type"]),
            models.Index(fields=["gene_bio_type"]),
            models.Index(fields=["reference"]),
            models.Index(fields=["programme", "cancer_type"]),
            models.Index(fields=["cancer_type", "gene_bio_type"]),
            models.Index(fields=["programme", "obs_type"]),
            models.Index(fields=["programme", "obs_type", "gene_bio_type"]),
            models.Index(fields=["obs_type", "gene_bio_type"]),
        ]

    def __str__(self):
        return self.dataset


class DatasetAxisFinalProject(models.Model):
    class Source(models.TextChoices):
        TCGA = "TCGA", "TCGA"
        TIMEDB = "TIMEDB", "TIMEDB"

    class Module(models.TextChoices):
        MODULE2 = "module2", "Module 2"
        MODULE3 = "module3", "Module 3"

    class GroupType(models.TextChoices):
        NONE = "none", "None"
        COMMON = "common", "Common"
        GRADE = "grade", "Grade"
        STAGE = "stage", "Stage"

    source = models.CharField(
        max_length=20,
        choices=Source.choices,
        db_index=True,
    )

    module = models.CharField(
        max_length=20,
        choices=Module.choices,
        db_index=True,
    )

    dataset_name = models.CharField(
        max_length=128,
        db_index=True,
    )

    group_type = models.CharField(
        max_length=20,
        choices=GroupType.choices,
        default=GroupType.NONE,
        db_index=True,
    )

    group_by = models.CharField(
        max_length=128,
        blank=True,
        default="",
        db_index=True,
    )

    annotation_dir_name = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    annotation_file_prefix = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    axis_final_file_name = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    axis_final_file_path = models.TextField(
        blank=True,
        default="",
    )

    file_sha256 = models.CharField(
        max_length=64,
        blank=True,
        default="",
        db_index=True,
    )

    row_count = models.PositiveIntegerField(default=0)

    is_active = models.BooleanField(
        default=True,
        db_index=True,
    )

    imported_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "dataset_axis_final_project"
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "source",
                    "module",
                    "dataset_name",
                    "group_type",
                    "group_by",
                ],
                name="uniq_dataset_axis_project",
            ),
        ]
        indexes = [
            models.Index(
                fields=["source", "module", "dataset_name"],
                name="idx_dataset_axis_project_base",
            ),
            models.Index(
                fields=["dataset_name", "group_type"],
                name="idx_dataset_axis_project_group",
            ),
        ]

    def __str__(self):
        if self.group_type != self.GroupType.NONE:
            return (
                f"{self.source}:{self.module}:"
                f"{self.dataset_name}:{self.group_type}:{self.group_by}"
            )

        return f"{self.source}:{self.module}:{self.dataset_name}"


class DatasetAxisFinalOccurrence(models.Model):
    project = models.ForeignKey(
        DatasetAxisFinalProject,
        on_delete=models.CASCADE,
        related_name="axis_occurrences",
    )

    row_index = models.PositiveIntegerField()

    axis_signature = models.CharField(
        max_length=1024,
        db_index=True,
        help_text="axis_type|miRNA|mRNA|lncRNA|circRNA",
    )

    axis_id = models.CharField(
        max_length=512,
        db_index=True,
    )

    axis_type = models.CharField(
        max_length=128,
        blank=True,
        default="",
        db_index=True,
    )

    axis_regulation = models.CharField(
        max_length=64,
        blank=True,
        default="",
        db_index=True,
    )

    miRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    mRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    lncRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    circRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    mRNA_log2FC = models.FloatField(null=True, blank=True)
    mRNA_regulation = models.CharField(
        max_length=32,
        blank=True,
        default="",
        db_index=True,
    )

    miRNA_log2FC = models.FloatField(null=True, blank=True)
    miRNA_regulation = models.CharField(
        max_length=32,
        blank=True,
        default="",
        db_index=True,
    )

    lncRNA_log2FC = models.FloatField(null=True, blank=True)
    lncRNA_regulation = models.CharField(
        max_length=32,
        blank=True,
        default="",
        db_index=True,
    )

    circRNA_log2FC = models.FloatField(null=True, blank=True)
    circRNA_regulation = models.CharField(
        max_length=32,
        blank=True,
        default="",
        db_index=True,
    )

    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "dataset_axis_final_occurrence"
        constraints = [
            models.UniqueConstraint(
                fields=["project", "axis_signature"],
                name="uniq_project_axis_signature",
            ),
        ]
        indexes = [
            models.Index(
                fields=["axis_signature"],
                name="idx_axis_occ_signature",
            ),
            models.Index(
                fields=["project", "axis_signature"],
                name="idx_axis_occ_project_sig",
            ),
            models.Index(
                fields=["miRNA", "mRNA"],
                name="idx_axis_occ_mirna_mrna",
            ),
            models.Index(
                fields=["miRNA", "mRNA", "lncRNA"],
                name="idx_axis_occ_lncrna",
            ),
            models.Index(
                fields=["miRNA", "mRNA", "circRNA"],
                name="idx_axis_occ_circrna",
            ),
        ]

    def __str__(self):
        return f"{self.project_id}:{self.axis_signature}"


class AxisSignatureProjectIndex(models.Model):
    """
    Precomputed index:
        axis_signature -> reference projects

    This table can be rebuilt from DatasetAxisFinalOccurrence.
    """

    axis_signature = models.CharField(
        max_length=1024,
        db_index=True,
    )

    project = models.ForeignKey(
        DatasetAxisFinalProject,
        on_delete=models.CASCADE,
        related_name="axis_signature_indexes",
    )

    occurrence = models.ForeignKey(
        DatasetAxisFinalOccurrence,
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="signature_index_rows",
    )

    source = models.CharField(
        max_length=20,
        db_index=True,
    )

    module = models.CharField(
        max_length=20,
        db_index=True,
    )

    dataset_name = models.CharField(
        max_length=128,
        db_index=True,
    )

    group_type = models.CharField(
        max_length=20,
        db_index=True,
    )

    group_by = models.CharField(
        max_length=128,
        blank=True,
        default="",
        db_index=True,
    )

    axis_id = models.CharField(
        max_length=512,
        blank=True,
        default="",
    )

    axis_type = models.CharField(
        max_length=128,
        blank=True,
        default="",
    )

    axis_regulation = models.CharField(
        max_length=64,
        blank=True,
        default="",
    )

    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "axis_signature_project_index"
        constraints = [
            models.UniqueConstraint(
                fields=["axis_signature", "project"],
                name="uniq_axis_signature_project_index",
            ),
        ]
        indexes = [
            models.Index(
                fields=["axis_signature"],
                name="idx_axis_sig_project_sig",
            ),
            models.Index(
                fields=["axis_signature", "source"],
                name="idx_axis_sig_project_source",
            ),
            models.Index(
                fields=["axis_signature", "dataset_name"],
                name="idx_axis_sig_project_dataset",
            ),
            models.Index(
                fields=["axis_signature", "group_type"],
                name="idx_axis_sig_project_group",
            ),
            models.Index(
                fields=["dataset_name", "group_type"],
                name="idx_axis_sig_dataset_group",
            ),
        ]

    def __str__(self):
        return (
            f"{self.axis_signature} -> "
            f"{self.source}:{self.module}:"
            f"{self.dataset_name}:{self.group_type}:{self.group_by}"
        )


class AxisRecurrentSummary(models.Model):
    """
    Precomputed recurrent summary for one structural ceRNA axis.

    Rebuildable from AxisSignatureProjectIndex.
    """

    axis_signature = models.CharField(
        max_length=1024,
        unique=True,
        db_index=True,
    )

    axis_type = models.CharField(
        max_length=128,
        blank=True,
        default="",
        db_index=True,
    )

    miRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    mRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    lncRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    circRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    project_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )

    dataset_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )

    tcga_project_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )

    timedb_project_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )

    regulation_pattern_count = models.PositiveIntegerField(
        default=0,
        help_text=(
            "Number of distinct non-empty axis_regulation patterns."
        ),
    )

    dominant_axis_regulation = models.CharField(
        max_length=64,
        blank=True,
        default="",
        db_index=True,
    )

    dominant_regulation_count = models.PositiveIntegerField(
        default=0,
    )

    regulation_consistent = models.BooleanField(
        default=False,
        db_index=True,
        help_text=(
            "True when all non-empty occurrences have the same "
            "axis_regulation pattern."
        ),
    )

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "axis_recurrent_summary"
        indexes = [
            models.Index(
                fields=["-project_count", "axis_signature"],
                name="idx_axis_recur_project_count",
            ),
            models.Index(
                fields=["-dataset_count", "axis_signature"],
                name="idx_axis_recur_dataset_count",
            ),
            models.Index(
                fields=["axis_type", "-project_count"],
                name="idx_axis_recur_type_count",
            ),
            models.Index(
                fields=["regulation_consistent", "-project_count"],
                name="idx_axis_recur_reg_count",
            ),
            models.Index(
                fields=["miRNA", "mRNA"],
                name="idx_axis_recur_mirna_mrna",
            ),
        ]

    def __str__(self):
        return (
            f"{self.axis_signature}: "
            f"{self.project_count} projects"
        )


class AxisDatasetSource(models.TextChoices):
    TCGA = "TCGA", "TCGA"
    TIMEDB = "TIMEDB", "TIMEDB"


class AxisModule(models.TextChoices):
    MODULE2 = "module2", "Module 2"
    MODULE3 = "module3", "Module 3"


class AxisGroupType(models.TextChoices):
    NONE = "none", "None"
    OTHER = "other", "Other"
    GRADE = "grade", "Grade"
    STAGE = "stage", "Stage"


class AxisResultKind(models.TextChoices):
    AXIS_FINAL = "axis_final", "Axis Final"
    SPONGE = "sponge", "Sponge"


class CanonicalAxisType(models.TextChoices):
    MRNA_MIRNA = "mRNA-miRNA", "mRNA-miRNA"
    MRNA_MIRNA_LNCRNA = (
        "mRNA-miRNA-lncRNA",
        "mRNA-miRNA-lncRNA",
    )
    MRNA_MIRNA_CIRCRNA = (
        "mRNA-miRNA-circRNA",
        "mRNA-miRNA-circRNA",
    )


class AxisDatasetContext(models.Model):
    """
    One dataset/grouping context that may contain multiple Axis result types.

    Module 2:
        dataset_name must reference a DatasetMetadata.dataset value such as:
            TCGA_ACC_mRNA

    Module 3:
        dataset_name references a DatasetMetadata.dataset value such as:
            GSE20194
    """

    dataset_source = models.CharField(
        max_length=20,
        choices=AxisDatasetSource.choices,
        db_index=True,
    )

    module = models.CharField(
        max_length=20,
        choices=AxisModule.choices,
        db_index=True,
    )

    # The database column remains dataset_name, while the Django field is a
    # real foreign key to DatasetMetadata.dataset.
    dataset_metadata = models.ForeignKey(
        "DatasetMetadata",
        to_field="dataset",
        db_column="dataset_name",
        on_delete=models.PROTECT,
        related_name="axis_dataset_contexts",
    )

    group_type = models.CharField(
        max_length=20,
        choices=AxisGroupType.choices,
        default=AxisGroupType.NONE,
        db_index=True,
    )

    group_by = models.CharField(
        max_length=128,
        blank=True,
        default="",
        db_index=True,
    )

    annotation_dir_name = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    annotation_file_prefix = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    is_active = models.BooleanField(
        default=True,
        db_index=True,
    )

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "axis_dataset_context"

        constraints = [
            models.UniqueConstraint(
                fields=[
                    "dataset_source",
                    "module",
                    "dataset_metadata",
                    "group_type",
                    "group_by",
                ],
                name="uniq_axis_dataset_context",
            ),
        ]

        indexes = [
            models.Index(
                fields=[
                    "dataset_source",
                    "module",
                    "dataset_metadata",
                ],
                name="idx_axis_ctx_source_mod",
            ),
            models.Index(
                fields=[
                    "dataset_metadata",
                    "group_type",
                ],
                name="idx_axis_ctx_dataset_grp",
            ),
        ]

    @property
    def dataset_name(self) -> str:
        """
        Return the referenced DatasetMetadata.dataset value.

        Since dataset_metadata uses to_field='dataset', the foreign-key ID is
        the dataset name itself rather than DatasetMetadata's integer PK.
        """
        return str(self.dataset_metadata_id or "")

    def clean(self):
        super().clean()

        errors = {}
        dataset_name = self.dataset_name

        if self.module == AxisModule.MODULE2:
            if not dataset_name.endswith("_mRNA"):
                errors["dataset_metadata"] = (
                    "Module 2 context must reference a "
                    "DatasetMetadata.dataset ending with '_mRNA'."
                )

            if self.group_type != AxisGroupType.NONE:
                errors["group_type"] = (
                    "Module 2 context must use group_type='none'."
                )

            if self.group_by:
                errors["group_by"] = (
                    "Module 2 context must use an empty group_by value."
                )

        if errors:
            raise ValidationError(errors)

    def save(self, *args, **kwargs):
        # Context rows are low-volume, so ordinary writes enforce model-level
        # validation. QuerySet.update() and bulk_create() still bypass this.
        self.full_clean()

        return super().save(*args, **kwargs)

    def __str__(self):
        base = (
            f"{self.dataset_source}:"
            f"{self.module}:"
            f"{self.dataset_name}"
        )

        if self.group_type != AxisGroupType.NONE:
            return (
                f"{base}:"
                f"{self.group_type}:"
                f"{self.group_by}"
            )

        return base


class AxisResultArtifact(models.Model):
    """
    One imported Axis result file under an AxisDatasetContext.

    The same context may contain multiple result artifacts, such as:
        - Axis Final
        - Sponge
    """

    context = models.ForeignKey(
        AxisDatasetContext,
        on_delete=models.CASCADE,
        related_name="result_artifacts",
    )

    result_kind = models.CharField(
        max_length=32,
        choices=AxisResultKind.choices,
        db_index=True,
    )

    file_name = models.CharField(
        max_length=255,
    )

    file_path = models.TextField()

    file_sha256 = models.CharField(
        max_length=64,
        db_index=True,
    )

    row_count = models.PositiveIntegerField(
        default=0,
    )

    schema_version = models.CharField(
        max_length=32,
        blank=True,
        default="v1",
    )

    is_active = models.BooleanField(
        default=True,
        db_index=True,
    )

    imported_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "axis_result_artifact"

        constraints = [
            models.UniqueConstraint(
                fields=[
                    "context",
                    "result_kind",
                ],
                condition=Q(is_active=True),
                name="uniq_active_axis_artifact",
            ),
        ]

        indexes = [
            models.Index(
                fields=[
                    "context",
                    "result_kind",
                    "is_active",
                ],
                name="idx_axis_art_ctx_kind",
            ),
        ]

    def __str__(self):
        return (
            f"{self.context_id}:"
            f"{self.result_kind}:"
            f"{self.file_name}"
        )


class CanonicalAxis(models.Model):
    """
    Source-independent structural identity of one ceRNA Axis.

    axis_signature format:
        miRNA|mRNA|lncRNA|circRNA

    Examples:
        hsa-miR-210|ESPL1|CDCA3|
        hsa-miR-210|ESPL1||hsa_circ_000001
    """

    # Records the algorithm/schema version without embedding the version into
    # axis_signature.
    signature_version = models.PositiveSmallIntegerField(
        default=2,
    )

    axis_key = models.CharField(
        max_length=64,
        unique=True,
        help_text=(
            "SHA-256 of the normalized structural signature: "
            "miRNA|mRNA|lncRNA|circRNA."
        ),
    )

    axis_signature = models.TextField(
        help_text="miRNA|mRNA|lncRNA|circRNA",
    )

    axis_type = models.CharField(
        max_length=64,
        choices=CanonicalAxisType.choices,
        db_index=True,
    )

    miRNA = models.CharField(
        max_length=255,
        db_index=True,
    )

    mRNA = models.CharField(
        max_length=255,
        db_index=True,
    )

    lncRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    circRNA = models.CharField(
        max_length=255,
        blank=True,
        default="",
        db_index=True,
    )

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "canonical_axis"

        constraints = [
            models.CheckConstraint(
                condition=(
                    Q(miRNA__gt="")
                    & Q(mRNA__gt="")
                ),
                name="axis_need_mirna_mrna",
            ),
            models.CheckConstraint(
                condition=(
                    Q(lncRNA="")
                    | Q(circRNA="")
                ),
                name="axis_one_cerna_kind",
            ),
        ]

        indexes = [
            models.Index(
                fields=[
                    "miRNA",
                    "mRNA",
                ],
                name="idx_canon_axis_mir_mrna",
            ),
            models.Index(
                fields=[
                    "miRNA",
                    "mRNA",
                    "lncRNA",
                ],
                name="idx_canon_axis_lncrna",
            ),
            models.Index(
                fields=[
                    "miRNA",
                    "mRNA",
                    "circRNA",
                ],
                name="idx_canon_axis_circrna",
            ),
        ]

    def __str__(self):
        return self.axis_signature


class AxisObservation(models.Model):
    """
    One CanonicalAxis row observed in one imported result artifact.

    Structural identity is stored in CanonicalAxis. This model only records
    where the Axis was observed and its source-file row information.
    """

    artifact = models.ForeignKey(
        AxisResultArtifact,
        on_delete=models.CASCADE,
        related_name="observations",
    )

    axis = models.ForeignKey(
        CanonicalAxis,
        on_delete=models.PROTECT,
        related_name="observations",
    )

    row_index = models.PositiveIntegerField()

    source_axis_id = models.CharField(
        max_length=1024,
        blank=True,
        default="",
        db_index=True,
    )

    source_axis_type = models.CharField(
        max_length=128,
        blank=True,
        default="",
    )

    # Store only source columns that do not yet have stable typed fields.
    extra_data = models.JSONField(
        default=dict,
        blank=True,
    )

    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "axis_observation"

        constraints = [
            models.UniqueConstraint(
                fields=[
                    "artifact",
                    "row_index",
                ],
                name="uniq_axis_artifact_row",
            ),
            models.UniqueConstraint(
                fields=[
                    "artifact",
                    "axis",
                ],
                name="uniq_axis_artifact_axis",
            ),
        ]

        indexes = [
            models.Index(
                fields=[
                    "axis",
                    "artifact",
                ],
                name="idx_axis_obs_axis_art",
            ),
        ]

    @property
    def result_kind(self) -> str:
        return self.artifact.result_kind

    def __str__(self):
        return f"{self.artifact_id}:{self.axis_id}"


class AxisFinalEvidence(models.Model):
    """
    Axis Final-specific evidence associated with one AxisObservation.
    """

    observation = models.OneToOneField(
        AxisObservation,
        on_delete=models.CASCADE,
        primary_key=True,
        related_name="axis_final_evidence",
    )

    axis_regulation = models.CharField(
        max_length=64,
        blank=True,
        default="",
        db_index=True,
    )

    mRNA_log2FC = models.FloatField(
        null=True,
        blank=True,
    )

    mRNA_regulation = models.CharField(
        max_length=32,
        blank=True,
        default="",
        db_index=True,
    )

    miRNA_log2FC = models.FloatField(
        null=True,
        blank=True,
    )

    miRNA_regulation = models.CharField(
        max_length=32,
        blank=True,
        default="",
        db_index=True,
    )

    lncRNA_log2FC = models.FloatField(
        null=True,
        blank=True,
    )

    lncRNA_regulation = models.CharField(
        max_length=32,
        blank=True,
        default="",
        db_index=True,
    )

    circRNA_log2FC = models.FloatField(
        null=True,
        blank=True,
    )

    circRNA_regulation = models.CharField(
        max_length=32,
        blank=True,
        default="",
        db_index=True,
    )

    class Meta:
        db_table = "axis_final_evidence"

    def clean(self):
        super().clean()

        if (
            self.observation_id
            and self.observation.artifact.result_kind
            != AxisResultKind.AXIS_FINAL
        ):
            raise ValidationError({
                "observation": (
                    "AxisFinalEvidence requires an "
                    "axis_final result artifact."
                ),
            })

    def __str__(self):
        return f"axis_final:{self.observation_id}"


class SpongeEvidence(models.Model):
    """
    Sponge-specific statistical evidence associated with one AxisObservation.
    """

    observation = models.OneToOneField(
        AxisObservation,
        on_delete=models.CASCADE,
        primary_key=True,
        related_name="sponge_evidence",
    )

    cor = models.FloatField(
        null=True,
        blank=True,
    )

    pcor = models.FloatField(
        null=True,
        blank=True,
    )

    mscor = models.FloatField(
        null=True,
        blank=True,
    )

    class Meta:
        db_table = "sponge_evidence"

        indexes = [
            models.Index(
                fields=["cor"],
                name="idx_sponge_cor",
            ),
            models.Index(
                fields=["pcor"],
                name="idx_sponge_pcor",
            ),
            models.Index(
                fields=["mscor"],
                name="idx_sponge_mscor",
            ),
        ]

    def clean(self):
        super().clean()

        if (
            self.observation_id
            and self.observation.artifact.result_kind
            != AxisResultKind.SPONGE
        ):
            raise ValidationError({
                "observation": (
                    "SpongeEvidence requires a "
                    "sponge result artifact."
                ),
            })

    def __str__(self):
        return f"sponge:{self.observation_id}"


class AxisContextPresence(models.Model):
    """
    Rebuildable context-level index for matching and recurrence queries.

    One row means that one CanonicalAxis is present in one AxisDatasetContext,
    regardless of how many result artifacts contain the Axis.
    """

    context = models.ForeignKey(
        AxisDatasetContext,
        on_delete=models.CASCADE,
        related_name="axis_presences",
    )

    axis = models.ForeignKey(
        CanonicalAxis,
        on_delete=models.CASCADE,
        related_name="context_presences",
    )

    observation_count = models.PositiveIntegerField(
        default=0,
    )

    axis_final_observation_count = models.PositiveIntegerField(
        default=0,
    )

    sponge_observation_count = models.PositiveIntegerField(
        default=0,
    )

    has_axis_final = models.BooleanField(
        default=False,
        db_index=True,
    )

    has_sponge = models.BooleanField(
        default=False,
        db_index=True,
    )

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "axis_context_presence"

        constraints = [
            models.UniqueConstraint(
                fields=[
                    "context",
                    "axis",
                ],
                name="uniq_axis_context_presence",
            ),
        ]

        indexes = [
            models.Index(
                fields=[
                    "axis",
                    "context",
                ],
                name="idx_axis_presence_axis_ctx",
            ),
            models.Index(
                fields=[
                    "axis",
                    "has_axis_final",
                ],
                name="idx_axis_presence_final",
            ),
            models.Index(
                fields=[
                    "axis",
                    "has_sponge",
                ],
                name="idx_axis_presence_sponge",
            ),
        ]

    def __str__(self):
        return f"{self.axis_id}@{self.context_id}"


class AxisStructureRecurrentSummary(models.Model):
    axis = models.OneToOneField(
        CanonicalAxis,
        on_delete=models.CASCADE,
        primary_key=True,
        related_name="recurrent_summary",
    )

    # Primary recurrence metrics
    context_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )
    dataset_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )

    # Distinct datasets by source
    tcga_dataset_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )
    timedb_dataset_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )

    # Context counts by source
    tcga_context_count = models.PositiveIntegerField(default=0)
    timedb_context_count = models.PositiveIntegerField(default=0)

    # Context counts by module
    module2_context_count = models.PositiveIntegerField(default=0)
    module3_context_count = models.PositiveIntegerField(default=0)

    # Result coverage
    axis_final_context_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )
    sponge_context_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )
    both_result_context_count = models.PositiveIntegerField(
        default=0,
    )

    summary_version = models.PositiveSmallIntegerField(default=1)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "axis_structure_recur_summary"

        indexes = [
            models.Index(
                fields=[
                    "-context_count",
                    "axis",
                ],
                name="idx_axis_recur_ctx_count",
            ),
            models.Index(
                fields=[
                    "-dataset_count",
                    "axis",
                ],
                name="idx_axis_recur_ds_count",
            ),
        ]

    def __str__(self):
        return (
            f"{self.axis_id}:"
            f"{self.context_count} contexts"
        )


class AxisFinalRecurrentSummary(models.Model):
    """
    Rebuildable Axis Final-specific regulation summary.
    """

    axis = models.OneToOneField(
        CanonicalAxis,
        on_delete=models.CASCADE,
        primary_key=True,
        related_name="axis_final_recurrent_summary",
    )

    context_count = models.PositiveIntegerField(
        default=0,
        db_index=True,
    )

    observation_count = models.PositiveIntegerField(
        default=0,
    )

    regulation_pattern_count = models.PositiveIntegerField(
        default=0,
    )

    dominant_axis_regulation = models.CharField(
        max_length=64,
        blank=True,
        default="",
        db_index=True,
    )

    dominant_regulation_count = models.PositiveIntegerField(
        default=0,
    )

    regulation_consistent = models.BooleanField(
        default=False,
        db_index=True,
    )

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "axis_final_recur_summary"

        indexes = [
            models.Index(
                fields=[
                    "regulation_consistent",
                    "-context_count",
                ],
                name="idx_axis_final_recur_reg",
            ),
        ]

    def __str__(self):
        return (
            f"{self.axis_id}:"
            f"{self.context_count} contexts"
        )
