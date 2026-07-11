from django.db import models


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
