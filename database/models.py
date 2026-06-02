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
