import gzip
import csv
from pathlib import Path

bed_path = Path(r"E:\Projects\ceRNAxis\data\circRNA_interactions\circRNA_miRNA_interaction_sorted.bed.gz")
csv_path = Path(r"E:\Projects\ceRNAxis\data\circRNA_interactions\circRNA_miRNA_interaction.csv")

with gzip.open(bed_path, "rt", encoding="utf-8") as fin, open(
    csv_path,
    "w",
    newline="",
    encoding="utf-8-sig",
) as fout:
    writer = csv.DictWriter(
        fout,
        fieldnames=[
            "miRNA",
            "ceRNA",
            "species",
            "database",
            "type",
        ],
    )

    writer.writeheader()

    total_lines = 0
    written_rows = 0
    skipped_rows = 0

    for line in fin:
        total_lines += 1
        line = line.strip()

        if not line or line.startswith("#"):
            skipped_rows += 1
            continue

        fields = line.split("\t")

        if len(fields) < 5:
            skipped_rows += 1
            continue

        circRNA_location = fields[3]
        miRNA_name = fields[4]

        writer.writerow({
            "miRNA": miRNA_name,
            "ceRNA": circRNA_location,
            "species": "Homo sapiens",
            "database": "circNet2.0",
            "type": "miRNA-circRNA",
        })

        written_rows += 1

print(f"Input file: {bed_path}")
print(f"Output file: {csv_path}")
print(f"Total lines: {total_lines}")
print(f"Written rows: {written_rows}")
print(f"Skipped rows: {skipped_rows}")
