import gzip

bed_path = "E:\\Projects\\ceRNAxis\\data\\circRNA_interactions\\circRNA_miRNA_interaction_sorted.bed.gz"

total_lines = 0
data_lines = 0
empty_lines = 0
comment_lines = 0

with gzip.open(bed_path, "rt", encoding="utf-8") as f:
    for line in f:
        total_lines += 1

        line = line.rstrip("\n")

        if not line:
            empty_lines += 1
        elif line.startswith("#"):
            comment_lines += 1
        else:
            data_lines += 1

print(f"Total lines: {total_lines}")
print(f"Data lines: {data_lines}")
print(f"Comment lines: {comment_lines}")
print(f"Empty lines: {empty_lines}")
