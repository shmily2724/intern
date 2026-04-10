cat > fasta_to_arp.py <<'PY'
#!/usr/bin/env python3
import sys
import csv
from collections import OrderedDict, defaultdict

if len(sys.argv) != 4:
    print("Usage: python3 fasta_to_arp.py samples.fasta metadata.tsv output.arp")
    sys.exit(1)

fasta_file = sys.argv[1]
meta_file = sys.argv[2]
out_file = sys.argv[3]

seqs = OrderedDict()
current_id = None
current_seq = []

with open(fasta_file, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                seqs[current_id] = "".join(current_seq).upper()
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)

if current_id is not None:
    seqs[current_id] = "".join(current_seq).upper()

if not seqs:
    raise ValueError("No sequences found in FASTA.")

lengths = {len(s) for s in seqs.values()}
if len(lengths) != 1:
    raise ValueError(f"Sequences are not aligned: found multiple lengths {sorted(lengths)}")

cleaned = OrderedDict()
for sid, seq in seqs.items():
    seq = seq.replace("N", "?")
    bad = set(seq) - set("ACGT?-")
    if bad:
        raise ValueError(f"Sequence {sid} contains unsupported characters: {sorted(bad)}")
    cleaned[sid] = seq
seqs = cleaned

id_to_region = {}
with open(meta_file, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    required = {"ID", "Region"}
    if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
        raise ValueError("metadata.tsv must contain tab-separated columns: ID and Region")

    for row in reader:
        sid = row["ID"].strip()
        region = row["Region"].strip()
        if sid and region:
            id_to_region[sid] = region

missing_meta = [sid for sid in seqs if sid not in id_to_region]
if missing_meta:
    raise ValueError(
        "These FASTA IDs are missing in metadata.tsv:\n" + "\n".join(missing_meta[:20])
    )

groups = defaultdict(list)
for sid, seq in seqs.items():
    region = id_to_region[sid]
    groups[region].append((sid, seq))

region_order = ["North", "South", "Central"]
remaining = [r for r in groups.keys() if r not in region_order]
region_order = [r for r in region_order if r in groups] + sorted(remaining)

with open(out_file, "w", encoding="utf-8", newline="\n") as out:
    out.write('[Profile]\n')
    out.write('Title="Kinh mtDNA"\n')
    out.write(f'NbSamples={len(region_order)}\n')
    out.write('GenotypicData=0\n')
    out.write('DataType=DNA\n')
    out.write('LocusSeparator=NONE\n')
    out.write("MissingData='?'\n\n")

    out.write('[Data]\n')
    out.write('[[Samples]]\n\n')

    for region in region_order:
        samples = groups[region]
        out.write(f'SampleName="{region}"\n')
        out.write(f'SampleSize={len(samples)}\n')
        out.write('SampleData={\n')
        for sid, seq in samples:
            out.write(f'{sid} 1 {seq}\n')
        out.write('}\n\n')

    out.write('[[Structure]]\n')
    out.write('StructureName="Kinh_3pops"\n')
    out.write('NbGroups=1\n')
    out.write('Group={\n')
    for region in region_order:
        out.write(f'"{region}"\n')
    out.write('}\n')

print(f"Wrote {out_file}")
PY
