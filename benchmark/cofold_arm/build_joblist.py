import csv, json, sys

RES = '/data/p_csb_meiler/huntek1/benchmark/alleles/resolved_alleles.csv'
REF = '/data/p_csb_meiler/huntek1/benchmark/refs/reference_manifest.csv'
SEQ = '/home/huntek1/Data/MHC_database/build/allele_seq.info'

# load allele -> mhc_seq
allele_to_seq = {}
with open(SEQ) as f:
    r = csv.DictReader(f)
    for row in r:
        seq = row['Sequence'].strip()
        for name in row['Allele_Names'].split():
            allele_to_seq[name.strip()] = seq

# load reference manifest keyed by pdb_id
ref = {}
with open(REF) as f:
    for row in csv.DictReader(f):
        ref[row['pdb_id']] = row

rows = []
skipped = []
missing_seq = []
with open(RES) as f:
    for row in csv.DictReader(f):
        if row['usable_for_rosetta'] != 'yes':
            continue
        pdb_id = row['pdb_id']
        allele = row['resolved_allele']
        r = ref.get(pdb_id)
        if r is None:
            skipped.append((pdb_id, 'no reference_manifest entry'))
            continue
        pep = r['peptide_seq']
        if 'X' in pep:
            skipped.append((pdb_id, f'peptide has non-standard residue X ({pep})'))
            continue
        mhc_seq = allele_to_seq.get(allele)
        if mhc_seq is None:
            missing_seq.append((pdb_id, allele))
            continue
        rows.append({
            'pdb_id': pdb_id,
            'allele': allele,
            'peptide_seq': pep,
            'mhc_seq': mhc_seq,
            'heavy_chain_id': r['heavy_chain_id'],
            'b2m_chain_id': r['b2m_chain_id'],
            'peptide_chain_id': r['peptide_chain_id'],
        })

print(f"Total usable_for_rosetta==yes: 174", file=sys.stderr)
print(f"Skipped (peptide X or no ref): {len(skipped)}", file=sys.stderr)
for s in skipped:
    print(f"  SKIP {s[0]}: {s[1]}", file=sys.stderr)
print(f"Missing mhc_seq lookup: {len(missing_seq)}", file=sys.stderr)
for m in missing_seq:
    print(f"  NOSEQ {m[0]}: {m[1]}", file=sys.stderr)
print(f"Final job count: {len(rows)}", file=sys.stderr)

with open('/data/p_csb_meiler/huntek1/benchmark/cofold_arm/joblist.tsv', 'w') as out:
    out.write("idx\tpdb_id\tallele\tpeptide_seq\tmhc_seq\n")
    for i, r in enumerate(rows, start=1):
        out.write(f"{i}\t{r['pdb_id']}\t{r['allele']}\t{r['peptide_seq']}\t{r['mhc_seq']}\n")

with open('/data/p_csb_meiler/huntek1/benchmark/cofold_arm/joblist_full.json', 'w') as out:
    json.dump(rows, out, indent=2)
