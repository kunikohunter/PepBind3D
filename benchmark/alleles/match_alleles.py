#!/usr/bin/env python3
"""
Assign a specific MHC allele to each of 223 benchmark crystal structures by
SEQUENCE (not by parsing the free-text allele annotation), and characterise
the ones whose text is unresolvable.

Outputs (written under /data/p_csb_meiler/huntek1/benchmark/alleles/):
  resolved_alleles.csv   - one row per pdb_id with the sequence-based call
  match_report.txt       - human-readable summary / category breakdown
"""
import csv, re, os, sys
from collections import Counter, defaultdict

REFS = "/data/p_csb_meiler/huntek1/benchmark/refs"
MANIFEST = os.path.join(REFS, "reference_manifest.csv")
CIF_DIR = os.path.join(REFS, "cif")
ALLELE_SEQ_INFO = "/home/huntek1/Data/MHC_database/build/allele_seq.info"
OUTDIR = "/data/p_csb_meiler/huntek1/benchmark/alleles"

# 3-letter -> 1-letter amino acid code (for parsing _entity_poly_seq as a fallback)
AA3to1 = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'
}

# Canonical alpha1/alpha2 domain references, taken verbatim from HLA_db.py's
# trim_seq() method (canon_domains dict) so our trimming matches its behaviour.
CANON_DOMAINS = {
 "A": "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQIMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQLRAYLDGTCVEWLRRYLENGKETLQRT",
 "B": "GSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQSMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGECVEWLRRYLENGKDKLERA",
 "C": "CSHSMRYFDTAVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRETQKYKRQAQADRVSLRNLRGYYNQSEDGSHTLQRMSGCDLGPDGRLLRGYDQSAYDGKDYIALNEDLRSWTAADTAAQITQRKLEAARAAEQLRAYLEGTCVEWLRRYLENGKETLQRA",
 "H2-Db": "GPHSMRYFETAVSRPGLEEPRYISVGYVDNKEFVRFDSDAENPRYEPRAPWMEQEGPEYWERETQKAKGQEQWFRVSLRNLLGYYNQSAGGSHTLQQMSGCDLGSDWRLLRGYLQFAYEGRDYIALNEDLKTWTAADMAAQITRRKWEQSGAAEHYKAYLEGECVEWLHRYLKNGNATLLRT",
 "H2-Kb": "GPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRT",
}

def sliding_trim(seq, canon_domain):
    """Reimplementation of HLA_db.py MHCdatabase.trim_seq()'s sliding-window
    core: find the offset of `seq` that maximises character-identity overlap
    against `canon_domain`, then return the substring of `seq` aligned under
    the canonical domain window, plus the match score."""
    len_diff = len(seq) - len(canon_domain)
    if len_diff <= 0:
        # Sequence is not longer than the canonical domain: HLA_db.py returns it unchanged.
        score = sum(1 for a, b in zip(seq, canon_domain) if a == b)
        return seq, score

    canon_win = list(canon_domain) + [0] * len(seq)
    seq_win = [0] * len(canon_domain) + list(seq)
    best_score = -1
    best_offset = 0
    n = len(canon_win)
    for ioffset in range(n):
        score = 0
        for i, aa in enumerate(seq_win):
            c = canon_win[i]
            if aa != 0 and aa == c:
                score += 1
        if score > best_score:
            best_score = score
            best_offset = ioffset
        # roll canon_win right by one position (equivalent to python list pop/insert in orig)
        canon_win.pop()
        canon_win = [0] + canon_win

    # Reconstruct the winning alignment directly (avoids storing all alignments, which is O(n^2) memory)
    canon_win2 = list(canon_domain) + [0] * len(seq)
    for _ in range(best_offset):
        canon_win2.pop()
        canon_win2 = [0] + canon_win2
    seq_win2 = [0] * len(canon_domain) + list(seq)
    newseq = "".join(seq_win2[i] for i in range(len(canon_win2))
                      if canon_win2[i] != 0 and seq_win2[i] != 0)
    return newseq, best_score

def best_trim(seq):
    """Try all canonical domains (human A/B/C + murine H2-Db/H2-Kb) and return
    the trimmed alpha1/alpha2 sequence and locus-guess for whichever canon
    reference gives the highest raw identity score. This generalises
    HLA_db.py's per-locus trim_seq() to loci not already known from text."""
    results = {}
    for name, canon in CANON_DOMAINS.items():
        trimmed, score = sliding_trim(seq, canon)
        results[name] = (trimmed, score, len(canon))
    best_name = max(results, key=lambda k: results[k][1] / results[k][2])
    trimmed, score, canon_len = results[best_name]
    return best_name, trimmed, score / canon_len

def parse_cif_entities(path):
    """Parse _entity_poly loop from an mmCIF file. Returns dict:
    strand_id (individual chain letter) -> one-letter sequence (canonical, no newlines)."""
    text = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"_entity_poly\.pdbx_target_identifier\s*\n(.*?)\n#\s*\n", text, re.S)
    if not m:
        return {}
    block = m.group(1)
    # Records are either single-line 'id type nstd nstd seq seq_can strand target'
    # or multi-line with ';'-delimited seq fields. Parse robustly by scanning tokens
    # and ';'-quoted blocks in order.
    chain_to_seq = {}
    # Tokenize preserving ;...; multi-line quoted blocks
    tokens = []
    lines = block.split("\n")
    i = 0
    cur = []
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            # multi-line quoted value; gather until a line that is exactly ';'
            val_lines = [line[1:]]
            i += 1
            while i < len(lines) and lines[i].strip() != ";":
                val_lines.append(lines[i])
                i += 1
            tokens.append("".join(val_lines))
            i += 1
            continue
        else:
            # split respecting simple space separation; strand ids w/ commas have no spaces
            parts = line.split()
            tokens.extend(parts)
            i += 1
    # Now tokens is a flat stream; each entity record has 8 fields:
    # entity_id type nstd_linkage nstd_monomer seq seq_can strand_id target_identifier
    # 'polypeptide(L)' is a single quoted token already split correctly since no internal spaces except that phrase itself is one token without spaces. Good.
    idx = 0
    while idx + 7 < len(tokens) + 1 and idx < len(tokens):
        # entity_id must look like an integer
        if not re.match(r"^\d+$", tokens[idx]):
            idx += 1
            continue
        try:
            entity_id, etype, nstd1, nstd2, seq, seq_can, strand, target = tokens[idx:idx+8]
        except ValueError:
            break
        if etype.strip("'\"").startswith("polypeptide"):
            seq_clean = seq_can.replace("\n", "").strip()
            for ch in strand.split(","):
                chain_to_seq[ch] = seq_clean
        idx += 8
    return chain_to_seq

def load_allele_seq_info(path):
    """Return list of (sequence, [allele names])."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append((r["Sequence"], r["Allele_Names"].split()))
    return rows

def identity(a, b):
    """Percent identity between two equal-or-near-length strings, aligned at
    position 0 (both already trimmed to alpha1/alpha2 domain length ~180)."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    matches = sum(1 for x, y in zip(a, b) if x == y)
    return 100.0 * matches / n

MAX_SHIFT = 3

def identity_best_shift(query, ref, max_shift=MAX_SHIFT):
    """%identity between query and ref, allowing query's N-terminus to be
    offset by up to max_shift residues relative to ref. Needed because some
    deposited chains are missing 1-2 N-terminal residues that our canon-domain
    sliding trim does not itself correct for (e.g. a chain modeled starting at
    'HSMRYF...' instead of 'SHSMRYF...' loses an extra residue relative to the
    reference alpha1/alpha2 domain, and the two sequences would otherwise be
    frame-shifted by one position for their entire length)."""
    best = identity(query, ref)
    for shift in range(1, max_shift + 1):
        best = max(best, identity(query[shift:], ref), identity(query, ref[shift:]))
    return best

def best_allele_match(query180, allele_rows):
    best_pct = -1.0
    best_alleles = []
    for seq, names in allele_rows:
        pct = identity_best_shift(query180, seq)
        if pct > best_pct:
            best_pct = pct
            best_alleles = names
        elif pct == best_pct:
            best_alleles = best_alleles + names
    return best_pct, best_alleles

def classify_unresolved(locus, allele_text, best_name, best_pct):
    """Assign a category to a target that did not resolve to a confident
    classical HLA-A/B/C call."""
    txt = allele_text.lower()
    if locus == "murine" or best_name in ("H2-Db", "H2-Kb"):
        return "non-human (murine H-2)"
    if "rat" in txt or "rt1" in txt:
        return "non-human (rat RT1)"
    if re.search(r"hla-?e\b", txt) or re.search(r"hla-?f\b", txt) or re.search(r"hla-?g\b", txt) \
       or "mr1" in txt or "cd1" in txt or locus == "other-nonclassical":
        return "non-classical human (HLA-E/F/G, MR1, CD1)"
    if "chimeric" in txt or "single-chain" in txt or "single chain" in txt or "sctcr" in txt \
       or "fusion" in txt or "tag" in txt or "trimer" in txt or "construct" in txt:
        return "engineered/chimeric/fusion construct"
    if locus in ("HLA-A", "HLA-B", "HLA-C") or re.search(r"hla-[abc]", txt):
        # Depositor's own metadata already calls this a classical A/B/C locus;
        # only the subtype digits are missing from the free text.
        return "classical HLA, depositor gave no subtype"
    if "mhc class i" in txt:
        # Depositor text is generic ("MHC class I antigen") and locus is not
        # even pinned to A/B/C by the manifest (locus == unclassified). Use
        # the sequence identity itself to decide: a decent identity to the
        # nearest classical allele (even short of the 90% "usable" bar) means
        # this is very likely a classical HLA that is just weakly resolved by
        # sequence (e.g. a divergent construct or partial chain); a very low
        # identity means this heavy chain does not really resemble classical
        # HLA-A/B/C at all and is better described as an unidentified/likely
        # non-classical or highly divergent MHC-I molecule.
        if best_pct >= 60:
            return "classical HLA, depositor gave no subtype (weak sequence match)"
        return "unidentified: text says 'MHC class I antigen' but sequence does not resemble classical HLA-A/B/C (likely non-classical or highly divergent)"
    return "other/unclear"

def main():
    manifest_rows = list(csv.DictReader(open(MANIFEST)))
    allele_rows = load_allele_seq_info(ALLELE_SEQ_INFO)
    print(f"Loaded {len(manifest_rows)} manifest rows, {len(allele_rows)} allele_seq.info entries", file=sys.stderr)

    out_rows = []
    n_exact = 0
    n_nearest = defaultdict(int)  # bucket by identity band
    n_unresolved = 0
    category_counts = Counter()
    category_examples = defaultdict(list)
    extraction_fail = []

    EXACT_THRESHOLD = 100.0
    NEAREST_THRESHOLD = 90.0  # below this we call it "no confident match" -> unresolved by sequence too

    for row in manifest_rows:
        pdb_id = row["pdb_id"]
        locus = row["locus"]
        allele_text = row["allele_text"]
        heavy_chain_id = row["heavy_chain_id"]
        cif_path = os.path.join(CIF_DIR, f"{pdb_id}.cif")

        out = {
            "pdb_id": pdb_id,
            "allele_text_original": allele_text,
            "resolved_allele": "",
            "percent_identity": "",
            "match_type": "unresolved",
            "category": "",
            "usable_for_rosetta": "no",
        }

        if not os.path.exists(cif_path):
            out["category"] = "cif file missing"
            extraction_fail.append((pdb_id, "cif missing"))
            out_rows.append(out)
            continue

        chain_to_seq = parse_cif_entities(cif_path)
        seq = chain_to_seq.get(heavy_chain_id)
        if seq is None:
            extraction_fail.append((pdb_id, f"heavy_chain_id {heavy_chain_id} not found among entities {list(chain_to_seq.keys())}"))
            out["category"] = "heavy chain sequence extraction failed"
            out_rows.append(out)
            continue

        if not (250 <= len(seq) <= 400):
            extraction_fail.append((pdb_id, f"heavy chain length {len(seq)} outside expected 250-400"))
        if not re.search(r"SHSMRYF|SLRYF|SHSLRYF", seq[:20]):
            extraction_fail.append((pdb_id, f"heavy chain does not start SHSMRYF-like: {seq[:20]}"))

        best_canon_name, trimmed, frac = best_trim(seq)

        # Chain rescue: several manifest rows flag "multiple candidate heavy-chain
        # entities" and the recorded heavy_chain_id sometimes turns out to be a
        # TCR alpha/beta chain instead of the MHC heavy chain (confirmed by
        # inspection: these give a near-zero identity fraction against every
        # canonical alpha1/alpha2 domain, ~0.09-0.10, vs ~0.7-1.0 for genuine
        # heavy chains). When that happens, search the file's other polypeptide
        # entities (length 200-450) for one that aligns far better, and use it
        # instead. This is purely a chain-identification fix based on the
        # structure itself; it does not touch the manifest or any files outside
        # the alleles/ working directory.
        rescue_note = ""
        if frac < 0.35:
            candidates = []
            for ch, cseq in chain_to_seq.items():
                if ch == heavy_chain_id or not (200 <= len(cseq) <= 450):
                    continue
                cname, ctrimmed, cfrac = best_trim(cseq)
                candidates.append((cfrac, ch, cname, ctrimmed, cseq))
            if candidates:
                candidates.sort(reverse=True)
                cfrac, ch, cname, ctrimmed, cseq = candidates[0]
                if cfrac > 0.5 and cfrac > frac:
                    rescue_note = (f"heavy_chain_id in manifest ({heavy_chain_id}) looks like a "
                                    f"non-MHC chain (alpha1/alpha2 alignment fraction {frac:.2f}); "
                                    f"reassigned to chain {ch} (alignment fraction {cfrac:.2f})")
                    extraction_fail.append((pdb_id, rescue_note))
                    seq, best_canon_name, trimmed, frac = cseq, cname, ctrimmed, cfrac
        # allele_seq.info reference sequences omit the very first residue (G/C) of
        # the canon_domains strings used above (181 vs 182 aa). Whether `trimmed`
        # itself still carries that leading residue depends on whether the
        # deposited construct included it (some heavy chains in the CIF already
        # start at "SHSMRYF..." with no leading G/C, e.g. 24ST chain C, so the
        # sliding-window trim returns a 181-aa result already flush with
        # allele_seq.info). Rather than assume either way, try both alignments
        # and keep whichever scores higher.
        if best_canon_name in ("H2-Db", "H2-Kb"):
            # Murine: do not bother matching against the human allele_seq.info
            best_pct, best_alleles = 0.0, []
        else:
            pct_a, alleles_a = best_allele_match(trimmed, allele_rows)
            query180 = trimmed[1:] if len(trimmed) >= 181 else trimmed
            pct_b, alleles_b = best_allele_match(query180, allele_rows)
            if pct_a >= pct_b:
                best_pct, best_alleles = pct_a, alleles_a
            else:
                best_pct, best_alleles = pct_b, alleles_b

        if best_pct >= EXACT_THRESHOLD:
            match_type = "exact"
            n_exact += 1
            resolved = best_alleles[0]
            category = "classical HLA (resolved by sequence)"
            usable = "yes"
        elif best_pct >= NEAREST_THRESHOLD:
            match_type = "nearest"
            band = f"{int(best_pct//5)*5}-{int(best_pct//5)*5+5}%"
            n_nearest[band] += 1
            resolved = best_alleles[0] if best_alleles else ""
            category = "classical HLA (resolved by sequence, nearest-neighbour)"
            usable = "yes"
        else:
            match_type = "unresolved"
            n_unresolved += 1
            resolved = ""
            category = classify_unresolved(locus, allele_text, best_canon_name, best_pct)
            category_counts[category] += 1
            if len(category_examples[category]) < 4:
                category_examples[category].append((pdb_id, allele_text, best_canon_name, round(best_pct, 1)))
            usable = "no"

        out.update({
            "resolved_allele": resolved,
            "percent_identity": f"{best_pct:.1f}",
            "match_type": match_type,
            "category": category,
            "usable_for_rosetta": usable,
        })
        out_rows.append(out)

    # Write CSV
    with open(os.path.join(OUTDIR, "resolved_alleles.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id", "allele_text_original", "resolved_allele",
                                           "percent_identity", "match_type", "category", "usable_for_rosetta"])
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    # Report
    with open(os.path.join(OUTDIR, "match_report.txt"), "w") as f:
        def p(*a):
            print(*a)
            print(*a, file=f)
        p(f"Total targets: {len(manifest_rows)}")
        p(f"Exact (100%) sequence match to a classical HLA-A/B/C allele: {n_exact}")
        nearest_total = sum(n_nearest.values())
        p(f"Nearest-neighbour match (90-99.9% identity) to classical HLA-A/B/C: {nearest_total}")
        for band in sorted(n_nearest):
            p(f"   {band} identity: {n_nearest[band]}")
        p(f"Unresolved (no >=90% match to any classical HLA-A/B/C allele): {n_unresolved}")
        p("")
        p("Category breakdown of the unresolved set:")
        for cat, cnt in category_counts.most_common():
            p(f"  {cnt:3d}  {cat}")
            for ex in category_examples[cat]:
                p(f"        e.g. {ex[0]}: allele_text='{ex[1]}' best_canon={ex[2]} best_pct_vs_HLA={ex[3]}")
        p("")
        if extraction_fail:
            p(f"Extraction anomalies flagged ({len(extraction_fail)}):")
            for pdb_id, msg in extraction_fail[:40]:
                p(f"  {pdb_id}: {msg}")
            if len(extraction_fail) > 40:
                p(f"  ... and {len(extraction_fail)-40} more")

if __name__ == "__main__":
    main()
