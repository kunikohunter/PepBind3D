#!/usr/bin/env python3
"""
Build reference manifest for pMHC-I benchmark structures.
- Selects rows from expanded_candidates.csv (release_date > 2024-01-01, peptide 8-15 aa)
- Downloads mmCIF from RCSB into refs/cif/
- Uses RCSB Data API to classify entities (heavy chain / b2m / peptide) and count ASU copies
- Parses the downloaded CIF's atom_site records for the peptide chain to check backbone
  (N, CA, C) completeness and residue-numbering gaps
- Writes refs/reference_manifest.csv
"""
import csv, os, sys, time, json, re
from datetime import date
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests

BASE = "/data/p_csb_meiler/huntek1/benchmark/refs"
CIF_DIR = os.path.join(BASE, "cif")
SRC_CSV = "/data/p_csb_meiler/huntek1/altpred/rcsb_expand/expanded_candidates.csv"
MANIFEST = os.path.join(BASE, "reference_manifest.csv")
LOG = os.path.join(BASE, "build_log.txt")

os.makedirs(CIF_DIR, exist_ok=True)

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "benchmark-refs-build/1.0 (research; contact via institution)"})

def log(msg):
    print(msg)
    with open(LOG, "a") as f:
        f.write(msg + "\n")

def load_selection():
    rows = list(csv.DictReader(open(SRC_CSV)))
    sel = []
    for r in rows:
        try:
            d = date.fromisoformat(r["initial_release_date"])
        except Exception:
            continue
        try:
            plen = int(r["peptide_len"])
        except Exception:
            continue
        if d > date(2024, 1, 1) and 8 <= plen <= 15:
            sel.append(r)
    return sel

def api_get(url, retries=4, timeout=20):
    for i in range(retries):
        try:
            r = SESSION.get(url, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            elif r.status_code == 404:
                return None
            else:
                time.sleep(1.0 + i)
        except Exception:
            time.sleep(1.0 + i)
    return None

def download_cif(pdb_id):
    out = os.path.join(CIF_DIR, f"{pdb_id}.cif")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out, "skipped_existing"
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    for i in range(4):
        try:
            r = SESSION.get(url, timeout=30)
            if r.status_code == 200 and len(r.text) > 100:
                with open(out, "w") as f:
                    f.write(r.text)
                return out, "downloaded"
            else:
                time.sleep(1.5)
        except Exception:
            time.sleep(1.5)
    return None, "failed"

HEAVY_KEYWORDS = ["class i", "histocompatibility", "alpha chain", "hla-", "h-2", "mhc class i"]
B2M_KEYWORDS = ["microglobulin", "b2m", "beta-2", "beta 2"]

def classify_entities(entry_json, pdb_id):
    """Return dict entity_id -> role info, using polymer_entity API."""
    pe_ids = entry_json.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids") or []
    entities = {}
    for pid in pe_ids:
        pe = api_get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{pid}")
        if pe is None:
            continue
        desc = (pe.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description") or ""
        seq = (pe.get("entity_poly", {}) or {}).get("pdbx_seq_one_letter_code_can") or ""
        seq = seq.replace("\n", "")
        seqlen = len(seq)
        ids = pe.get("rcsb_polymer_entity_container_identifiers", {}) or {}
        auth_asym = ids.get("auth_asym_ids") or []
        dl = desc.lower()
        role = "other"
        if 8 <= seqlen <= 15:
            role = "peptide"
        elif 90 <= seqlen <= 135 and any(k in dl for k in B2M_KEYWORDS):
            role = "b2m"
        elif 150 <= seqlen <= 400 and any(k in dl for k in HEAVY_KEYWORDS):
            role = "heavy"
        elif 90 <= seqlen <= 135:
            role = "b2m_candidate"
        elif 150 <= seqlen <= 400:
            role = "heavy_candidate"
        entities[pid] = dict(desc=desc, seq=seq, seqlen=seqlen, auth_asym=auth_asym, role=role)
    return entities

SOLVENT_IONS = {"HOH", "NA", "CL", "MG", "CA", "ZN", "K", "SO4", "PO4", "GOL", "EDO", "PEG",
                 "CD", "MN", "CU", "NI", "CO", "BR", "I", "ACT", "TRS", "BME", "DMS", "IPA",
                 "MPD", "FMT", "EPE", "IOD", "PG4", "1PE", "MES", "HEPES", "BOG", "P6G", "PGE"}

def parse_atom_site_for_chain(cif_path, chain_id):
    """Parse ATOM records for a given auth_asym_id (label_asym / auth chain).
    Returns dict: residue_key -> set(atom names), ordered list of (resnum, icode), and
    a set of non-solvent HETATM residues seen on this chain (possible modified residues).
    """
    residues = {}  # (resnum, icode) -> set(atom_name)
    order = []
    resname_map = {}
    hetero_nonsolvent = set()
    header = []
    in_loop = False
    col_idx = {}
    try:
        with open(cif_path) as f:
            lines = f.readlines()
    except Exception:
        return None, None

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.strip() == "_atom_site.group_PDB" or (line.startswith("_atom_site.") and not header):
            # start collecting header block
            header = []
            j = i
            while j < n and lines[j].startswith("_atom_site."):
                header.append(lines[j].strip())
                j += 1
            col_idx = {name.split(".", 1)[1]: idx for idx, name in enumerate(header)}
            i = j
            # now read data rows until next loop_/'#'/'_' tag
            while i < n:
                dl = lines[i]
                if dl.startswith("#") or dl.startswith("loop_") or (dl.startswith("_") and not dl.startswith("_atom_site.")):
                    break
                if dl.strip() == "":
                    i += 1
                    continue
                if dl.startswith("_atom_site."):
                    # shouldn't happen after header collected, but guard
                    i += 1
                    continue
                parts = dl.split()
                if len(parts) < len(header):
                    i += 1
                    continue
                try:
                    group = parts[col_idx.get("group_PDB", 0)]
                    auth_asym_id = parts[col_idx["auth_asym_id"]]
                    auth_seq_id = parts[col_idx["auth_seq_id"]]
                    atom_id = parts[col_idx["label_atom_id"]]
                    comp_id = parts[col_idx["label_comp_id"]]
                    icode = parts[col_idx.get("pdbx_PDB_ins_code", col_idx.get("auth_seq_id"))] if "pdbx_PDB_ins_code" in col_idx else "?"
                    altloc = parts[col_idx["label_alt_id"]] if "label_alt_id" in col_idx else "."
                except (KeyError, IndexError):
                    i += 1
                    continue
                if auth_asym_id == chain_id and group == "ATOM":
                    key = (auth_seq_id, icode)
                    if key not in residues:
                        residues[key] = set()
                        order.append(key)
                        resname_map[key] = comp_id
                    residues[key].add(atom_id.strip('"'))
                elif auth_asym_id == chain_id and group == "HETATM" and comp_id not in SOLVENT_IONS:
                    # non-standard/modified residue covalently part of the peptide chain
                    # (e.g. crosslinked or capped residues) - still include in backbone check
                    hetero_nonsolvent.add(comp_id)
                    key = (auth_seq_id, icode)
                    if key not in residues:
                        residues[key] = set()
                        order.append(key)
                        resname_map[key] = comp_id
                    residues[key].add(atom_id.strip('"'))
                i += 1
            continue
        i += 1
    if not order:
        return None, None, hetero_nonsolvent
    return residues, order, hetero_nonsolvent

def analyze_peptide_backbone(cif_path, chain_id, expected_len):
    residues, order, hetero_nonsolvent = parse_atom_site_for_chain(cif_path, chain_id)
    if residues is None:
        return dict(n_residues_observed=0, fully_resolved=False, gap_notes="no atom_site records found for chain")
    n_obs = len(order)
    missing_bb = []
    for key in order:
        atoms = residues[key]
        need = {"N", "CA", "C"}
        if not need.issubset(atoms):
            missing_bb.append(f"{key[0]}{'' if key[1] in ('?','.','') else key[1]}")
    # check numbering gaps (auth_seq_id as int where possible)
    nums = []
    for key in order:
        try:
            nums.append(int(key[0]))
        except ValueError:
            pass
    gap = False
    if nums:
        nums_sorted = sorted(set(nums))
        expected_range = set(range(nums_sorted[0], nums_sorted[-1] + 1))
        if set(nums_sorted) != expected_range:
            gap = True
    fully_resolved = (n_obs == expected_len) and (not missing_bb) and (not gap)
    notes = []
    if n_obs != expected_len:
        notes.append(f"observed {n_obs}/{expected_len} residues")
    if missing_bb:
        notes.append(f"backbone atom(s) missing at residue(s): {','.join(missing_bb)}")
    if gap:
        notes.append("gap in residue numbering")
    if hetero_nonsolvent:
        notes.append(f"non-standard/HETATM residue(s) on peptide chain: {','.join(sorted(hetero_nonsolvent))}")
    return dict(n_residues_observed=n_obs, fully_resolved=fully_resolved, gap_notes="; ".join(notes) if notes else "")

def process_entry(row):
    pdb_id = row["pdb_id"]
    result = dict(row)
    result["download_status"] = None
    result["notes"] = []

    entry_json = api_get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    if entry_json is None:
        result["notes"].append("entry API lookup failed")
        entry_json = {}

    # resolution / release date cross check from API (fallback to csv values)
    info = entry_json.get("rcsb_entry_info", {}) if entry_json else {}
    res_combined = info.get("resolution_combined")
    resolution = res_combined[0] if res_combined else row.get("resolution")
    status = entry_json.get("pdbx_database_status", {}) if entry_json else {}
    release_date = row.get("initial_release_date")

    entities = classify_entities(entry_json, pdb_id) if entry_json else {}

    heavy = [e for e in entities.values() if e["role"] == "heavy"]
    b2m = [e for e in entities.values() if e["role"] == "b2m"]
    pep = [e for e in entities.values() if e["role"] == "peptide"]
    heavy_cand = [e for e in entities.values() if e["role"] == "heavy_candidate"]
    b2m_cand = [e for e in entities.values() if e["role"] == "b2m_candidate"]

    if not heavy and heavy_cand:
        heavy = heavy_cand
        result["notes"].append("heavy chain identified by length only (description lacked expected keywords)")
    if not b2m and b2m_cand:
        b2m = b2m_cand
        result["notes"].append("b2m identified by length only (description lacked expected keywords)")

    n_pep_entities = len(pep)
    if n_pep_entities == 0:
        result["notes"].append("NO PEPTIDE ENTITY IDENTIFIED (8-15aa)")
    elif n_pep_entities > 1:
        result["notes"].append(f"multiple peptide-length entities ({n_pep_entities}); ambiguous")

    if not heavy:
        result["notes"].append("MHC heavy chain not identified")
    elif len(heavy) > 1:
        result["notes"].append(f"multiple candidate heavy-chain entities ({len(heavy)})")
    if not b2m:
        result["notes"].append("beta-2-microglobulin not identified")
    elif len(b2m) > 1:
        result["notes"].append(f"multiple candidate b2m entities ({len(b2m)})")

    heavy_chain_id = heavy[0]["auth_asym"][0] if heavy and heavy[0]["auth_asym"] else ""
    b2m_chain_id = b2m[0]["auth_asym"][0] if b2m and b2m[0]["auth_asym"] else ""
    if heavy and not (260 <= heavy[0]["seqlen"] <= 285):
        result["notes"].append(f"heavy chain entity length atypical ({heavy[0]['seqlen']} aa; expected ~270-280) - possible fusion/tag/truncated construct")
    if b2m and not (95 <= b2m[0]["seqlen"] <= 105):
        result["notes"].append(f"b2m entity length atypical ({b2m[0]['seqlen']} aa; expected ~99) - possible fusion/tag construct")
    peptide_chain_id = pep[0]["auth_asym"][0] if pep and pep[0]["auth_asym"] else ""
    peptide_seq_observed = pep[0]["seq"] if pep else ""

    n_copies_heavy = len(heavy[0]["auth_asym"]) if heavy else 0
    n_copies_pep = len(pep[0]["auth_asym"]) if pep else 0
    n_copies = max(n_copies_heavy, n_copies_pep, 1)
    if n_copies > 1:
        result["notes"].append(f">1 copy of complex in ASU (heavy chains: {heavy[0]['auth_asym'] if heavy else []}, peptide chains: {pep[0]['auth_asym'] if pep else []}); manifest reports FIRST copy only ({peptide_chain_id})")

    # verify csv peptide sequence matches observed peptide entity sequence
    if pep and row.get("peptide") and peptide_seq_observed and row["peptide"] != peptide_seq_observed:
        result["notes"].append(f"csv peptide '{row['peptide']}' != API peptide entity seq '{peptide_seq_observed}'")

    # download cif
    cif_path, dstatus = download_cif(pdb_id)
    result["download_status"] = dstatus
    if cif_path is None:
        result["notes"].append("DOWNLOAD FAILED")

    fully_resolved = False
    bb_notes = ""
    n_obs = ""
    if cif_path and peptide_chain_id:
        try:
            plen = int(row["peptide_len"])
        except Exception:
            plen = len(peptide_seq_observed) if peptide_seq_observed else 0
        analysis = analyze_peptide_backbone(cif_path, peptide_chain_id, plen)
        fully_resolved = analysis["fully_resolved"]
        bb_notes = analysis["gap_notes"]
        n_obs = analysis["n_residues_observed"]
        if bb_notes:
            result["notes"].append(bb_notes)
    elif cif_path and not peptide_chain_id:
        result["notes"].append("cannot assess backbone: no peptide chain id identified")

    result.update(dict(
        resolution=resolution,
        locus=row.get("locus"),
        allele_text=row.get("allele_text"),
        peptide_seq=row.get("peptide"),
        peptide_len=row.get("peptide_len"),
        heavy_chain_id=heavy_chain_id,
        b2m_chain_id=b2m_chain_id,
        peptide_chain_id=peptide_chain_id,
        n_copies_in_asu=n_copies,
        n_residues_observed=n_obs,
        fully_resolved_backbone=fully_resolved,
    ))
    result["notes"] = "; ".join(result["notes"])
    return result

def main():
    sel = load_selection()
    log(f"Selected {len(sel)} entries, {len(set(r['peptide'] for r in sel))} distinct peptides")
    results = []
    with ThreadPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(process_entry, row): row["pdb_id"] for row in sel}
        done_ct = 0
        for fut in as_completed(futs):
            pdb_id = futs[fut]
            try:
                r = fut.result()
            except Exception as e:
                r = dict(pdb_id=pdb_id, notes=f"EXCEPTION: {e}")
            results.append(r)
            done_ct += 1
            if done_ct % 20 == 0:
                log(f"...{done_ct}/{len(sel)} processed")
            time.sleep(0.02)

    # write manifest
    cols = ["pdb_id", "resolution", "release_date", "locus", "allele_text", "peptide_seq",
            "peptide_len", "heavy_chain_id", "b2m_chain_id", "peptide_chain_id",
            "n_copies_in_asu", "notes"]
    results.sort(key=lambda r: r.get("pdb_id", ""))
    with open(MANIFEST, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in results:
            w.writerow([
                r.get("pdb_id", ""),
                r.get("resolution", ""),
                r.get("release_date", r.get("initial_release_date", "")),
                r.get("locus", ""),
                r.get("allele_text", ""),
                r.get("peptide_seq", ""),
                r.get("peptide_len", ""),
                r.get("heavy_chain_id", ""),
                r.get("b2m_chain_id", ""),
                r.get("peptide_chain_id", ""),
                r.get("n_copies_in_asu", ""),
                r.get("notes", ""),
            ])
    # also dump full raw json for downstream stats use
    with open(os.path.join(BASE, "reference_manifest_full.json"), "w") as f:
        json.dump(results, f, indent=1, default=str)
    log(f"Wrote manifest with {len(results)} rows to {MANIFEST}")

if __name__ == "__main__":
    main()
