"""
Which HLA class I supertypes the released alleles cover.

The Editor asked for "a more complete list of HLA subtype". Subtype is answered
by the allele count (37 -> 95 four-digit alleles). Supertype is a different
question that HLA specialists will ask, and the paper had never addressed it:
supertypes group alleles by shared peptide-binding motif, so coverage of them
says whether the resource spans the binding chemistry or just re-samples one
corner of it.

Reference classification: Sidney J, Peters B, Frahm N, Brander C, Sette A.
"HLA class I supertypes: a revised and updated classification." BMC Immunol.
2008;9:1, doi:10.1186/1471-2172-9-1, Additional File 1, checked out into
analysis/reference/sidney2008_hla_supertypes.csv (945 HLA-A and HLA-B alleles).

Two limits of that scheme, both of which must be stated rather than glossed:
  * it covers HLA-A and HLA-B only, so every HLA-C allele here is outside it,
    not "unclassified" by it;
  * 181 of its 945 alleles are themselves "Unclassified", a real category
    meaning no supertype assignment, not missing data.

Usage:
    python3 supertype_coverage.py --out-dir <dir>
Self-test:
    python3 supertype_coverage.py --self-test
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import RELEASE_METADATA  # noqa: E402

REFERENCE = Path(__file__).resolve().parent / "reference" / "sidney2008_hla_supertypes.csv"


def to_sidney(allele: str) -> str:
    """'A*02:01' -> 'A*0201', the four-digit form the reference table uses."""
    a = allele[4:] if allele.startswith("HLA-") else allele
    return a.replace(":", "")


def load_reference(path: Path = REFERENCE) -> pd.DataFrame:
    return pd.read_csv(path, comment="#")


def coverage(alleles, ref: pd.DataFrame) -> dict:
    """Supertype coverage of an iterable of alleles in 'A*02:01' form."""
    m = dict(zip(ref["allele"], ref["supertype"]))
    rows = []
    for a in sorted(set(alleles)):
        locus = a[0]
        st = m.get(to_sidney(a))
        if locus == "C":
            st = "HLA-C (outside the scheme)"
        elif st is None:
            st = "not in the reference table"
        rows.append({"allele": a, "locus": locus, "supertype": st})
    df = pd.DataFrame(rows)
    # A supertype counts as represented only if a real supertype label lands on
    # it; "Unclassified" is an assignment but not a supertype.
    real = df[~df["supertype"].isin(
        ["Unclassified", "HLA-C (outside the scheme)", "not in the reference table"])]
    return {
        "n_alleles": int(len(df)),
        "n_with_supertype": int(len(real)),
        "supertypes_represented": sorted(real["supertype"].unique().tolist()),
        "n_supertypes_represented": int(real["supertype"].nunique()),
        "per_supertype": real["supertype"].value_counts().sort_index().to_dict(),
        "n_unclassified": int((df["supertype"] == "Unclassified").sum()),
        "n_hla_c": int((df["supertype"] == "HLA-C (outside the scheme)").sum()),
        "n_absent_from_reference": int((df["supertype"] == "not in the reference table").sum()),
    }, df


def self_test() -> None:
    ref = pd.DataFrame({
        "allele": ["A*0201", "A*0301", "B*0702", "A*0102"],
        "supertype": ["A02", "A03", "B07", "Unclassified"],
    })
    assert to_sidney("A*02:01") == "A*0201"
    assert to_sidney("HLA-A*02:01") == "A*0201"

    res, df = coverage(["A*02:01", "A*03:01", "B*07:02", "A*01:02", "C*04:01", "B*99:99"], ref)
    assert res["n_alleles"] == 6, res
    assert res["n_supertypes_represented"] == 3, res
    assert res["supertypes_represented"] == ["A02", "A03", "B07"], res
    assert res["n_unclassified"] == 1, res           # A*01:02 is assigned Unclassified
    assert res["n_hla_c"] == 1, res                  # C*04:01 is outside the scheme
    assert res["n_absent_from_reference"] == 1, res  # B*99:99 is in neither
    # An allele the scheme calls Unclassified must not inflate the supertype count.
    assert "Unclassified" not in res["supertypes_represented"]
    print("Self-test PASSED: four-digit conversion is correct, and Unclassified, "
          "HLA-C and absent alleles are each counted separately from a real "
          "supertype assignment.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.out_dir:
        raise SystemExit("--out-dir is required")
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    self_test()
    print()
    ref = load_reference()
    md = pd.read_csv(RELEASE_METADATA, low_memory=False)
    res, df = coverage(md["allele"].unique(), ref)

    print(f"released alleles: {res['n_alleles']}")
    print(f"  with a Sidney 2008 supertype : {res['n_with_supertype']}")
    print(f"  Unclassified by that scheme  : {res['n_unclassified']}")
    print(f"  HLA-C, outside the scheme    : {res['n_hla_c']}")
    print(f"  absent from the table        : {res['n_absent_from_reference']}")
    SINGLE = ["A01", "A02", "A03", "A24", "B07", "B08", "B27", "B44", "B58", "B62"]
    have = [s for s in SINGLE if s in res["per_supertype"]]
    mixed = sorted(s for s in res["per_supertype"] if s not in SINGLE)
    print(f"\nsupertypes represented: {len(have)} of 10 -> {', '.join(have)}")
    if mixed:
        print(f"  plus the mixed categories: {', '.join(mixed)}")
    for st, n in sorted(res["per_supertype"].items()):
        print(f"  {st:<10} {n:>3} alleles")

    df.to_csv(out / "supertype_coverage_per_allele.csv", index=False)
    with open(out / "supertype_coverage.json", "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nwrote supertype_coverage.{{csv,json}} to {out}")


if __name__ == "__main__":
    main()
