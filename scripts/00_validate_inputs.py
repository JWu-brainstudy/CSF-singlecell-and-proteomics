from _bootstrap import pd, setup, write_table
from nmcode.io import sha256


EXPECTED = {
    1: {"CellName", "orig.ident", "Disease group", "Type", "SubType_260501"},
    2: {"cellnames", "cc_aa_alignment", "cc_aa_alignment_size", "chain_pairing"},
    3: {"cellnames", "TCR_dataset", "SubType_260501", "clone_id", "clone_id_size"},
    4: {"ProteinAccessions", "Accession", "PG.Genes", "PG.ProteinDescriptions"},
}


def main():
    _, repo, out = setup("Validate published analysis inputs", "00_validation")
    for number, expected in EXPECTED.items():
        observed = set(repo.source(number, nrows=5).columns)
        missing = expected - observed
        if missing:
            raise ValueError(f"Source Data {number} is missing: {sorted(missing)}")
    rows = [{"file": str(path.relative_to(repo.root)), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in repo.files()]
    write_table(pd.DataFrame(rows), out / "input_manifest.tsv")


if __name__ == "__main__":
    main()
