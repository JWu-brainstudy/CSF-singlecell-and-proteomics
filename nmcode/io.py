import hashlib
from pathlib import Path

import pandas as pd


class DataRepository:
    def __init__(self, root):
        self.root = Path(root)
        self.source_root = self.root / "Supplementary Source data"
        self.table_root = self.root / "Supplementary_Tables"

    def source(self, number, **kwargs):
        path = self.source_root / f"Supplementary Source Data {number}.tsv"
        return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)

    def table_path(self, number):
        matches = sorted(self.table_root.glob(f"ST{number:02d}_*.xlsx"))
        if number == 20:
            matches = [path for path in matches if "v0.3" in path.name]
        if len(matches) != 1:
            raise FileNotFoundError(f"Expected one current ST{number:02d} workbook; found {len(matches)}")
        return matches[0]

    def table(self, number, sheet=0):
        return read_table(self.table_path(number), sheet)

    def sheets(self, number):
        path = self.table_path(number)
        return {name: read_table(path, name) for name in pd.ExcelFile(path).sheet_names}

    def files(self):
        sources = sorted(self.source_root.glob("Supplementary Source Data *.tsv"))
        tables = [self.table_path(number) for number in range(1, 40)]
        return sources + tables


def read_table(path, sheet=0):
    raw = pd.read_excel(path, sheet_name=sheet, header=None)
    scan = raw.iloc[: min(12, len(raw))]
    populated = scan.notna().sum(axis=1)
    unique = scan.apply(lambda row: row.dropna().astype(str).nunique(), axis=1)
    header = (populated + unique / 100).idxmax()
    columns = raw.iloc[header].fillna("").astype(str).str.strip()
    table = raw.iloc[header + 1 :].copy()
    table.columns = _unique_columns(columns)
    table = table.dropna(how="all").reset_index(drop=True)
    return table.loc[:, [column for column in table.columns if column]]


def _unique_columns(columns):
    seen = {}
    result = []
    for column in columns:
        if not column:
            result.append("")
            continue
        seen[column] = seen.get(column, 0) + 1
        result.append(column if seen[column] == 1 else f"{column}_{seen[column]}")
    return result


def write_table(table, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
