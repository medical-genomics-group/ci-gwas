"""Preparation of phenotype data as input to mps.
"""

from dataclasses import dataclass
import numpy as np
import pandas as pd

@dataclass
class PhenotypesFile:
    filepath: str
    columns: list[str]

    def file_num_cols(self) -> int:
        return len(pd.read_csv(self.filepath, sep=" ", nrows=1, header=None).columns)

    def has_header(self) -> bool:
        first_row = pd.read_csv(self.filepath, sep=" ", nrows=1, header=None)
        return (first_row[0] == "FID").bool() and (first_row[1] == "IID").bool()

    def is_single_phen(self) -> bool:
        return self.file_num_cols() == 3

    def df(self, fam: pd.DataFrame) -> pd.DataFrame:
        assert self.has_header(), "No or incorrect header in pheno file."
        if self.is_single_phen():
            res = pd.read_csv(self.filepath, sep=" ")
            res.columns = ["FID", "IID", *self.columns]
        else:
            res = pd.read_csv(self.filepath, sep=" ", usecols=self.columns)
        res = res.set_index("IID").reindex(index=fam["IID"]).reset_index()
        res.drop(columns=["IID", "FID"], inplace=True)
        return res

@dataclass 
class FamFile:
    filepath: str

    def df(self) -> pd.DataFrame:
        return pd.read_csv(self.filepath, sep=" ", header=None, names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phen'])

def is_standardized(df: pd.DataFrame) -> bool:
    return np.all(np.abs(df.iloc[:, 2:].std(axis=0).values - 1) < 0.1) and np.all(np.abs(df.iloc[:, 2:].mean(axis=0).values) < 0.1)

def merge_phenos(phenos: list[PhenotypesFile], fam: FamFile) -> pd.DataFrame:
    fam_df = fam.df()
    dfs = []

    for p in phenos:
        curr_df = p.df(fam=fam_df)
        assert is_standardized(curr_df), f"data in {p.filepath} seems not precisely standardized"
        dfs.append(curr_df)

    return pd.concat([fam_df[["FID", "IID"]], *dfs], axis=1)

def make_merged_pheno_file(phenos: list[PhenotypesFile], fam: FamFile, outfile: str):
    merge_phenos(phenos, fam).to_csv(outfile, sep="\t", index=False, na_rep="nan")
