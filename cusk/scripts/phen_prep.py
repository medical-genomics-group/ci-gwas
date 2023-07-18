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
        return len(
            pd.read_csv(self.filepath, sep=" ", nrows=1, header=None).columns)

    def canonical_id_colnames(self) -> list[str]:
        first_row = pd.read_csv(self.filepath, sep=" ", nrows=1, header=None)
        field1 = first_row[0].values[0]
        field2 = first_row[1].values[0]
        if self.is_fid_col(field1) and self.is_iid_col(field2):
            return ["FID", "IID"]
        elif self.is_iid_col(field1) and self.is_fid_col(field2):
            return ["IID", "FID"]
        else:
            raise ValueError(f"Header of {self.filepath} is invalid")

    def set_canonical_id_colnames(self, df: pd.DataFrame):
        new_cn = self.canonical_id_colnames()
        old_cn = df.columns
        df.rename(columns={
            old_cn[0]: new_cn[0],
            old_cn[1]: new_cn[1]
        },
                  inplace=True)

    def set_phen_colnames(self, df: pd.DataFrame):
        if self.is_single_phen():
            df.rename(columns={df.columns[2]: self.columns[0]}, inplace=True)

    def phen_colnames(self) -> list[str]:
        if self.is_single_phen():
            return ["PHEN"]
        else:
            return self.columns

    def id_colnames(self) -> list[str]:
        first_row = pd.read_csv(self.filepath, sep=" ", nrows=1, header=None)
        field1 = first_row[0].values[0]
        field2 = first_row[1].values[0]
        return [field1, field2]

    def is_iid_col(self, colname: str) -> bool:
        return colname.upper() == "IID" or colname.upper() == "EID"

    def is_fid_col(self, colname: str) -> bool:
        return colname.upper() == "FID"

    def is_single_phen(self) -> bool:
        return self.file_num_cols() == 3

    def df(self, fam: pd.DataFrame) -> pd.DataFrame:
        res = pd.read_csv(self.filepath, sep=" ")
        self.set_canonical_id_colnames(res)
        self.set_phen_colnames(res)
        res = res[[*self.canonical_id_colnames(), *self.columns]]
        res = res.set_index("IID").reindex(index=fam["IID"]).reset_index()
        res.drop(columns=["IID", "FID"], inplace=True)
        return res


@dataclass
class FamFile:
    filepath: str

    def df(self) -> pd.DataFrame:
        return pd.read_csv(
            self.filepath,
            sep=" ",
            header=None,
            names=["FID", "IID", "Father", "Mother", "Sex", "Phen"],
        )


def is_standardized(df: pd.DataFrame) -> bool:
    return np.all(
        np.abs(df.iloc[:, 2:].std(axis=0).values - 1) < 0.1) and np.all(
            np.abs(df.iloc[:, 2:].mean(axis=0).values) < 0.1)


def merge_phenos(phenos: list[PhenotypesFile], fam: FamFile) -> pd.DataFrame:
    fam_df = fam.df()
    dfs = []

    for p in phenos:
        curr_df = p.df(fam=fam_df)
        assert is_standardized(
            curr_df), f"data in {p.filepath} seems not precisely standardized"
        dfs.append(curr_df)

    return pd.concat([fam_df[["FID", "IID"]], *dfs], axis=1)


def make_merged_pheno_file(phenos: list[PhenotypesFile], fam: FamFile,
                           outfile: str):
    merged_df = merge_phenos(phenos, fam)
    print(f"Writing cols: {merged_df.columns}")
    merged_df.to_csv(outfile, sep="\t", index=False, na_rep="nan")
