import pandas as pd
import numpy as np
from molmass import *

"""
    
    Molmass data taken from the python molmass package
"""
class ElementalAnalysis:

    def __init__(self, data, name, as_wt=True):
        """
        represents one elemental analysis.
        If data comes as wt.-% it is recalculated to mol.-%
        :param data: dict analyses of the sample that is represented {"Element" : amount}
        :param name: str name of the sample that is represented
        :param as_wt: data comes as wt.-%
        """
        if as_wt:
            self.to_mol(data)
        else:
            self.mol = data

    @staticmethod
    def to_mol(input):

        for element, wt in input.items():
            print(molmass.molmass(element))

        return input


class PhaseAnalysis:

    def __init__(self, elemental):
        pass


class SemFileParser():

    def __init__(self, path:str):
        """
        Reads a tabular file, each row representing one elemental analysis.
        A Rows object creates and stores an ElementalAnalysis object for each row.

        :param path: path to file containing sem edx data in the following form
                     Name_col  |  Element_1  |  ...  |Element_i
                     ----------|-------------|-------|---------
                     Sample i  |    wt.-%    |  ...  | wt.%
        """
        *filename, ext = path.split(".")
        if ext in ["xls", "xlsx"]:
            df = pd.read_excel("".join(filename)+f".{ext}", engine="openpyxl", index_col=0)
        else:
            raise ValueError(f"Format {ext} not supported yet")

        columns = [column for column in df.columns if column in ELEMENTS]
        data = df[columns].to_numpy()

        self.rows = Rows(columns)
        for name, row_data in zip(df.index, data):
            self.rows.add(row_data, name)


class Rows(dict):

    def __init__(self, columns):
        self.columns = columns

    def add(self, data, name):
        data = dict(zip(self.columns, data))
        self[name].append(ElementalAnalysis(data, self.columns))

    def __missing__(self, key):
        self[key] = list()
        return self[key]


class RemFile():
    def __init__(self):
        pass

class CompoundAnalysis:
    def init(self):
        pass


if __name__ == "__main__":
    file = SemFileParser("/Users/manuel/OneDrive - Montanuniversität Leoben/Arbeit/Projekte/2021 Green-a-Industry/Versuchsdaten/2020_ZnO-Reduktion/Analysen/allinone.xlsx")
