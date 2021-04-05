import pandas as pd
import numpy as np
from molmass import *

"""
    
    Molmass data taken from the python molmass package
"""
class ElementalAnalysis:

    def __init__(self, data, name=None, as_wt=True):
        """
        represents one elemental analysis.
        If data comes as wt.-% it is recalculated to mol.-%
        :param data: dict analyses of the sample that is represented {"Element" : amount}
        :param name: str name of the sample that is represented
        :param as_wt: data comes as wt.-%
        """
        self.name = name

        if as_wt:
            self.to_mol(data)
        else:
            self.mol = data

    @staticmethod
    def to_mol(input):
        # @todo implement
        for element, wt in input.items():
            print(molmass(element))

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

        (filename, ext), *_ = re.findall("(.*)\.([A-Za-z0-9]+$)", path)

        if ext in ["xls", "xlsx"]:
            df = pd.read_excel(f"{filename}.{ext}", engine="openpyxl", index_col=0)
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


class CompoundParser(dict):

    def __init__(self, elemental_analysis, compounds, ignore):
        """
        Check if ignore elements are real elements!
        Check if compounds can be parsed to elements

        Parameters
        ----------
        elemental_analysis basis from which a compound analysis is calculated
        compounds
        ignore elements that are not limited by the pool. e.g. C and O analyses
               are not accurate which is why they are not taken into considerations.
        """
        if type(elemental_analysis) is not ElementalAnalysis:
            raise ValueError(f"Elemental analysis must be of type ElementalAnalysis. "
                             f"Got {type(elemental_analysis)}")
        self.ignore = ignore

        self.mol_pool = elemental_analysis

        for compound in compounds:
            compound = parse_compound(compound)
            # @todo implement, refactor. All the simple stuff regarding splitting of strings into compounds has to be
            # @todo outsourced from analysis (because we can reuse this very often)
            pass




if __name__ == "__main__":
    file = SemFileParser("../example/data/sem.xlsx")

    data = {"Ca" : 10, "Si": 5}
    recalc = CompoundParser(ElementalAnalysis(data), ["Ca2F", "K2O", "CaO", "SiO2", "Al2O3", "ZnO", "FeO"], ["O", "C"])
