import pandas as pd
import numpy as np
from chemical import *

class Analysis(dict):

    def __init__(self):
        super().__init__()

    def normalize(self):
        total = sum(self.values())

        for key, old_value in self.items():
            self[key] = old_value/total


class ElementalAnalysis(Analysis):

    def __init__(self, data, name=None, as_wt=True):
        """
        represents one elemental analysis.
        If data comes as wt.-% it is recalculated to mol.-%.
        After loading, analysis is normalized to 1
        :param data: dict analyses of the sample that is represented {"Element" : amount}
        :param name: str name of the sample that is represented
        :param as_wt: data comes as wt.-%
        """
        super().__init__()
        self.name = name

        if as_wt:
            self.from_wt(data)
        else:
            self.from_mol(data)

        self.normalize()

    def from_wt(self, input):
        for element, wt in input.items():
            self[element] = wt/Element(element).mass

        return input

    def from_mol(self, input):
        for element, mol in input.items():
            self[element] = mol


class PhaseAnalysis(Analysis):

    def __init__(self, elemental_analysis, phases, ignore=[]):
        self.element_pool = elemental_analysis
        self.ignore = ignore

        # when recalculating we consider phases in their order
        # check which element of the phase is limiting
        # e.g. CaO requires 1mol Ca and 1mol O. If 2 mol Ca and 3 mol O are available,
        # then Ca is limiting. We can create 2 mol CaO,
        # with 1 mol O remaining in the element pool
        # if one of the elements is in the "ignore" lookup, we don't subtract it
        #
        # after all phases have been produced, we check if all mols from the
        # mol pool have been harvested. If not an exception is thrown
        # Except for elements from the ignore list, they don't need to be harvested
        for phase in phases:
            phase = Molecule(phase)
            self.build(phase)

        self.check_harvest()

        self.normalize()

    def build(self, phase):
        # check how much moles we can build
        possible_moles = None
        limit = self.build_limit(phase)

        for element, count in phase.items():
            if element in self.ignore:
                continue

            self.element_pool[element] -= limit*count

        self[str(phase)] = limit

    def check_harvest(self):
        if (total := sum(self.element_pool.values())) != 0:
            raise ValueError(f"Please mark unused values as ignore. "
                             f"There are still {total} moles in the element pool."
                             f"{self.element_pool}")

    def build_limit(self, phase):
        possible_moles = None
        for element, count in phase.items():
            if element in self.ignore:
                continue
            try:
                possible_from_element = self.element_pool[element]/count
            except KeyError:
                possible_from_element = 0
                print(f"Element {element} not found in elemental analysis")

            if possible_moles is None:
                possible_moles = possible_from_element
            if possible_from_element < possible_moles:
                possible_moles = possible_from_element

        return possible_moles

class SemFileParser:

    def __init__(self, path: str):
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


if __name__ == "__main__":
    #file = SemFileParser("../example/data/sem.xlsx")
    #recalc = CompoundParser(ElementalAnalysis(data), ["Ca2F", "K2O", "CaO", "SiO2", "Al2O3", "ZnO", "FeO"], ["O", "C"])


    import unittest

    class MyTestCase(unittest.TestCase):

        def setUp(self):
            self.data = {"Ca": 10, "Si": 5, "Al": 1}

        def test_elemental_analysis(self):
            analysis = ElementalAnalysis(self.data)
            self.assertAlmostEqual(analysis["Ca"], 0.537, 2)
            self.assertAlmostEqual(analysis["Si"], 0.383, 2)
            self.assertAlmostEqual(analysis["Al"], 0.080, 2)

        def test_elemental_analysis_mol(self):
            analysis = ElementalAnalysis(self.data, as_wt=False)
            self.assertAlmostEqual(analysis["Ca"], 0.625, 2)


        def test_simple_recalculation(self):
            analysis = ElementalAnalysis(self.data)
            phase_analysis = PhaseAnalysis(elemental_analysis=analysis, phases=["CaO", "SiO2", "Al2O3", "SiAl"], ignore=["O"])
            # todo get manually calculated data to check!
            self.assertAlmostEqual(phase_analysis["CaO"], 0, 2)
            self.assertAlmostEqual(phase_analysis["SiO2"], 0, 2)
            self.assertAlmostEqual(phase_analysis["Al2O3"], 0, 2)
            self.assertAlmostEqual(phase_analysis["SiAl"], 0)

    unittest.main()
