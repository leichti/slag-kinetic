import pandas as pd
import numpy as np
from chemical import *


class Analysis(dict):

    def __init__(self):
        super().__init__()
        self.name = ""

    def normalize(self):
        total = sum(self.values())

        if total == 0:
            raise ValueError(f"Invalid Analysis {self.name}: {self}")

        for key, old_value in self.items():
            self[key] = old_value/total


class ElementalAnalysis(Analysis):

    def __init__(self, data, name=None, as_wt=True, drop=[]):
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

        for element in drop:
            del(self[element])
        self.normalize()

    def from_wt(self, input):
        for element, wt in input.items():
            wt = 0 if np.isnan(wt) else wt
            self[element] = wt/Element(element).mass

        return input

    def from_mol(self, input):
        for element, mol in input.items():
            mol = 0 if np.isnan(mol) else mol
            self[element] = mol


class PhaseConstructor():
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

    def __init__(self, phases, ignore=[]):
        self.ignore = ignore
        self.phases = phases
        self.element_pool = ElementalAnalysis

    def parse(self, element_analysis: ElementalAnalysis):
        self.element_pool = element_analysis.copy()
        phase_analysis = PhaseAnalysis()

        for phase in self.phases:
            phase = Molecule(phase)
            self.build(phase, phase_analysis)

        self.check_harvest()
        phase_analysis.normalize()

        return phase_analysis

    def build(self, phase, phase_analysis):
        limit = self.build_limit(phase)

        for element, count in phase.items():
            if element in self.ignore:
                continue

            self.element_pool[element] -= limit*count

        phase_analysis[str(phase)] = limit

    def build_limit(self, phase):
        # check how much moles we can build
        possible_moles = None
        for element, count in phase.items():
            if element in self.ignore:
                continue
            try:
                possible_from_element = self.element_pool[element] / count
            except KeyError:
                possible_from_element = 0
                print(f"Element {element} not found in elemental analysis")

            if possible_moles is None:
                possible_moles = possible_from_element
            if possible_from_element < possible_moles:
                possible_moles = possible_from_element

        return possible_moles

    def check_harvest(self):
        ignored_sum = 0
        for element in self.ignore:
            ignored_sum += self.element_pool[element]

        if (total := sum(self.element_pool.values())-ignored_sum) != 0:
            raise ValueError(f"Please mark unused values as ignore. "
                             f"There are still {total} moles in the element pool."
                             f"{self.element_pool}")


class PhaseAnalysis(Analysis):

    def __init__(self):
        super().__init__()


class FileParser:

    def __init__(self, path: str, **kwargs):

        (self.filename, self.ext), *_ = re.findall("(.*)\.([A-Za-z0-9]+$)", path)
        self.df = None

        if self.ext in ["xls", "xlsx"]:
            self.load_excel()
        else:
            raise ValueError(f"Format {self.ext} not supported yet")


class ElementalAnalysesFile(FileParser):
    def __init__(self, path):
        """
        Reads a tabular file, each row representing one elemental analysis.
        A Rows object creates and stores an ElementalAnalysis object for each row.

        :param path: path to file containing sem edx data in the following form
                     Name_col  |  Element_1  |  ...  |Element_i
                     ----------|-------------|-------|---------
                     Sample i  |    wt.-%    |  ...  | wt.%
        """
        super().__init__(path)

        self.elements = only_elements(self.df.columns)
        self.data = self.df[self.elements].to_numpy()

    def load_excel(self):
        self.df = pd.read_excel(f"{self.filename}.{self.ext}", engine="openpyxl", index_col=0)

    def phased(self, phase_constructor: PhaseConstructor):
        rows = AnalysesOrganizer(phase_constructor.phases)

        for name, row_data in zip(self.df.index, self.data):
            row_data = ElementalAnalysis(dict(zip(self.elements, row_data)), name)
            row_data = phase_constructor.parse(row_data)
            rows.add(row_data, name)

        return rows

    def elemental(self):
        rows = AnalysesOrganizer(self.elements)
        for name, row_data in zip(self.df.index, self.data):
            row_data = ElementalAnalysis(row_data, name)
            rows.add(row_data, name)

        return rows


class InfoOrganizer(FileParser):

    def __init__(self, path, lookup_col, read_cols):
        super().__init__(path)

        if lookup_col not in self.df.columns:
            raise KeyError(f"Given file {path} has no lookup column {lookup_col}")

        for column in read_cols:
            if column not in self.df.columns:
                raise KeyError(f"Given file {path} misses column {column}")

        self.df.index = self.df[lookup_col]
        self.data = self.df[read_cols].to_dict("index")

    def load_excel(self):
        self.df = pd.read_excel(f"{self.filename}.{self.ext}", engine="openpyxl")

    def get(self, key):
        if key in self.data.keys():
            return self.data[key]

        return self.data["default"]


class AnalysesOrganizer(dict):

    def __init__(self, columns):
        self.columns = columns

    def __missing__(self, sample):
        self[sample] = SampleOrganizer()
        return self[sample]

    def add(self, data, name):
        self[name].append(data)

    def informal(self, info_organizer):
        # info organizer has a get method that returns a dict with information
        # for the given index/key
        for sample in self.keys():
            self[sample].informal = info_organizer.get(sample)
        pass


class SampleOrganizer(dict):

    def __init__(self):
        # @todo implement as numpy array! Now, data comes in dict-form..
        # @todo or maybe implement a generator
        self.same_sample = list()
        self.informal = dict()

    def __missing__(self, key):
        if key in self.informal.keys():
            return self.informal[key]

        total = 0
        for i, sample in enumerate(self.same_sample, 1):
            total += sample[key]

        return total/i

    def append(self, data):
        self.same_sample.append(data)


if __name__ == "__main__":
    #file = SemFileParser("../example/data/sem.xlsx")


    import unittest

    class MyTestCase(unittest.TestCase):

        def setUp(self):
            self.data = {
                "C": 16.6, "O": 35.65, "F": 7.14,
                "Na": 0.3, "Mg": 0.3, "Al": 5.2,
                "Si": 11.97, "K": 0.47, "Ca": 21.87,
                "Fe": 0.12, "Zn": 0.39}
            self.analysis = ElementalAnalysis(self.data, drop=["C"])

        def test_elemental_analysis(self):
            analysis = self.analysis
            self.assertAlmostEqual(analysis["O"], 0.58419, 4)
            self.assertAlmostEqual(analysis["Si"], 0.111739334, 4)
            self.assertAlmostEqual(analysis["F"], 0.098531525, 4)


        def test_elemental_analysis_mol(self):
            analysis = ElementalAnalysis(self.data, as_wt=False)
            #self.assertAlmostEqual(analysis["Ca"], 0.625, 2)


        def test_simple_recalculation(self):
            analysis = self.analysis
            phase_analysis = PhaseConstructor(
                                           phases=["CaF2", "CaO",
                                                   "SiO2", "Al2O3",
                                                   "MgO", "Na2O",
                                                   "K2O", "FeO",
                                                   "ZnO"],
                                           ignore=["O"]).parse(analysis)

            self.assertAlmostEqual(phase_analysis["CaO"], 0.324884446, 4)
            self.assertAlmostEqual(phase_analysis["CaF2"], 0.170635578, 4)
            self.assertAlmostEqual(phase_analysis["Al2O3"], 0.087503434, 4)
            self.assertAlmostEqual(phase_analysis["FeO"], 0.001951261, 4)

    class MyFileTestCase(unittest.TestCase):

        def setUp(self):
            phase_constructor = PhaseConstructor(phases=["CaF2", "CaO", "SiO2",
                                                         "Al2O3", "MgO", "Na2O",
                                                         "K2O", "FeO", "ZnO"],
                                                 ignore=["O", "C"])
            self.data = ElementalAnalysesFile("../example/data/sem.xlsx").phased(phase_constructor)

        def test_loading(self):
            self.assertAlmostEqual(self.data["V29P12"]["CaO"], 0.31590732, 5)
            self.assertAlmostEqual(self.data["V29P12"]["CaF2"], 0.170899339, 5)
            self.assertAlmostEqual(self.data["V29P12"]["SiO2"], 0.392776275, 5)

        def test_informal(self):
            info_organizer = InfoOrganizer("../example/data/info.xlsx", "sample", ["time", "initial mass"])
            self.data.informal(info_organizer)
            self.assertAlmostEqual(self.data["V29P12"]["time"], 75, 5)
            self.assertAlmostEqual(self.data["V29P12"]["initial mass"], 2200, 5)
            self.assertAlmostEqual(self.data["V28P0"]["time"], 0, 5)
            self.assertAlmostEqual(self.data["V18P1"]["time"], 0, 5)
            self.assertAlmostEqual(self.data["V18P1"]["initial mass"], 0, 5)


    unittest.main()
