import pandas as pd
import numpy as np
from chemical import *


class Analysis(dict):

    def __init__(self):
        """
            Base class of analysis based on python's dict
            inheriting basic functionality like normalize
            to PhaseAnalysis and ElementalAnalysis.
        """
        super().__init__()
        self.name = ""

    def normalize(self, to=1):
        """
            Normalizes the analysis to a specified value

        Parameters
        ----------
        to : float the sum of all values in the analysis after normalization
        """
        total = sum(self.values())

        if total == 0:
            raise ValueError(f"Invalid Analysis {self.name}: {self}")

        for key, old_value in self.items():
            self[key] = old_value/total*to


class ElementalAnalysis(Analysis):

    def __init__(self, data, name=None, as_wt=True, drop=[]):
        """
        ElementalAnalysis represents one chemical analysis in elemental form.
        If data comes as wt.-% it is recalculated to mol.-%.
        After loading, analysis is normalized to 1
        @todo add possibility to transform between mol-% and wt.-%

        Parameters
        ----------
        data : dict analyses of the sample that is represented {"Element" : amount}
        name : str name of the sample. Can be helpful while debugging
        as_wt : boolean defines if parsed data comes as wt.-% or mol.-%
        drop : list elements that will be dropped during initialization. May be used
        for elements that are analyzed inaccurate like  C or O in the EDX
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

    def from_wt(self, data):
        """
        Parses data and assumes it to be given as wt.-%
        Empty data-values (nan values) are set to 0. Otherwise
        functionality like normalizing may be broken

        Parameters
        ----------
        data : dict or dict-like
        """
        for element, wt in data.items():
            wt = 0 if np.isnan(wt) else wt
            self[element] = wt/Element(element).mass

    def from_mol(self, data):
        """
        Parses data and assumes it to be given as mol.-%.
        Empty data-values (nan values) are set to 0. Otherwise
        functionality like normalizing may be broken

        Parameters
        ----------
        data :
        """
        for element, mol in data.items():
            mol = 0 if np.isnan(mol) else mol
            self[element] = mol


class PhaseConstructor:

    def __init__(self, phases, ignore=[]):
        """
        PhaseConstructor mathematically derives chemical compounds from an elemental analysis.
        Building compounds consumes molecules from the element pool (which is an ElementalAnalysis)
        The order in which the phases are given defines which compounds are build first (first->first).

        Parameters
        ----------
        phases : list of compound names that will be builded.
        ignore : list of elements that will not get consumed during phase-building.
                 This is helpful for inaccurate analyzed elements like C and O in the EDX
        """
        self.ignore = ignore
        self.phases = phases
        self.element_pool = ElementalAnalysis

    def parse(self, element_analysis, normalize=1):
        """
        Parses a PhaseAnalysis based on the defined phases in the phases attribute

        Parameters
        ----------
        element_analysis : ElementalAnalysis used as element pool to determine how much of
                           each chemical compound can be build
        normalize : float result is normalized to the given number (default 1)

        Returns
        -------
        PhaseAnalysis instance with the calculated composition
        """
        self.element_pool = element_analysis.copy()
        phase_analysis = PhaseAnalysis()

        for phase in self.phases:
            phase = Molecule(phase)
            self.build(phase, phase_analysis)

        self.check_harvest()
        if normalize:
            phase_analysis.normalize()

        return phase_analysis

    def build(self, phase, phase_analysis):
        """
        Builds the chemical compound  phase and subtracts required elements from the element pool
        (except for elements that are flagged as ignored)

        Parameters
        ----------
        phase : Molecule
            representing one chemical compound (e.g. Al2O3, CaO or SiO2)
        phase_analysis : PhaseAnalysis
            Phase analysis to which the phase is added
        """
        limit = self.build_limit(phase)

        for element, count in phase.items():
            if element in self.ignore:
                continue

            self.element_pool[element] -= limit*count

        phase_analysis[str(phase)] = limit

    def build_limit(self, phase):
        """
        Calculates maximum moles that can be build for the given phase based on the remaining element pool

        Parameters
        ----------
        phase : Molecule
            Phase from which the limit is calculated

        Returns
        -------
        float result of calculation
        """
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
        """
        Short check if the element pool was completely used. If the element pool
        still contains (not ignored) elements, compounds have been left behind.
        """
        ignored_sum = 0
        for element in self.ignore:
            ignored_sum += self.element_pool[element]

        if (total := sum(self.element_pool.values())-ignored_sum) != 0:
            raise ValueError(f"Please mark unused values as ignore. "
                             f"There are still {total} moles in the element pool."
                             f"{self.element_pool}")


class PhaseAnalysis(Analysis):

    def __init__(self):
        """
        Inherits from Analysis base class. Not a lot of use yet.
        Helps determining if an analysis is either elemental or phase-based.
        """
        super().__init__()


class FileParser:

    def __init__(self, path: str):
        """
        Determines the file type by the extension of the given path and
        calls based on that the relevant loading functions implemented in
        child classes of FileParser.

        Parameters
        ----------
        path : str
            relative or absolute path to file that should be loaded
        """
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
        Name_col  |  Element_1  |  ...  |Element_i
        ----------|-------------|-------|---------
        Sample i  |    wt.-%    |  ...  | wt.%

        Parameters
        ----------
        path : str
            path to the file
        """
        super().__init__(path)

        self.elements = only_elements(self.df.columns)
        self.data = self.df[self.elements].to_numpy()

    def load_excel(self):
        """
            Loads data from an excel file. Usually called by the Parent's __init__() method
        """
        self.df = pd.read_excel(f"{self.filename}.{self.ext}", engine="openpyxl", index_col=0)

    def phased(self, phase_constructor):
        """
        Calculates phase analyses for all analyses in the file storing the rows
        in an AnalysesOrganizer object.

        Parameters
        ----------
        phase_constructor : PhaseConstructor

        Returns
        -------
        The phase analyses as AnalysesOrganizer object
        """
        rows = AnalysesOrganizer(phase_constructor.phases)

        for name, row_data in zip(self.df.index, self.data):
            row_data = ElementalAnalysis(dict(zip(self.elements, row_data)), name)
            row_data = phase_constructor.parse(row_data)
            rows.add(row_data, name)

        return rows

    def elemental(self):
        """
        Calculates elemental analyses for all analyses in the file storing the rows
        in an AnalysesOrganizer object.

        Returns
        -------
        The elemental analyses as AnalysesOrganizer object
        """
        rows = AnalysesOrganizer(self.elements)
        for name, row_data in zip(self.df.index, self.data):
            row_data = ElementalAnalysis(row_data, name)
            rows.add(row_data, name)

        return rows


class InfoOrganizer(FileParser):

    def __init__(self, path, lookup_col, read_cols):
        """
        Reads additional information for analyses from a lookup table
        Parameters
        ----------
        path :
        lookup_col :
        read_cols :
        """
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
