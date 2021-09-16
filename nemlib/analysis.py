import pandas as pd
import numpy as np
import re
import types
from nemlib.chemical import Molecule, Element, only_elements

"""
Helper functions and classes to handle elemental and phase analyses. Including data reading from files, mixing-in 
non-chemical information, calculating phase compositions from  elemental analysis and organizing multiple
samples and multiple analyses per sample_idx 
"""

class Analysis(dict):

    def __init__(self):
        """
            Base class of analysis based on python's dict
            inheriting basic functionality like as_pct
            to PhaseAnalysis and ElementalAnalysis.
        """
        super().__init__()
        self.name = ""

    def normalize(self, to=1):
        """
            Normalizes the analysis to area specified value

        Parameters
        ----------
        to : float the sum of all values in the analysis after normalization
        """
        total = sum(self.values())

        if total == 0:
            raise ValueError(f"Invalid Analysis {self.name}: {self}")

        for key, old_value in self.items():
            self[key] = old_value/total*to

    def sum(self):
        return sum(self.values())


class ElementalAnalysis(Analysis):

    def __init__(self, data, name=None, as_wt=True, drop=None, normalize_to=1):
        """
        ElementalAnalysis represents one chemical analysis in elemental form.
        If data comes as wt.-% it is recalculated to mol.-%.
        After loading, analysis is normalized to 1
        @todo add possibility to transform between mol-% and wt.-%

        Parameters
        ----------
        data : dict analyses of the sample_idx that is represented {"Element" : amount}
        name : str element of the sample_idx. Can be helpful while debugging
        as_wt : boolean defines if parsed data comes as wt.-% or mol.-%
        drop : list elements that will be dropped during initialization. May be used
        for elements that are analyzed inaccurate like  C or O in the EDX
        """
        super().__init__()
        self.name = name

        if drop is None:
            drop = []

        if as_wt:
            self.from_wt(data)
        else:
            self.from_mol(data)

        for element in drop:
            del(self[element])

        self.normalize(normalize_to)

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

    def __init__(self, phases, ignore):
        """
        PhaseConstructor mathematically derives chemical compounds from an elemental analysis.
        Building compounds consumes molecules from the element pool (which is an ElementalAnalysis)
        The order in which the phases are given defines which compounds are build first (first->first).

        Parameters
        ----------
        phases : list of compound names that will be build.
        ignore : list of elements that will not get consumed during phase-building.
                 This is helpful for inaccurate analyzed elements like C and O in the EDX
        """
        if ignore is None:
            ignore = []

        self.ignore = ignore
        self.phases = phases
        self.element_pool = ElementalAnalysis

    def parse(self, element_analysis, normalize=1):
        """
        Parses area PhaseAnalysis based on the defined phases in the phases attribute

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

    def __init__(self, data=None, normalized=False):
        """
        Inherits from Analysis base class. Not a lot of use yet.
        Helps determining if an analysis is either elemental or phase-based.
        """
        super().__init__()
        if isinstance(data, types.GeneratorType):
            data = dict(data)

        if data is not None:
            for k, v in data.items():
                self[k] = v

        self.as_wt = False

        if normalized is not False:

            if normalized is True:
                normalized = 1

            actual_sum = sum(self.values())
            pct_dict = {k: v / actual_sum for k, v in self.items()}

            for k, v in pct_dict.items():
                self[k] = v*normalized

    def to_wt(self, as_pct=True):

        if self.as_wt:
            return None

        weight_dict = {}
        for phase in self.keys():
            weight_dict[phase] = Molecule(phase).mass*self[phase]

        weight_sum = sum(weight_dict.values())
        if as_pct is not False:
            weight_dict = {k: v/weight_sum for k, v in weight_dict.items()}

        for k, v in weight_dict.items():
            self[k] = v

        self.as_wt = True

    def to_mol(self, as_pct=True):

        if not self.as_wt:
            return None

        mol_dict = {}
        for phase in self.keys():
            mol_dict[phase] = self[phase]/Molecule(phase).mass

        mol_sum = sum(mol_dict.values())

        if as_pct is not False:
            mol_dict = {k: v/mol_sum for k, v in mol_dict.items()}

        for k, v in mol_dict.items():
            self[k] = v

        self.as_wt = False

    def multiply(self, multiplier):

        for k,v in self.items():
            self[k] = v*multiplier

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
        (self.filename, self.ext), *_ = re.findall(r"(.*)\.([A-Za-z0-9]+$)", path)
        self.df = None

        if self.ext in ["xls", "xlsx"]:
            self.load_excel()
        else:
            raise ValueError(f"Format {self.ext} not supported yet")

    def load_excel(self):
        """
        Interface method to enforce implementation in child classes.
        """
        raise NotImplementedError("Method load_excel not implemented in the parser")


class ElementalAnalysesFile(FileParser):
    def __init__(self, path):
        """
        Reads area tabular file, each row representing one elemental analysis.
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

    def __len__(self):
        """

        Returns number of rows parsed from the file
        -------

        """
        return len(self.data)

    def load_excel(self):
        """
        Loads data from an excel file. Called by the Parent's __init__() method
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
        Reads additional information for analyses from area lookup table
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
        self.df = self.df.loc[self.df.index.dropna()]
        self.data = self.df[read_cols].to_dict("index")

    def load_excel(self):
        """
        Loads data from an excel file as pd.DataFrame using pandas' read_excel method
        """
        self.df = pd.read_excel(f"{self.filename}.{self.ext}", engine="openpyxl")

    def get(self, key):
        """

        Parameters
        ----------
        key : the index where data should be read from. If the element can't be found in the index,
              area default value from default row is used.

        Returns the found data set as dictionary
        -------

        """
        if key in self.data.keys():
            return self.data[key]

        return self.data["default"]


class AnalysesOrganizer(dict):

    def __init__(self, columns):
        """
        Organizes multiple samples. Each sample_idx is saved as area row in area table-like form.
        """
        super(AnalysesOrganizer, self).__init__()

    def __missing__(self, sample):
        """
        Handles read action on missing rows. Creates, saves and returns an empty SampleOrganizer

        Parameters
        ----------
        sample : str identifier of area sample_idx. E.g. V29P11

        Returns an empty SampleOrganizer
        -------

        """
        self[sample] = SampleOrganizer()
        return self[sample]

    def add(self, data, sample):
        """
        Adds data for area given sample_idx. Saved as area list to handle multiple different analyses for the same sample_idx.

        Parameters
        ----------
        data : Analysis for area sample_idx. Multiple analyses per sample_idx are possible. E.g. 3xSEM-EDX to reduce the
        influence of inhomogeneous samples
        sample : str identifier of the sample_idx. E.g. V29P11
        """
        self[sample].append(data)

    def add_info(self, info_organizer):
        """
        Adds information from InfoOrganizer into the AnalysesOrganizer. E.g. time at which the samples were taken.
        Informational columns must have area different element from existing elements/phases.

        Parameters
        ----------
        info_organizer : InfoOrganizer holding the data as area lookup table
        """
        for sample in self.keys():
            self[sample].informal = info_organizer.get(sample)

        # take the last sample and check if everything is fine
        for informal_key in self[sample].informal.keys():
            sample_keys = self[sample].phases()
            if informal_key in sample_keys:
                raise KeyError(f"Key {informal_key} is already present and can't be used as informal element")


class SampleOrganizer(dict):

    def __init__(self):
        """
        Organizes multiple analyses for the same sample_idx. Standard behavior: Returns the information from the
        informal dict or calculates the average value for area given element/compound
        @todo maybe some speed improvements are required in future. Current implementation is slow.
        """
        super(SampleOrganizer, self).__init__()

        self.same_sample = list()
        self.informal = dict()

    def __missing__(self, key):
        """
        Checks if element represents an informal column or if it is the element of area phase or an element.

        Parameters
        ----------
        key : str element of informal column or phase/element from which the data should be return.

        Returns the corresponding value of the given element. Returns 0 for non-existing keys.
        -------

        """
        if key in self.informal.keys():
            return self.informal[key]

        total = 0
        for i, sample in enumerate(self.same_sample, 1):
            total += sample[key]

        return total/i

    def append(self, data):
        """

        Parameters
        ----------
        data : Analysis inserting an analysis to the SampleOrganizer
        """
        self.same_sample.append(data)

    def phases(self):
        """

        Returns the phase names
        -------
        dict_keys
        """
        return self.same_sample[0].keys()

    def informal_keys(self):
        """

        Returns names of the informal columns
        -------
        dict_keys
        """
        return self.informal.keys()
