"""
Only data handling and organization.
No chemical/physical logic.
"""

from nemlib.analysis import PhaseAnalysis
from nemlib.modelhelper import TracerCompoundsCalculator
import numpy as np


class Selector(dict):
    def __init__(self):
        """
        Parent class for Selectors. Currently a placeholder since we only have the VxPySelector
        """
        super().__init__()


class VxPySelector(Selector):
    def __init__(self, selected_id):
        """
        Selects all samples from trial Vx where the sample_idx identifier column has a VxPy as naming convention.

        Parameters
        ----------
        selected_id : str the trial to select
        """
        super().__init__()
        self.id = selected_id

    def __call__(self, data, *args, **kwargs):
        """
        Extracts the samples of the trial to select.

        Parameters
        ----------
        data : SampleOrganizer from which samples of one trial should be extracted. Keys must have the VxPy naming convention
        args : args
        kwargs : kwargs

        Returns dictionary containing only the samples for the selected trial
        -------

        """
        selected = {k.replace(self.id, ""): v for k, v in data.items() if k.startswith(self.id)}

        for k, v in selected.items():
            self[k] = v

        return self


class SelectorTrialParser:

    def __init__(self, file_data):
        """
        SelectorTrialParser reorganizes a VxPySelector selection for a better handling
        during kinetic modelling. Instead of calling trial[sample_id][compound],
        the reorganization allows to call trial[compound] returning a list with
        all samples. The keyword "id" is preserved.

        Parameters
        ----------
        file_data : selection from an Selector
        """
        self.data = file_data

    def parse(self, trial):
        """
        Inserts reorganized data into the given Trial.

        Parameters
        ----------
        trial : Trial representing one experiment
        """
        trial.id = self.data.id

        for k, analysis in self.data.items():
            self.construct(trial, analysis)

            trial["id"].append(k)

        expected_length = len(analysis.phases()) + len(analysis.informal_keys()) + 1
        real_length = len(trial.keys())

        if expected_length != real_length:
            raise KeyError("Couldn't parse trial. Don't use id as informal keys."
                           "Check if informal keys have the same"
                           "value as phases from the analysis file.")

    @staticmethod
    def construct(trial, analysis):
        """
        Helper method to construct the trial from a phase|element analysis. Takes care of the informal columns.

        Parameters
        ----------
        trial : Trial to assemble.
        analysis : SampleOrganizer from which we extract data.
        """
        trial.phases = analysis.phases()
        trial.sample_info = analysis.informal_keys()
        for phase in analysis.phases():
            trial[phase].append(analysis[phase])

        for info_key in analysis.informal_keys():
            trial[info_key].append(analysis[info_key])


class TrialInfo(dict):
    def __init__(self, data=None):
        """
        Organizes information regarding the trial in dict form. E.g. date, initial initial_mass, temperature and the like.

        Parameters
        ----------
        data : dict
        """
        if data is None:
            data = {}

        super(TrialInfo, self).__init__()
        for key in data.keys():
            self[key] = data[key]


class Trial(dict):
    def __init__(self):
        """
        Base class of an experiment. Trial organizes data in table-like form. The columns are phase or element names,
        the rows are the corresponding values saved in a list.
        """
        super().__init__()
        self.id = None
        self.info = TrialInfo()
        self.sample_info = []
        self.phases = []
        self.model_compounds = None

    def __missing__(self, column):
        """
        If a column doesn't exist, an empty list is created and return.

        Parameters
        ----------
        column :

        Returns
        -------

        """
        if column in self.info.keys():
            return self.info[column]

        self[column] = list()
        return self[column]

    def drop_multiple(self, sample_ids):
        if type(sample_ids) is not list:
            raise TypeError("sample_ids must be of type list")

        sample_ids.sort()
        subtract = 0

        for id in sample_ids:
            self.drop(id-subtract)
            subtract += 1

    def drop(self, sample_id):
        """
        Drops a row. E.g. a sample_idx is contaminated and can't be used for evaluation. Or sampling starts too early.

        Parameters
        ----------
        sample_id : int or list the index in the list (starting with 0).
        """
        for key in self.keys():
            self[key].pop(sample_id)

    def add_info(self, info_organizer):
        """
        Adds trial information. E.g. the temperature, date or initial initial_mass.

        Parameters
        ----------
        info_organizer : InfoOrganizer lookup table from with information is added
        """
        info = TrialInfo(info_organizer.get(self.id))

        for key in info.keys():
            if key in self.keys():
                raise KeyError(f"Can't add trial info column {key} that is already present as sample_idx column.")

        self.info = info

    def sample(self, sample_idx, limit_to=None):

        if limit_to is None:
            for k in self.keys():
                yield k, self[k][sample_idx]
            return

        for k in limit_to:
            yield k, self[k][sample_idx]

    def sample_analysis(self, sample_idx, limit_to=None, normalized=False):
        if limit_to is None:
            limit_to = list(self.phases)
        return PhaseAnalysis(self.sample(sample_idx, limit_to), normalized=normalized)

    def sample_analyses(self, limit_to=None, normalized=False):
        first_key = next(iter(self.phases))
        for idx, _ in enumerate(self[first_key]):
            yield self.sample_analysis(idx, limit_to=limit_to, normalized=normalized)

    def concentration_course(self, compound):
        trial_analyses = self.sample_analyses(limit_to=None, normalized=100)

        for analysis in trial_analyses:
            yield analysis[compound]

    def multiply(self, factor):
        """

        Parameters
        ----------
        factor : the factor every analysis of all samples are multiplied with

        """
        for key in self.phases:
            column = self[key]
            self[key] = [v*factor for v in column]

    def multiply_sample(self, factor, sample):
        """

        Parameters
        ----------
        factor : the factor the sample is multiplied with
        sample : idx of the sample that will be multiplied with the factor
        """
        for key in self.phases:
            self[key][sample] *= factor

    def sw_sum(self, limit_to=[], exclude=[]):
        """
        Calculates the row-wise (sample-wise) sum of the phase analyses.

        Returns list containig sample-wise sum of the analyses (note that "analyses" can also
        mean absolute mole numbers)
        -------

        Parameters
        ----------
        exclude : substances to be excluded from summarization
        limit_to : substances to be summed up

        """
        first_key = next(iter(self.keys()))
        sample_sum = [0]*len(self[first_key])

        for key in self.phases:
            if key not in limit_to and limit_to:
                continue

            if key in exclude:
                continue

            for idx, value in enumerate(self[key]):
                v = self[key][idx]
                v2 = sample_sum[idx]
                sample_sum[idx] += self[key][idx]

        return sample_sum

    def t_shift(self, to_shift):
        for i, v in enumerate(self[self.time_col]):
            self[self.time_col][i] = v + to_shift

    @property
    def t_min(self):
        return min(self[self.time_col])

    @property
    def t_max(self):
        return max(self[self.time_col])

    @property
    def t_min_max(self):
        return self.t_min, self.t_max

    @property
    def t(self):
        return self[self.time_col]

class SlagReductionTrial(Trial):
    def __init__(self, data):
        """

        Parameters
        ----------
        data : Selector data that must be parsed with a corresponding parser.

        """
        super().__init__()
        parser = None
        self.model = None
        self.mol_course = None
        self.fitted = False
        self.parameters = []
        self.time_col = None
        self.t_eval = None
        self._inert_moles = None

        if isinstance(data, Selector):
            parser = SelectorTrialParser(data)

        if parser is None:
            raise("Couldn't find a parser for the given selector {0}".format(type(data)))

        parser.parse(self)

    def setup_model(self, model, dimensions, time_col="time"):
        """

        Parameters
        ----------
        time_col : the column name for stored time data
        model : KineticModel a callable kinetic model class.
        dimensions : TrialDimensions mapping the dimensions of the trial
        """
        self.time_col = time_col
        self.t_eval = np.arange(min(self[self.time_col]), max(self[self.time_col]))
        self.model = model(self, dimensions)

    def setup_mol_course(self, inert_compounds, mass_column="inert mass"):
        """
        Multiplies all mol analyses based on the reference mass and the
        inert substances to generate a mol course instead of mol analyses.
        Since the Trial object is not copied, all changes are applied to self

        Parameters
        ----------
        inert_compounds : inert compounds in the slag
        mass_column : column holding the information of the experimental's weighted slag mass
        """
        TracerCompoundsCalculator(inert_compounds, inert_mass=self[mass_column]).determine(self)

    def fit(self):
        self.model.optimized_parameters = self.model.fit()
        self.fitted = True

    def model_y(self, t=False, pct=False, which=0):
        """
        Function to get values vs time from the fitted model. Useful for drawing

        Parameters
        ----------
        t :

        Returns
        -------

        """
        if self.fitted is False:
            self.fit()

        if t is False:
            t = self.t_eval

        if pct is True:
            y = self.model.y_pct(t, *self.model.optimized_parameters)
        else:
            y = self.model.y(t, *self.model.optimized_parameters)

        if type(which) is int:
            return t, y[which]

        return t, y

    def y(self, compound, pct=False):
        result = self[compound]

        if pct is True:
            sw_sum = self.sw_sum()
            result = np.divide(self[compound], sw_sum) * 100

        return result

    def xy(self, compound, pct=False):
        return self[self.time_col], self.y(compound, pct)

    @property
    def inert_moles(self):
        """
        Every sample should have the same concentration of inert moles.
        However, random deviations of the analyses generate some noise.
        Therefore, the average amount of inert substances of all samples is taken.

        Returns mol amount of inert substances
        -------

        """
        if self._inert_moles is None:
            v = self.sw_sum(exclude=["ZnO", "FeO"])
            self._inert_moles = np.average(v)

        return self._inert_moles


