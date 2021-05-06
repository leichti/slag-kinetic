"""
Only data handling and organization.
No chemical/physical logic.
"""
from nemlib.analysis import PhaseAnalysis


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
        Organizes information regarding the trial in dict form. E.g. date, initial mass, temperature and the like.

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
        pass

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

    def drop(self, sample_id):
        """
        Drops a row. E.g. a sample_idx is contaminated and can't be used for evaluation. Or sampling starts too early.

        Parameters
        ----------
        sample_id : int the index in the list (starting with 0).
        """
        for key in self.keys():
            self[key].pop(sample_id)

    def add_info(self, info_organizer):
        """
        Adds trial information. E.g. the temperature, date or initial mass.

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

    def sample_analysis(self, sample_idx, limit_to=None):
        if limit_to is None:
            limit_to = list(self.phases)
        return PhaseAnalysis(self.sample(sample_idx, limit_to))

    def sample_analyses(self, limit_to=None):
        first_key = next(iter(self.phases))
        for idx, _ in enumerate(self[first_key]):
            yield self.sample_analysis(idx, limit_to)

    def multiply(self, factor):
        for key in self.phases:
            column = self[key]
            self[key] = [v*factor for v in column]

    def multiply_sample(self, factor, sample):
        for key in self.phases:
            self[key][sample] *= factor


class SlagReductionTrial(Trial):
    def __init__(self, data):
        """

        Parameters
        ----------
        data : Selector data that must be parsed with a corresponding parser.

        """
        super().__init__()
        parser = None

        if isinstance(data, Selector):
            parser = SelectorTrialParser(data)

        if parser is None:
            raise("Couldn't find a parser for the given selector {0}".format(type(data)))

        parser.parse(self)
