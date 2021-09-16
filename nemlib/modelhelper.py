import copy
import math


class TrialDimensions:

    def __init__(self, crucible_inner_diameter, mass, density, option="shell"):
        """
        Represents the experimental dimensions, holds important surface areas
        required for the kinetic modelling.

        Parameters
        ----------
        crucible_inner_diameter : inner diameter in mm.
                                  Used to calculate the bottom area
        mass : slag mass in g. Used to calculated slag height
        density : slag density in g/cm3. Used to calculate slag height
        """
        self.crucible_inner_diameter = crucible_inner_diameter
        self.slag_height = mass/(self.bottom*density/1000)
        self.option = option

    @property
    def bottom(self):
        """

        Returns
        -------
        area : float
            bottom area in mm2 (this is usually the metal-slag surface area)
        """
        return pow(self.crucible_inner_diameter, 2)*math.pi/4

    @property
    def shell(self):
        """

        Returns
        -------
        area : float
            shell area in mm2 (this is usually the crucible-slag surface area)
        """
        return self.crucible_inner_diameter*math.pi*self.slag_height

    @property
    def area(self):
        """
        Returns
        -------
        area : float
            depending on the option setting, the shell, bottom or sum of
            both areas area returned (in mm2)
        """

        if self.option == "shell":
            return self.shell

        if self.option == "bottom":
            return self.bottom

        if self.option == "both":
            return self.shell + self.bottom

        raise ValueError(f"Invalid option value {self.option}")

    @area.setter
    def area(self, value):
        """

        Parameters
        ----------
        value : str
            Option describing the behavior of the area property.
            'bottom' returns the bottom area
            'shell' returns the shell area
            'both' returns the sum of shell and bottom
        """
        self.option = value


class MolCourseCalculator:
    def __init__(self):
        pass

    def determine(self, trial):
        raise NotImplementedError("Method determine not implemented")


class TracerCompoundsCalculator(MolCourseCalculator):

    def __init__(self, inert_phases, inert_mass):
        """
        Uses references phases (inert substances) and the initial mass to
        transform the phase analyses into an mol courses which are required for
        accurate kinetic modelling. TracerCompoundsCalculator has to be setup once,
        then it can just be applied to instances of "Trial" by calling the
        determine method. One problem with this is that this approach is somehow
        overweighting the first sample analysis. However, since most substances have
        similar mol weights, the possible error is small.

        Parameters
        ----------
        inert_phases : Phases that do not react meaning they behave inert during the experiment
        inert_mass : the mass of inert compounds in the first sample.
        """
        super(TracerCompoundsCalculator, self).__init__()
        self.inert_phases = inert_phases
        self.mass = inert_mass

    def determine(self, trial):
        reference_moles = self.reference_moles(trial)

        for i, sample in enumerate(trial.sample_analyses(self.inert_phases)):
            sample_inert_fraction = sample.sum()
            sample_moles = reference_moles/sample_inert_fraction
            trial.multiply_sample(sample_moles, i)

        return trial

    def reference_moles(self, trial):
        reference_sample = trial.sample_analysis(0, self.inert_phases)

        reference_sample.to_wt()

        reference_sample.multiply(self.mass)
        reference_sample.to_mol(as_pct=False)

        return reference_sample.sum()


class ModelCompounds:
    def __init__(self, compound_list: list, trial=None):
        """

        Parameters
        ----------
        compound_list: list
            list of str describing the compounds that are used for kinetic
            modelling in the differential equations
            (in the KineticModel.model() method)
            e.g. ["ZnO", "FeO"]
        trial : SlagReductionTrial
            Experimental representation of a trial building the
            basis for the fitting procedure of KineticModel
        """
        self.compounds = compound_list
        self.trial = trial

    @property
    def name_first(self):
        """

        Returns
        -------
        compound : str
            Name of the first compound in the compound list. E.g. helpful to
            keep the KineticModel.model() method general
        """
        return self.compounds[0]

    @property
    def name_second(self):
        """

        Returns
        -------
        compound : str
            Name of the first compound in the compound list. E.g. helpful to
            keep the KineticModel.model() method general
        """
        return self.compounds[1]

    def name_x(self, x):
        """

        Returns
        -------
        compound : str
            Name of the x's compound in the compound list.
            BE AWARE, x STARTS WITH 0.
            So name_x(1) gives the second compound's name
            E.g. helpful to keep the KineticModel.model() method general
        """
        return self.compounds[x]

    @property
    def first(self):
        """

        Returns
        -------
        values : list
            The experimental values of the first compound from
            the trial-representation class. e.g. SlagReductionTrial
            E.g. helpful to keep the KineticModel.model() method general
        """
        return self.trial[self.name_first]

    @property
    def second(self):
        """

         Returns
         -------
         values : list
             The experimental values of the second compound from
             the trial-representation class. e.g. SlagReductionTrial
             E.g. helpful to keep the KineticModel.model() method general
         """
        return self.trial[self.name_second]

    def x(self, x, pct=False):
        """

         Returns
         -------
         values : list
             The experimental values of the x's compound from
             the trial-representation class. e.g. SlagReductionTrial
             E.g. helpful to keep the KineticModel.model() method general

             BE AWARE, x STARTS WITH 0
         """
        compound = self.name_x(x)
        return self.trial.y(compound, pct=pct)

    @property
    def len(self):
        """

         Returns
         -------
         length : int
             Number of independent compounds used in the kinetic model
             Required to distinguish between y0 and args of model
             E.g. the ZnOFeModel uses ZnO and FeO -> returns 2, so the
             first two values in KineticModel.optimized_parameters property
             are ZnO_0 and FeO_0

         """
        return len(self.compounds)


class TrialAssembler:

    def __init__(self, trial_info=None, inert_compounds=None,
                 mass_col="inert mass", droplets=None,
                 model=None, dimensions=None, slag_density=0,
                 model_compounds=None):

        if inert_compounds is None:
            inert_compounds = []

        if droplets is None:
            droplets = []

        self.trial_info = trial_info
        self.inert_compounds = inert_compounds
        self.mass_col = mass_col
        self.droplets = droplets
        self._model = model
        self.dimensions = dimensions
        self.slag_density = slag_density
        self.model_compounds = model_compounds

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        self._model = model

    def do(self, trial, exclude=None, model=None):
        if exclude is None:
            exclude = self.droplets

        if model is None:
            model = self.model

        trial.model_compounds = copy.copy(self.model_compounds)
        trial.model_compounds.trial = trial

        trial.add_info(self.trial_info)
        trial.drop_multiple(exclude)
        trial.setup_mol_course(self.inert_compounds, self.mass_col)
        trial.setup_model(model, self.dimensions(190, trial[self.mass_col], self.slag_density))

        return trial
