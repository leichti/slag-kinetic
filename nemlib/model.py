class SubstanceAmount():
    def __init__(self):
        pass


class MolCourseCalculator:
    def __init__(self):
        pass

    def determine(self):
        raise NotImplementedError("Method determine not implemented")


class TracerCompoundsCalculator(MolCourseCalculator):

    def __init__(self, reference_phases):
        super(TracerCompoundsCalculator, self).__init__()
        self.reference_phases = reference_phases

        # define inert phases where initial mass is known (e.g. CaO, SiO2, Al2O3 added at the beginning)
        # todo we already have moles here.
        # todo This mol sum is assumed to be stable over the experiment.
        # todo In this way we can calc the mol course over our samples. From the mole course we can calculate moles of all oxides (also ZnO and FeO)

        # calculate wt.-% from mol.-% for phases present at initial conditions
        # calculate mass of each compound present at initial conditions
        # calculate mol from mass for compounds present at initial conditions

        # use mol fraction in all analyses from this experiment to calculate amount for all substances

    def determine(self, trial):
        reference_moles = self.reference_moles(trial)
        print(reference_moles)
        for i, sample in enumerate(trial.sample_analyses(self.reference_phases)):
            sample_inert_fraction = sample.sum()
            sample_moles = reference_moles/sample_inert_fraction
            trial.multiply_sample(sample_moles, i)

        return trial

    def reference_moles(self, trial):
        reference_sample = trial.sample_analysis(0, self.reference_phases)

        reference_sample.to_wt()

        reference_sample.multiply(2200)
        reference_sample.to_mol(as_pct=False)

        return reference_sample.sum()


class KineticModel():
    def __init__(self):
        pass


class SingleCompoundModel(KineticModel):
    def __init__(self):
        super().__init__()