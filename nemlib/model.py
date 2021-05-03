class SubstanceAmount():
    def __init__(self):
        pass


class TracerCompoundsCalculator():

    def __init__(self, trial, inert_phases):

        initial_sample = trial.sample_analysis(0, inert_phases)

        # define inert phases where initial mass is known (e.g. CaO, SiO2, Al2O3 added at the beginning)
        initial_sample.to_wt()
        initial_sample.multiply(2200)
        initial_sample.to_mol(as_pct=False)
        print(initial_sample) # todo we already have moles here. This mol sum is assumed to be stable over the experiment. In this way we can calc the mol course over our samples. From the mole course we can calculate moles of all oxides (also ZnO and FeO)

        # calculate wt.-% from mol.-% for phases present at initial conditions
        # calculate mass of each compound present at initial conditions
        # calculate mol from mass for compounds present at initial conditions

        # use mol fraction in all analyses from this experiment to calculate amount for all substances
        pass


class KineticModel():
    def __init__(self):
        pass


class SingleCompoundModel(KineticModel):
    def __init__(self):
        super().__init__()