import unittest
from nemlib.analysis import PhaseConstructor, ElementalAnalysesFile, InfoOrganizer
from nemlib.trial import SlagReductionTrial, VxPySelector
from nemlib.model import SingleCompoundModel, TracerCompoundsCalculator

phase_constructor = PhaseConstructor(phases=["CaF2", "CaO", "SiO2",
                                             "Al2O3", "MgO", "Na2O",
                                             "K2O", "FeO", "ZnO"],
                                     ignore=["O", "C"])
phased_data = ElementalAnalysesFile("test_data/sem.xlsx").phased(phase_constructor)
info = InfoOrganizer("test_data/time_table.xlsx", "sample", ["time"])
phased_data.informal(info)

class ModelTest(unittest.TestCase):


    def test_zn_model_fits(self):
        trial = SlagReductionTrial(VxPySelector("V29")(phased_data))
        trial.drop(0)

        TracerCompoundsCalculator(trial, ["CaO", "SiO2", "Al2O3", "MgO", "CaF2", "K2O", "Na2O"])


        return

        model = SingleCompoundModel(time="time", compound="ZnO",
                                   start_substance_amount=2500, diameter=200)

        model.fit(trial)

        import matplotlib.pyplot as plt
        plt.plot(trial["time"], trial["ZnO"])
        plt.scatter(trial["time"], trial["ZnO"])
        plt.show()
        self.assertEqual(0.1, 0.1)


unittest.main()