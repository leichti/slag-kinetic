import unittest
from nemlib.analysis import PhaseConstructor, ElementalAnalysesFile, InfoOrganizer
from nemlib.trial import SlagReductionTrial, VxPySelector
from nemlib.model import TracerCompoundsCalculator

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
        self.assertAlmostEqual(trial["CaO"][11], 0.315907337)
        sample = trial.sample_analysis(11)
        self.assertAlmostEqual(sample["CaO"], 0.315907337)
        sample.to_wt(as_pct=True)
        self.assertAlmostEqual(sample["CaO"], 0.269576312)

        mol_course = TracerCompoundsCalculator(["CaO", "SiO2",
                                                "Al2O3", "MgO",
                                                "CaF2", "K2O",
                                                "Na2O"]).determine(trial)

        # excel says 10.65 because it calculates the reference moles from P12 not P1
        self.assertEqual(mol_course["CaO"][11], 10.566099085297882)

        mol_sum = mol_course.sample_analysis(11).sum()
        self.assertEqual(mol_sum, 33.44683028058214)


unittest.main()
