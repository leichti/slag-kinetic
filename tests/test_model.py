import unittest
from nemlib.analysis import PhaseConstructor, ElementalAnalysesFile, InfoOrganizer
from nemlib.trial import SlagReductionTrial, VxPySelector
from nemlib.model import TracerCompoundsCalculator, ZnOCCrucibleModel, TrialDimensions
import numpy as np

phase_constructor = PhaseConstructor(phases=["CaF2", "CaO", "SiO2",
                                             "Al2O3", "MgO", "Na2O",
                                             "K2O", "FeO", "ZnO"],
                                     ignore=["O", "C"])
phased_data = ElementalAnalysesFile("test_data/sem.xlsx").phased(phase_constructor)
info = InfoOrganizer("test_data/time_table.xlsx", "sample", ["time"])
phased_data.add_info(info)


class ModelTest(unittest.TestCase):

    def test_mol_course_calculation(self):
        trial = SlagReductionTrial(VxPySelector("V29")(phased_data))
        self.assertAlmostEqual(trial["CaO"][11], 0.315907337)
        sample = trial.sample_analysis(11)
        self.assertAlmostEqual(sample["CaO"], 0.315907337)
        sample.to_wt(as_pct=True)
        self.assertAlmostEqual(sample["CaO"], 0.269576312)

        mol_course = TracerCompoundsCalculator(["CaO", "SiO2",
                                                "Al2O3", "MgO",
                                                "CaF2", "K2O",
                                                "Na2O"], inert_mass=2200).determine(trial)

        # excel says 10.65 because it calculates the reference moles from P12 not P1
        self.assertEqual(mol_course["CaO"][11], 10.566099085297882)

        mol_sum = mol_course.sample_analysis(11).sum()
        self.assertEqual(mol_sum, 33.44683028058214)

    def test_mol_model_fit(self):
        trial = SlagReductionTrial(VxPySelector("V27")(phased_data))
        trial.drop(0)
        trial.drop(0)
        mol_course = TracerCompoundsCalculator(["CaO", "SiO2",
                                                "Al2O3", "MgO",
                                                "CaF2", "K2O",
                                                "Na2O"], inert_mass=2700).determine(trial)
        dimensions = TrialDimensions(270, 2700/(2500)*10)
        v27 = ZnOCCrucibleModel(mol_course, dimensions)
        t = np.arange(0, 300)

        parameters = v27.fit()

        # assertion here to prevent errors from future changes
        v27_zno_model = v27.y(t, *parameters)

        # test that model length is the same as time array length
        self.assertEqual(len(v27_zno_model), len(t))

        # test if model and experimental are similar
        self.assertAlmostEqual(v27_zno_model[0], v27.zno[0], places=0)

        # test that negative values get cut off
        self.assertAlmostEqual(v27_zno_model[-1], 0, places=2)


unittest.main()
