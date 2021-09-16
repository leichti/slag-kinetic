import unittest
from nemlib.trial import SlagReductionTrial, VxPySelector
from nemlib.analysis import PhaseConstructor, ElementalAnalysesFile, InfoOrganizer, PhaseAnalysis

phase_constructor = PhaseConstructor(phases=["CaF2", "CaO", "SiO2",
                                             "Al2O3", "MgO", "Na2O",
                                             "K2O", "FeO", "ZnO"],
                                     ignore=["O", "C"])
phased_data = ElementalAnalysesFile("test_data/sem.xlsx").phased(phase_constructor)
info = InfoOrganizer("test_data/info.xlsx", "sample", ["time", "initial mass"])
phased_data.add_info(info)


class SlagTrialTests(unittest.TestCase):

    def setUp(self):
        self.phased_data = phased_data
        self.trial = SlagReductionTrial(data=VxPySelector("V29")(phased_data))

    def test_selector(self):
        returned_data = VxPySelector("V29")(self.phased_data)

        self.assertEqual(returned_data["P1"], self.phased_data["V29P1"])
        self.assertEqual(len(returned_data.keys()), 12)

    def test_trial_sample_count(self):
        self.assertEqual(len(self.trial["CaO"]), 12)

    def test_all_compounds_parsed(self):
        should_contain = ["CaF2", "CaO", "SiO2", "Al2O3",
                          "MgO", "Na2O", "K2O", "FeO", "ZnO", "time"]
        parsed_list = list(self.trial.keys())
        for compound in should_contain:
            self.assertIn(compound, parsed_list)

    def test_trial_samples_load_correct(self):
        cao = self.trial["CaO"]
        self.assertEqual(cao[11], self.phased_data["V29P12"]["CaO"])
        self.assertEqual(cao[0], self.phased_data["V29P1"]["CaO"])

        time = self.trial["time"]
        self.assertEqual(time[11], 75)
        self.assertEqual(time[0], 0)

    def test_samples_multiplication(self):
        self.trial.multiply(100)
        cao = self.trial["CaO"]
        self.assertEqual(cao[11], self.phased_data["V29P12"]["CaO"]*100)
        self.assertEqual(cao[0], self.phased_data["V29P1"]["CaO"]*100)

        time = self.trial["time"]
        self.assertEqual(time[11], 75)
        self.assertEqual(time[0], 0)

    def test_sample_info_inserted(self):
        should_contain = ["time", "initial mass"]
        parsed_list = list(self.trial.keys())
        for compound in should_contain:
            self.assertIn(compound, parsed_list)

    def test_info_failure(self):
        failed_info = InfoOrganizer("test_data/info_failure.xlsx", "sample", ["time", "initial mass", "id"])
        self.phased_data = ElementalAnalysesFile("test_data/sem.xlsx").phased(phase_constructor)
        self.phased_data.add_info(failed_info)

        with self.assertRaises(KeyError) as context:
            self.trial = SlagReductionTrial(data=VxPySelector("V29")(self.phased_data))

        self.assertTrue("Couldn't parse trial." in str(context.exception))

    def test_sample_drop(self):
        trial_id = self.trial["id"]
        self.assertEqual(trial_id[0], "P1")

        self.trial.drop(0)
        trial_id = self.trial["id"]
        self.assertEqual(trial_id[0], "P2")

    def test_sample_drop_multiple(self):
        trial_id = self.trial["id"]
        self.assertEqual(trial_id[0], "P1")

        self.trial.drop_multiple([0,2])
        trial_id = self.trial["id"]
        self.assertEqual(trial_id[0], "P2")
        self.assertEqual(trial_id[1], "P4")

    def test_trial_info_failed(self):
        lookup = InfoOrganizer("test_data/trial_info.xlsx", "trial_id", ["initial mass"])
        with self.assertRaises(KeyError) as e:
            self.trial.add_info(lookup)

        self.assertTrue("Can't add trial info column initial mass that is already present as sample_idx column." in str(e.exception))

    def test_trial_info_added(self):
        lookup = InfoOrganizer("test_data/trial_info.xlsx", "trial_id", ["initial_mass"])
        self.trial.add_info(lookup)
        self.assertEqual(self.trial["initial_mass"], 2200)

    def test_select_single_trial(self):
        sample = PhaseAnalysis(self.trial.sample(11))

        self.assertAlmostEqual(sample["CaO"], 0.315907337)
        self.assertAlmostEqual(sample["SiO2"], 0.392776295)

        sample = PhaseAnalysis(self.trial.sample(11, ["CaO", "SiO2"]))
        self.assertAlmostEqual(sample.sum(), 0.708683632)


    def test_normalization(self):
        not_normalized = self.trial.sample_analysis(1)
        normalized_to_100 = self.trial.sample_analysis(1, normalized=100)

        self.assertEqual(not_normalized["CaO"]*100, normalized_to_100["CaO"])

unittest.main()
