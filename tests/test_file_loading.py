import unittest
from nemlib.analysis import (ElementalAnalysis,
                             PhaseConstructor,
                             ElementalAnalysesFile,
                             InfoOrganizer)


class AnalysisTestCase(unittest.TestCase):

    def setUp(self):
        self.data = {
            "C": 16.6, "O": 35.65, "F": 7.14,
            "Na": 0.3, "Mg": 0.3, "Al": 5.2,
            "Si": 11.97, "K": 0.47, "Ca": 21.87,
            "Fe": 0.12, "Zn": 0.39}

    def test_elemental_loading_as_wt_drop_element(self):
        analysis = ElementalAnalysis(self.data, drop=["C"])
        self.assertAlmostEqual(analysis["O"], 0.58419, 4)
        self.assertAlmostEqual(analysis["Si"], 0.111739334, 4)
        self.assertAlmostEqual(analysis["F"], 0.098531525, 4)

    def test_normalizing_to_100_when_initialize(self):
        analysis = ElementalAnalysis(self.data, drop=["C"], normalize_to=100)
        self.assertAlmostEqual(analysis["Si"], 11.1739334, 2)

    def test_elemental_loading_as_mol(self):
        analysis = ElementalAnalysis(self.data, as_wt=False)
        self.assertAlmostEqual(analysis["Ca"], 0.218678132, 2)
        self.assertAlmostEqual(analysis["Fe"], 0.00119988, 2)

    def test_exception_for_non_existing_element(self):
        analysis = ElementalAnalysis(self.data)
        with self.assertRaises(KeyError):
            print(analysis["Re"])

    def test_phase_constructor_from_elemental(self):
        analysis = ElementalAnalysis(self.data, drop=["C"])
        phase_analysis = PhaseConstructor(
            phases=["CaF2", "CaO",
                    "SiO2", "Al2O3",
                    "MgO", "Na2O",
                    "K2O", "FeO",
                    "ZnO"],
            ignore=["O"]).parse(analysis)

        self.assertAlmostEqual(phase_analysis["CaO"], 0.324884446, 7)
        self.assertAlmostEqual(phase_analysis["CaF2"], 0.170635578, 7)
        self.assertAlmostEqual(phase_analysis["Al2O3"], 0.087503434, 7)
        self.assertAlmostEqual(phase_analysis["FeO"], 0.001951261, 7)

    def test_phase_analysis_as_wt(self):
        analysis = ElementalAnalysis(self.data, drop=["C"])
        phase_analysis = PhaseConstructor(
            phases=["CaF2", "CaO",
                    "SiO2", "Al2O3",
                    "MgO", "Na2O",
                    "K2O", "FeO",
                    "ZnO"],
            ignore=["O"]).parse(analysis)

        phase_analysis.to_wt()

        self.assertAlmostEqual(phase_analysis["CaO"], 0.277593705, 6)
        self.assertAlmostEqual(phase_analysis["CaF2"], 0.202989364, 6)
        self.assertAlmostEqual(phase_analysis["Al2O3"], 0.135941863, 6)
        self.assertAlmostEqual(phase_analysis["FeO"], 0.002136, 6)

    def test_phase_analysis_double_conversion(self):
        analysis = ElementalAnalysis(self.data, drop=["C"])
        phase_analysis = PhaseConstructor(
            phases=["CaF2", "CaO",
                    "SiO2", "Al2O3",
                    "MgO", "Na2O",
                    "K2O", "FeO",
                    "ZnO"],
            ignore=["O"]).parse(analysis)

        self.assertAlmostEqual(phase_analysis["CaO"], 0.324884446, 6)
        phase_analysis.to_wt()
        self.assertAlmostEqual(phase_analysis["CaO"], 0.277593705, 6)
        phase_analysis.to_mol()
        self.assertAlmostEqual(phase_analysis["CaO"], 0.324884446, 6)

    def test_phase_analysis_as_wt_not_normalized(self):
        analysis = ElementalAnalysis(self.data, drop=["C"])
        phase_analysis = PhaseConstructor(
            phases=["CaF2", "CaO",
                    "SiO2", "Al2O3",
                    "MgO", "Na2O",
                    "K2O", "FeO",
                    "ZnO"],
            ignore=["O"]).parse(analysis)

        phase_analysis.to_wt(as_pct=False)

        self.assertAlmostEqual(phase_analysis["CaF2"], 13.32234017, 6)
        phase_analysis.multiply(10)
        self.assertAlmostEqual(phase_analysis["CaF2"], 133.2234017, 6)


class FileLoading(unittest.TestCase):
    file = ElementalAnalysesFile("test_data/sem.xlsx")

    def setUp(self):
        self.phase_constructor = PhaseConstructor(phases=["CaF2", "CaO", "SiO2",
                                                          "Al2O3", "MgO", "Na2O",
                                                          "K2O", "FeO", "ZnO"],
                                                  ignore=["O", "C"])

    def test_all_rows_columns_of_file_are_loaded(self):
        self.assertEqual(len(self.file), 1060)
        self.assertEqual(len(self.file.elements), 11)

    def test_accessing_and_calculating_mean(self):
        data = self.file.phased(self.phase_constructor)
        self.assertAlmostEqual(data["V29P12"]["CaO"], 0.31590732, 5)
        self.assertAlmostEqual(data["V29P12"]["CaF2"], 0.170899339, 5)
        self.assertAlmostEqual(data["V29P12"]["SiO2"], 0.392776275, 5)

    def test_information_is_added(self):
        data = self.file.phased(self.phase_constructor)

        # add information columns
        info_organizer = InfoOrganizer("test_data/info.xlsx", "sample_idx", ["time", "initial mass"])
        data.informal(info_organizer)

        # specified information
        self.assertAlmostEqual(data["V29P12"]["time"], 75, 5)
        self.assertAlmostEqual(data["V29P12"]["initial mass"], 2200, 5)

        # default information
        self.assertAlmostEqual(data["V28P0"]["time"], 0, 5)
        self.assertAlmostEqual(data["V18P1"]["time"], 0, 5)
        self.assertAlmostEqual(data["V18P1"]["initial mass"], 0, 5)

    def test_CaO_informal_raises(self):
        data = self.file.phased(self.phase_constructor)

        # add information columns
        info_organizer = InfoOrganizer("test_data/info_failure.xlsx", "sample_idx", ["time", "initial mass", "CaO"])

        with self.assertRaises(KeyError) as e:
            data.informal(info_organizer)

        self.assertTrue("Key CaO is already" in str(e.exception))


unittest.main()
