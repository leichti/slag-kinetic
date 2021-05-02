from nemlib.chemical import Molecule, Element


if __name__ == "__main__":
    import unittest

    class MoleculeTestCase(unittest.TestCase):

        def test_raise_non_element_in_string(self):
            with self.assertRaises(Exception) as context:
                self.assertTrue(Molecule("Fe3OG").mass)

            self.assertTrue("Symbol G does not represent" in str(context.exception))

        def test_raise_parsed_string_and_string_representation_do_not_match(self):
            with self.assertRaises(Exception) as context:
                self.assertTrue(Molecule("Fe2O3.3").mass)
            self.assertTrue("Only parsed" in str(context.exception))

        def test_raise_concentration_of_non_existing_element(self):
            with self.assertRaises(Exception) as context:
                Molecule("Fe2O3").concentration("Si")

            self.assertTrue("Si not in Fe2O3" in str(context.exception))

        def test_mass_calculation_works(self):
            self.assertAlmostEqual(Molecule("Fe2O3").mass, 159.7, 0)
            self.assertAlmostEqual(Molecule("Fe2OOO").mass, 159.7, 0)
            self.assertAlmostEqual(Molecule("CaF2").mass, 78.1, 0)
            self.assertAlmostEqual(Molecule("Fe").mass, 55.8, 0)

        def test_element_and_molecules_are_the_same(self):
            self.assertEqual(Element("Fe"), Molecule("Fe"))

        def test_molecule_cant_parsed_as_element(self):
            with self.assertRaises(ValueError) as context:
                self.assertTrue(Element("Fe2O3"))
            self.assertTrue("Please use Molecule to parse" in str(context.exception))

        def test_parsing_dict_string(self):
            fe2o3 = Molecule("Fe2O3")
            self.assertEqual(str(fe2o3), "Fe2O3")
            fe2o3_2 = Molecule({"Fe": 2, "O": 3})
            self.assertEqual(str(fe2o3), str(fe2o3_2))

        def test_name_preserved(self):
            fe2o3 = Molecule("FeFeOOO")
            self.assertEqual(fe2o3.name, "FeFeOOO")
            self.assertEqual(str(fe2o3), "Fe2O3")

        def test_string_magic_works(self):
            feo = Molecule("FeO")
            self.assertEqual(str(feo), "FeO")

            zn_ferrite = Molecule("ZnOFe2O3")
            self.assertEqual(str(zn_ferrite), "ZnO4Fe2")

        def test_concentration(self):
            fe2o3 = Molecule("Fe2O3")
            self.assertAlmostEqual(fe2o3.concentration("Fe"), 0.7, 0)

            fe = Element("Fe")
            self.assertAlmostEqual(fe.concentration("Fe2O3"), 0.7, 0)
            p = Element("P")
            self.assertAlmostEqual(p.conversion_factor("P2O5"), 2.29, 0)


    unittest.main()
