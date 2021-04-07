import re

ELEMENTS = {
    "H": 1.007941,
    "He": 4.002602,
    "Li": 6.94,
    "Be": 9.0121831,
    "B": 10.811,
    "C": 12.01074,
    "N": 14.006703,
    "O": 15.999405,
    "F": 18.998403163,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.3051,
    "Al": 26.9815385,
    "Si": 28.0855,
    "P": 30.973761998,
    "S": 32.0648,
    "Cl": 35.4529,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955908,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938044,
    "Fe": 55.845,
    "Co": 58.933194,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.63,
    "As": 74.921595,
    "Se": 78.971,
    "Br": 79.9035,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90584,
    "Zr": 91.224,
    "Nb": 92.90637,
    "Mo": 95.95,
    "Tc": 97.9072,
    "Ru": 101.07,
    "Rh": 102.9055,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.414,
    "In": 114.818,
    "Sn": 118.71,
    "Sb": 121.76,
    "Te": 127.6,
    "I": 126.90447,
    "Xe": 131.293,
    "Cs": 132.90545196,
    "Ba": 137.327,
    "La": 138.90547,
    "Ce": 140.116,
    "Pr": 140.90766,
    "Nd": 144.242,
    "Pm": 144.9128,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.92535,
    "Dy": 162.5,
    "Ho": 164.93033,
    "Er": 167.259,
    "Tm": 168.93422,
    "Yb": 173.054,
    "Lu": 174.9668,
    "Hf": 178.49,
    "Ta": 180.94788,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.966569,
    "Hg": 200.592,
    "Tl": 204.3834,
    "Pb": 207.2,
    "Bi": 208.9804,
    "Po": 208.9824,
    "At": 209.9871,
    "Rn": 222.0176,
    "Fr": 223.0197,
    "Ra": 226.0254,
    "Ac": 227.0278,
    "Th": 232.0377,
    "Pa": 231.03588,
    "U": 238.02891,
    "Np": 237.0482,
    "Pu": 244.0642,
    "Am": 243.0614,
    "Cm": 247.0704,
    "Bk": 247.0703,
    "Cf": 251.0796,
    "Es": 252.083,
    "Fm": 257.0951,
    "Md": 258.0984,
    "No": 259.101,
    "Lr": 262.1096,
    "Rf": 267.1218,
    "Db": 268.1257,
    "Sg": 271.1339,
    "Bh": 272.1383,
    "Hs": 270.1343,
    "Mt": 276.1516,
}


def only_elements(mixed_list):
    return [column for column in mixed_list if column in ELEMENTS]


class Molecule(dict):

    def __init__(self, input):
        super().__init__()
        self.name = ""

        if (type_check := type(input)) is str:
            self.parse_str(input)
        elif type_check is dict:
            self.parse_dict(input)
        else:
            raise TypeError("Input must be of str or dict type")

    def __missing__(self, key):
        self[key] = 0
        return self [key]

    def parse_str(self, input:str):
        result = re.findall("([A-Z][a-z]*)([0-9]*)", input)

        if not result:
            raise ValueError(f"Invalid input: {input} is neither an element nor a molecule")

        self.parse_tuple(result)

        if self.name != input:
            raise ValueError(f"Only parsed {self.name} from {input}")

        return self

    def parse_dict(self, input):
        self.parse_tuple(input.items())

    def parse_tuple(self, input):
        self.name = ""
        for symbol, count in input:
            self.name += symbol + str(count)
            if symbol not in ELEMENTS:
                raise ValueError(f"Symbol {symbol} does not represent an element")

            if not count:
                count = 1

            self[symbol] += int(count)

    @property
    def mass(self):
        molmass_sum = 0
        for element, count in self.items():
            molmass_sum += ELEMENTS[element] * count

        return molmass_sum

    def __str__(self):
        return self.dict_to_string(self)

    @staticmethod
    def dict_to_string(input_dict):
        return_str = ""
        for element, count in input_dict.items():
            if (count := str(count)) == "1":
                count = ""
            return_str += element+count

        return return_str

    def concentration(self, element):

        # make sure to get element as string
        element = str(element)

        if element not in self.keys():
            raise ValueError(f"{element} not in {self}")

        count = self[element]
        element_mass = Element(element).mass*count
        return element_mass/self.mass


class Element(Molecule):

    def __init__(self, name):
        super().__init__(name)

    def concentration(self, molecule):
        if type(molecule) is str:
            molecule = Molecule(molecule)

        return molecule.concentration(self)

    def conversion_factor(self, molecule):
        if type(molecule) is str:
            molecule = Molecule(molecule)

        return 1/self.concentration(molecule)


if __name__ == "__main__":
    import unittest

    def broken_1():
        Molecule("Fe3OG").mass

    def broken_2():
        Molecule("Fe2O3.3").mass

    def broken_3():
        Molecule("Fe2O3").concentration("Si")

    class MyTestCase(unittest.TestCase):
        def test_fails(self):

            with self.assertRaises(Exception) as context:
                broken_1()

            self.assertTrue("Symbol G does not represent" in str(context.exception))

            with self.assertRaises(Exception) as context:
                broken_2()
            self.assertTrue("Only parsed" in str(context.exception))

            with self.assertRaises(Exception) as context:
                broken_3()

            self.assertTrue("Si not in Fe2O3" in str(context.exception))

        def test_mass_calculation(self):
            self.assertAlmostEqual(Molecule("Fe2O3").mass, 159.7, 0)
            self.assertAlmostEqual(Molecule("Fe2OOO").mass, 159.7, 0)
            self.assertAlmostEqual(Molecule("CaF2").mass, 78.1, 0)
            self.assertAlmostEqual(Molecule("Fe").mass, 55.8, 0)

        def test_element_works(self):
            self.assertEqual(Element("Fe"), Molecule("Fe"))

        def test_parsing_dict_string(self):
            fe2o3 = Molecule("Fe2O3")
            self.assertEqual(str(fe2o3), "Fe2O3")
            fe2o3_2 = Molecule({"Fe":2, "O":3})
            self.assertEqual(str(fe2o3), str(fe2o3_2))

        def test_name_preserved(self):
            fe2o3 = Molecule("FeFeOOO")
            self.assertEqual(fe2o3.name, "FeFeOOO")
            self.assertEqual(str(fe2o3), "Fe2O3")

        def test_string_magic_works(self):
            feo = Molecule("FeO")
            self.assertEqual(str(feo), "FeO")

        def test_concentration(self):
            fe2o3 = Molecule("Fe2O3")
            self.assertAlmostEqual(fe2o3.concentration("Fe"), 0.7, 0)

            fe = Element("Fe")
            self.assertAlmostEqual(fe.concentration("Fe2O3"), 0.7, 0)
            p = Element("P")
            self.assertAlmostEqual(p.conversion_factor("P2O5"), 2.29, 0)

    unittest.main()

