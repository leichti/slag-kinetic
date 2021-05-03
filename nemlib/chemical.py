import re
"""
Working with elements and phases requires some helper functions and classes to 
get basic chemical data such as mol weight.
"""

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
    """
    Selects only element names from a list and returns them as a list.
    e.g:  ["Cr", "Fe", "time"] returns ["Cr", "Fe"]

    Parameters
    ----------
    mixed_list : list containing a mixture of element names and other values e.g: ["Cr", "Fe", "time"]

    Returns list containing values that refer to element names
    -------

    """
    return [column for column in mixed_list if column in ELEMENTS]


class Molecule(dict):

    def __init__(self, molecule):
        """
        Parses an molecule from str|dict|tuple. The molecule
        is is represented as dictionary. E.g. Fe2O3 becomes {"Fe": 2, "O": 3}
        Allows calculation of the molecule's mol mass.

        Parameters
        ----------
        molecule : str|dict|tuple representing a molecule
        """
        super().__init__()
        self.name = ""

        if (type_check := type(molecule)) is str:
            self.parse_str(molecule)
        elif type_check is dict:
            self.parse_dict(molecule)
        elif type_check is tuple:
            self.parse_tuple(molecule)
        else:
            raise TypeError("Input must be of str or dict type")

    def __missing__(self, element):
        """
        If column the element can't be found in the molecule's dictionary, therefore, 0 is return as occurrence count.

        Parameters
        ----------
        element : str element element

        Returns 0
        -------

        """
        self[element] = 0
        return self[element]

    def parse_str(self, molecule: str):
        """
        Assembles the molecule's dict by parsing a string. Does some validity testing afterwards.

        Parameters
        ----------
        molecule : str representation of the molecule
        -------

        """
        result = re.findall("([A-Z][a-z]*)([0-9]*)", molecule)

        if not result:
            raise ValueError(f"Invalid input: {molecule} is neither an element nor a molecule")

        self.parse_tuple(result)

        if self.name != molecule:
            raise ValueError(f"Only parsed {self.name} from {molecule}")

    def parse_dict(self, molecule):
        """
        Assembles the molecule's dict by parsing a dict. Uses the dictionary items and Molecule parse_tuple method
        Parameters
        ----------
        molecule : dict
        """
        self.parse_tuple(molecule.items())

    def parse_tuple(self, molecule):
        """
        Assembles the molecule's dict by parsing a nested tuple .
        Checks if all given keys are valid elements.

        Parameters
        ----------
        molecule : tuple representing the molecule. e.g. Fe2O3 is (("Fe", 2), ("O", 3))
        """
        self.name = ""
        for symbol, count in molecule:
            self.name += symbol + str(count)
            if symbol not in ELEMENTS:
                raise ValueError(f"Symbol {symbol} does not represent an element")

            if not count:
                count = 1

            self[symbol] += int(count)

    @property
    def mass(self):
        """

        Returns mol mass of the molecule
        -------

        """
        mol_mass_sum = 0
        for element, count in self.items():
            mol_mass_sum += ELEMENTS[element] * count

        return mol_mass_sum

    def __str__(self):
        """

        Returns String representation of the molecule. Caution: ZnOFe2O3 becomes ZnO4Fe3
        -------

        """
        return_str = ""
        for element, count in self.items():
            if (count := str(count)) == "1":
                count = ""
            return_str += element + count

        return return_str

    def concentration(self, element):
        """
        Calculates the relative mol fraction of an element in the molecule.

        Parameters
        ----------
        element : str|Element from which the mol fraction is calculated

        Returns the mol fraction of the element in the molecule
        -------

        """
        # make sure to get element as string
        element = str(element)

        if element not in self.keys():
            raise ValueError(f"{element} not in {self}")

        count = self[element]
        element_mass = Element(element).mass*count
        return element_mass/self.mass


class Element(Molecule):

    def __init__(self, element):
        """
        A simpler representation of Molecule. Also saved as dictionary, but with only 1 column that is the element name

        Parameters
        ----------
        element : str element name
        """
        super().__init__(element)

        if len(self) != 1:
            raise ValueError(f"Please use Molecule to parse {element}")

    def concentration(self, molecule):
        """
        Calculates the concentration of the element in a given molecule. Uses the concentration method from the
        Molecule class.

        Parameters
        ----------
        molecule : str|Molecule

        Returns the concentration of the element in the molecule
        -------

        """
        if type(molecule) is str:
            molecule = Molecule(molecule)

        return molecule.concentration(self)

    def conversion_factor(self, molecule):
        """
        Calculates the conversion factor of an element when recalculating it to a compound. E.g. P to P2O5 has a
        conversion factor of 2.29. Uses the Molecule's concentration method inverse value.

        Parameters
        ----------
        molecule : str|Molecule

        Returns the conversion factor for an element to a compound
        -------

        """
        if type(molecule) is str:
            molecule = Molecule(molecule)

        return 1/self.concentration(molecule)
