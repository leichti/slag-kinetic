from molmass.elements import ELEMENTS

ELEMENTS
print("ELEMENTS = {")
for element in ELEMENTS:
    print(f"\"{element.symbol}\": {element.mass},")
print("}")