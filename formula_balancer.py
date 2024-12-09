from collections import defaultdict
import numpy
import scipy
import fractions
import decimal


def parse_molecule(formula_string):
    """
    Returns the formula string split into elements and quantites.
    formula_string is the text representing the molecule, such as C6H12O6
    """
    elements = defaultdict(int)
    element = ""
    subscript = ""
    state = "expecting_upper"
    for char in formula_string:
        if state == "expecting_upper":
            assert char.isupper()
            element += char
            state = "expecting_letter_or_digit"
        elif state == "expecting_letter_or_digit":
            assert char.isalpha() or char.isdigit
            # We've just read the capitial letter of the element symbol
            if char.isupper():
                # No subscript means that the subscript defaults to 1 atom
                elements[element] += 1
                element = char
                subscript = ""
            elif char.islower():
                element += char
            elif char.isdigit():
                subscript += char
                state = "expecting_upper_or_digit"
        elif state == "expecting_upper_or_digit":
            # We've just read a digit of the molecular formula subscript
            assert char.isupper() or char.isdigit()
            if char.isupper():
                elements[element] += int(subscript)
                element = char
                subscript = ""
                state = "expecting_letter_or_digit"
            elif char.isdigit():
                subscript += char
    if element:
        # final element
        if not subscript:
            subscript = 1
        elements[element] += int(subscript)
        # cleanup not necessary here because we're returning
    return elements


def coeff_format(number):
    """Returns the input in coefficient format such that 1 becomes the empty string
    and anything else becomes the str of itself"""
    if number == 1:
        return ""
    else:
        return str(number)


def is_almost_whole(number):
    """Returns true if number is or almost is a whole number, with a small delta"""
    DELTA = 0.001
    number = decimal.Decimal(number)
    return abs(round(number)-number) < DELTA


row = input()
left_raw, right_raw = row.split(" -> ")

left_molecules = []
left_molecule_strings = []
for symbol in left_raw.split():
    if symbol == "+":
        continue
    left_molecules.append(parse_molecule(symbol))
    left_molecule_strings.append(symbol)

right_molecules = []
right_molecule_strings = []
for symbol in right_raw.split():
    if symbol == "+":
        continue
    right_molecules.append(parse_molecule(symbol))
    right_molecule_strings.append(symbol)

molecule_count = len(left_molecules) + len(right_molecules)

element_set = set()
for molecule in left_molecules:
    element_set = element_set.union(molecule.keys())
for molecule in right_molecules:
    element_set = element_set.union(molecule.keys())

element_count = len(element_set)
max_count = max(molecule_count, element_count)

matrix = numpy.zeros(shape=(element_count, molecule_count))

offset = len(left_molecules)
for i, element in enumerate(element_set):
    for j, molecule in enumerate(left_molecules):
        matrix[i][j] = molecule[element]
    for j, molecule in enumerate(right_molecules):
        matrix[i][offset+j] = -molecule[element]
        
# matrix is now a linear algebra matrix and we need to solve for the case
# when matrix equals the 0 vector
nullspace = scipy.linalg.null_space(matrix)
normalized_nullspace = nullspace / min((abs(i[0]) for i in nullspace))

# Now we need the answers to be in whole numbers
while not all((is_almost_whole(x[0]) for x in normalized_nullspace)):
    for i in normalized_nullspace:
        if not is_almost_whole(i[0]):
            normalized_nullspace *= fractions.Fraction.from_float(round(i[0], 5)).denominator
            break

normalized_nullspace = numpy.round(normalized_nullspace, decimals=0)

# Doing this to make sure the results are not negative
sign = int(normalized_nullspace[0] / abs(normalized_nullspace[0]))
normalized_nullspace *= sign

output = ""
for j, string in enumerate(left_molecule_strings):
    output += coeff_format(int(normalized_nullspace[j])) + string + " + "

# -3 to remove the final +
output = output[:-3] + " -> "

offset = len(left_molecules)
for j, string in enumerate(right_molecule_strings):
    output += coeff_format(int(normalized_nullspace[offset+j])) + string  + " + "

# -3 to remove the final +
output = output[:-3]
print(output)

