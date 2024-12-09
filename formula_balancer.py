"""A script to balance chemical formulae"""

from collections import defaultdict
import fractions
import decimal

import numpy
import scipy


def parse_molecule(formula_string): # noqa: E901
    """
    Returns a dict mapping from element symbol to number of atoms.
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
    """This function formats a coefficient for display. If the coefficient is 1,
    #it returns an empty string ( effectively hiding the 1), otherwise it returns
    #the coefficient as a string."""
    if number == 1:    # pylint: disable=no-else-return
        return ""
    else:
        return str(number)

TOLERANCE = 0.001
def is_almost_whole(number):
    """Returns true if number is or almost is a whole number, with a small delta"""
    number = decimal.Decimal(number)
    return abs(round(number)-number) < TOLERANCE


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

molecule_count = len(left_molecules) + len(right_molecules)  #pylint: disable=invalid-name

element_set = set()
for molecule in left_molecules:
    element_set = element_set.union(molecule.keys())
for molecule in right_molecules:
    element_set = element_set.union(molecule.keys())

element_count = len(element_set)
max_count = max(molecule_count, element_count)

matrix = numpy.zeros(shape=(element_count, molecule_count))

offset = len(left_molecules)     #pylint: disable=invalid-name
for i, element_for_matrix in enumerate(element_set):
    for j, molecule in enumerate(left_molecules):
        matrix[i][j] = molecule[element_for_matrix]
    for j, molecule in enumerate(right_molecules):
        matrix[i][offset + j] = -molecule[element_for_matrix]

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

output = ""    #pylint: disable=invalid-name
for j, string in enumerate(left_molecule_strings):
    output += coeff_format(int(normalized_nullspace[j])) + string + " + "

# -3 to remove the final +
output = output[:-3] + " -> "

offset = len(left_molecules)     #pylint: disable=invalid-name
for j, string in enumerate(right_molecule_strings):
    output += coeff_format(int(normalized_nullspace[offset + j])) + string + " + "

# -3 to remove the final +
output = output[:-3]
print(output)
