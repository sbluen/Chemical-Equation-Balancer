import unittest
from subprocess import run


class MyTestCase(unittest.TestCase):
    def test_water(self):
        _INPUT = "H2 + O2 -> H2O\n"
        _OUTPUT = "2H2 + O2 -> 2H2O\n"
        process = run(["python3", "formula_balancer.py"], input=_INPUT, capture_output=True, text=True)
        self.assertEqual(process.stdout, _OUTPUT)

    def test_photosynthesis(self):
        _INPUT = "CO2 + H2O -> C6H12O6 + O2\n"
        _OUTPUT = "6CO2 + 6H2O -> C6H12O6 + 6O2\n"
        process = run(["python3", "formula_balancer.py"], input=_INPUT, capture_output=True, text=True)
        self.assertEqual(process.stdout, _OUTPUT)

    def test_multiletter(self):
        _INPUT = "NaOH + H2CO3 -> Na2CO3 + H2O\n"
        _OUTPUT = "2NaOH + H2CO3 -> Na2CO3 + 2H2O\n"
        process = run(["python3", "formula_balancer.py"], input=_INPUT, capture_output=True, text=True)
        self.assertEqual(process.stdout, _OUTPUT)
        

if __name__ == '__main__':
    unittest.main()
