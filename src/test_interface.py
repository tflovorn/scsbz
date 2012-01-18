import unittest
import interface

class ValidateParamsTest(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(interface.validateParams([{}]), False)

    def test_full(self):
        fullParams = [{"zoneLength": 64, "t": 1.0, "tc": 0.1, "beta": 1.0, "x": 0.1, "J": 0.25}]
        self.assertEqual(interface.validateParams(fullParams), True)
        fullParams.append(fullParams[0])
        self.assertEqual(interface.validateParams(fullParams), True)

    def test_almost_full(self):
        almostFullParams = [{"t": 1.0, "tc": 0.1, "beta": 1.0, "x": 0.1, "J": 0.25}]
        self.assertEqual(interface.validateParams(almostFullParams), False)

if __name__ == "__main__":
    unittest.main()
