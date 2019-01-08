
import unittest

from snphwe import snphwe

class TestHWE(unittest.TestCase):
    def test_snphwe(self):
        ''' check snphwe gives expected p-values
        '''
        self.assertAlmostEqual(snphwe(500, 10, 5000), 0.65157189991)
        self.assertAlmostEqual(snphwe(1000, 20, 5000), 1.26598491e-5)
    
    def test_snphwe_odd_inputs(self):
        ''' check snphwe with odd inputs
        '''
        # should raise errors with odd inputs
        with self.assertRaises(ValueError):
            snphwe(0, 0, 0)
        with self.assertRaises(ValueError):
            snphwe(-5, 10, 1000)
    
    def test_snphwe_large_input(self):
        ''' check snphwe doesn't give errors with large sample sizes
        '''
        self.assertEqual(snphwe(200000, 200000, 200000), 0.0)
