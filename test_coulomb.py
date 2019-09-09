"""
Tests for the coulomb module
"""

import unittest
import coulomb as c


class TestCoreFunctions(unittest.TestCase):
    """ Test the core functions """

    def test_ned_to_normal_shear_tangential(self):
        """
        Test rotation to coordinate system of the fault plane
        """
                  #  nn ne nd en ee ed dn de dd
        stresses = [[1, 0, 0, 0, 0, 0, 0, 0, 0],  # 0
                    [1, 0, 0, 0, 0, 0, 0, 0, 0],  # 1
                    [1, 0, 0, 0, 0, 0, 0, 0, 0],  # 2
                    [0, 0, 0, 0, -1, 0, 0, 0, 0], # 3
                    [0, 0, 0, 0, 1, 0, 0, 0, 0],  # 4
                    [1, 0, 0, 0, 0, 0, 0, 0, 0],  # 5
                    [0, 0, 0, 0, -1, 0, 0, 0, 0]] # 6

        str_dip_raks = [[45, 90, 0],   # 0
                        [0,  90, 0],   # 1
                        [90, 45, 0],   # 2
                        [0,  45, 0],   # 3
                        [0,  45, 90],  # 4
                        [90, 45, 90],  # 5
                        [0,  45, -90]] # 6

        nor_shear_tans = [[-0.5, 0.5,  0.0],   # 0
                          [0.0,  1.0,  0.0],   # 1
                          [-0.5, 0.0,  0.5],   # 2
                          [0.5,  0.0, -0.5],   # 3
                          [-0.5, 0.5,  0],     # 4
                          [-0.5, 0.5,  0.0],   # 5 
                          [0.5,  -0.5,  0.0]]   # 6 <- fails

        ii = 0
        for stress, sdr, nst in zip(stresses, str_dip_raks, nor_shear_tans):
            print(ii)
            n, s, t = c.ned_to_normal_shear_tangential(
                stress, sdr[0], sdr[1], sdr[2])
            self.assertAlmostEqual(n, nst[0], delta=1e-15)
            self.assertAlmostEqual(s, nst[1], delta=1e-15)
            self.assertAlmostEqual(t, nst[2], delta=1e-15)
            ii += 1
