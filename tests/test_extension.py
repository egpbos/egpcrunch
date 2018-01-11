import egpcrunch as m
from unittest import TestCase
import numpy as np


class SPHDensityTest(TestCase):

    def test_density_SPH_simple(self):
        dens = m.density_SPH_simple(2, 3, [0.5, 1], [1, 2], [0, 0], 1)
        np.testing.assert_allclose(dens, np.array([[[0.13743393, 0.13743393],
                                                    [0.13743393, 0.13743393]],
                                                   [[0.01133323, 0.01133323],
                                                    [0.01133323, 0.01133323]]]))

    def test_density_SPH(self):
        N = 2
        L = 3.
        dx = L / N
        dens = m.density_SPH(N, N, N, L, L, L, dx, dx, dx, 0, 0, 0, [0.5, 1], [1, 2], [0, 0], [1, 2], True, 1)
        np.testing.assert_allclose(dens, np.array([[[0.14872685, 0.14872685],
                                                    [0.26357494, 0.26357494]],
                                                   [[0.01137354, 0.01137354],
                                                    [0.02262615, 0.02262615]]]))
