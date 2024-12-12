
import numpy as np

class Cosmology:
    def __init__(self,para):
        self._z = para['z']
        self._omega_m = para['omega_m']
        self._omega_l = para['omega_l']
        self._omega_k = para['omega_k']
        self._h = para['h']

        self._H0 = para['h0']
        self._w0 = para['w0']
        self._wa = para['wa']


    def H(self, a):
            return self._h*self._H0*np.sqrt(self.Esqr(a))

    def Esqr(self, a):
            return self._omega_m*pow(a, -3) + self._omega_k*pow(a, -2) + self._omega_l*pow(a, self.f_de(a))

    def f_de(self, a):
                  r"""Evolution parameter for the Dark Energy density.
                  Parameters
                  ----------
                  a : array_like
                  Scale factor
                  Returns
                  -------
                  f : ndarray, or float if input scalar
                  The evolution parameter of the Dark Energy density as a function
                  of scale factor
                  Notes
                  -----
                  For a given parametrisation of the Dark Energy equation of state,
                  the scaling of the Dark Energy density with time can be written as:
                  .. math::
                  \rho_{de}(a) \propto a^{f(a)}
                  (see :cite:`2005:Percival`) where :math:`f(a)` is computed as
                  :math:`f(a) = \frac{-3}{\ln(a)} \int_0^{\ln(a)} [1 + w(a^\prime)]
                  d \ln(a^\prime)`. In the case of Linder's parametrisation for the
                  dark energy in Eq. :eq:`linderParam` :math:`f(a)` becomes:
                  .. math::
                  f(a) = -3(1 + w_0) + 3 w \left[ \frac{a - 1}{ \ln(a) } - 1 \right]
                  """
                  # Just to make sure we are not diving by 0
                  epsilon = 0.000000001
                  return -3.0*(1.0+self._w0) + 3.0*self._wa*((a-1.0)/np.log(a-epsilon) - 1.0)
