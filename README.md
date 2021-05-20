# CONVOLVER

A basic gaussian convolver for a 2d spectroscopic or scattering signal, written in julia.

The user supplies (in the script for now) the data, and then the fwhm of the blur (in the same units as the time units), and then an optional pad value, which is defaulted at 4.

The basic algorithm for this relies on a simple approximation to the full convolution - namely, at time t_0 +- 4\sigma, the gaussian is effectively 0.

Therefore, the full gaussian vector which would normally be convolved with the signal can be approximated extremely well as a truncated gaussian, which reduces the number of calculations required significantly.

In the future, a few other useful features will be added (Adding a Lorentzian and Voigt function, adding the ability to use multiple fwhm for the Voigt, a more friendly user experience, etc.)
