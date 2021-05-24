# CONVOLVER

A basic gaussian convolver for a 2d spectroscopic or scattering signal, written in julia.

The user supplies (in the script for now) the data, and then the fwhm of the blur (in the same units as the time units), and then an optional pad value, which is defaulted at 4.

The basic algorithm for this relies on a simple approximation to the full convolution - namely, at time t_0 +- 4\sigma, the gaussian is effectively 0.

The user supplies a cutoff value, which is the smallest value that the kernel will contain. For non-periodic functions, this is the edge value.

Therefore, the full gaussian vector which would normally be convolved with the signal can be approximated extremely well as a truncated gaussian, which reduces the number of calculations required significantly.

Also included are Lorentzians, which don't work quite as well, as the function is not as local. You can see this by plotting a gaussian vs a lorentzian, of the same FWHM - the gaussian will be zero well before the lorentzian is, and so this speed up is no longer as effective. A box blur routine is also present.
