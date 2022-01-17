# CONVOLVER

A basic gaussian convolver for a 2d spectroscopic or scattering signal, written in julia and fortran.

The user supplies (in the script for now) the data, and then the fwhm of the blur (in the same units as the time units), and then an optional pad value, which is defaulted at 4.

The basic algorithm for this relies on a simple approximation to the full convolution - namely, at time t_0 +- 4\sigma, the gaussian is effectively 0.

The user supplies a cutoff value, which is the smallest value that the kernel will contain. For non-periodic functions, this is the edge value.

Therefore, the full gaussian vector which would normally be convolved with the signal can be approximated extremely well as a truncated gaussian, which reduces the number of calculations required significantly.

Also included are Lorentzians, which don't work quite as well, as the function is not as local. You can see this by plotting a gaussian vs a lorentzian, of the same FWHM - the gaussian will be zero well before the lorentzian is, and so this speed up is no longer as effective. A box blur routine is also present.

## Julia

One runs the julia program by running 
```
$ julia CONVOLVER.jl input_file
```

Feel free to modify as you wish

## Fortran

The fortran script must be compiled (e.g. using gnu compilers)
```
gfortran CONVOLVER.f90 -o CONVOLVER.o
```
and then can be run using
```
./CONVOLVER.o
```
The script will take you through the necessary arguments. In the fortran version of the code, the first line of the pdW signal must be the size of the file (number of rows/timesteps then number of columns/frequencies/q points)

# BLAS
 
If one wants maximum performance, one should compile the program using BLAS to utilise the `DDOT` routine.

This version is called `CONVOLVER_BLAS.f90`, and one should refer to the documentation of their BLAS modules to work it out. On my machine (using MKL and intel compilers), the compilation option is
```
ifort CONVOLVER_BLAS.f90  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl  -i8  -I"${MKLROOT}/include" -o CONVOLVER.o
```
after sourcing the compiler variables with 
```
source setvars.sh intel64
```
in the intel OneAPI directory.

You should see a fairly substantial increase in performance (10x increase with GNU compilers).
