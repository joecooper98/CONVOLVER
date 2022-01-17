! J. C. Cooper 2022
!
! Convolves a pump-probe style signal
! i.e. convolves in time, leaves q/nu/nm/whatever alone
!
! Works with truncated odd-sized kernels to make easier
! Provide a cutoff value, which is the smallest allowed values in the kernel
!
! Will pad the borders (t<0, t>final) with either zeros or the values at t=0, t=final
!
! Example - the kernel only reaches to the cutoff value to limit calculations
!
! 0 0 1 2 3 4 5 ...
! 1 4 6 4 1
!
! 0 0 1 2 3 4 5 ...
!   1 4 6 4 1
!
! 0 0 1 2 3 4 5 ...
!     1 4 6 4 1
!

PROGRAM CONVOLVER
   implicit none

   !initialise variables

   double precision :: fwhm, cutoff, dt, temp, kernsum
   double precision, dimension(:), allocatable :: kernel
   double precision, dimension(:, :), allocatable :: signal
   logical :: extendforward = .false., extendbackward = .false., lorz = .false., step = .false.
   integer :: nt, nq, i, j, length
   character(len=50) :: filename
   character :: inp
   double precision :: start, finish, tstart, tfinish

   WRITE (*, *) 'FWHM /fs ='
   READ (*, *) fwhm
   WRITE (*, *) 'Cutoff value ='
   READ (*, *) cutoff
   WRITE (*, *) 'dt / fs ='
   READ (*, *) dt
   WRITE (*, *) "Use Lorentzian (l) or Box (b)? default Gaussian"
   READ (*, *) inp
   WRITE (*, *) "File name?"
   READ (*, *) filename

   call cpu_time(tstart)


   if (inp .EQ. 'l') then
      lorz = .true.
   else if (inp .EQ. 'b') then
      step = .true.
   end if

   call cpu_time(start)
   write (*, *) ''
   write (*, *) ''
   write (*, *) "GENERAL INFORMATION"
   write (*, *) ''
   print '(" FWHM / fs       = ", ES10.3)', fwhm
   print '(" Cutoff          = ", ES10.3)', cutoff
   print '(" dt / fs         = ", ES10.3)', dt
   print '(" Data file       = ",a10)', trim(filename)

   !read in data. from file pdW with time as rows and q (or nm/nu) as columns
   OPEN (1, file=trim(filename))
   READ (1, *) nt, nq
   allocate (signal(nt, nq))
   read (1, *) ((signal(i, j), j=1, nq), i=1, nt)

   !read(1,*) signal
   !signal = transpose(signal)

   print '(" Points in time  = "I10)', nt
   print '(" Points in q     = "I10)', nq
   print '(" Points in total = "I10)', nt*nq

   call cpu_time(finish)

   print '(" Signal loading took ", f9.6, " s")', finish - start
   call cpu_time(start)
   write (*, *) ''
   write (*, *) ''
   write (*, *) "KERNEL INFORMATION"
   write (*, *) ''

   ! Takes kernel, and allocates it. Can be whatever you want here.
   ! Important to note that kernel must be odd-sized - massively improves numerical procedure.
   if (lorz) then
      CALL LORENTZIAN(kernel, fwhm, cutoff, dt)
      write (*, *) "Kernel type     = Lorentzian"
   else if (step) then
      CALL box(kernel, fwhm, dt)
      write (*, *) "Kernel type     =        Box"
   else
      CALL GAUSSIAN(kernel, fwhm, cutoff, dt)
      write (*, *) "Kernel type     =   Gaussian"
   end if

   kernsum=sum(kernel)

   print '(" Kernel size     = ", I10)', size(kernel)
   print '(" Kernel cutoff   = ", ES10.3)', cutoff
   print '(" Kernel sum      = ", ES10.3)', kernsum
   print '(" Kernel error    = ", ES10.3)', 1-kernsum
   !If the cutoff is too large you won't approximate the kernel well

   if (kernsum .LT. 0.95) then
      WRITE (*, *) "Warning! Kernel sum is low - consider decreasing cutoff!"
   end if

   call cpu_time(finish)
   
   print '(" Kernel generation took ", f9.6, " s")', finish - start
   
   write (*, *) ""
   write (*, *) ''
   write (*, *) "PROGRAM INFORMATION"
   write (*, *) ""

   if (extendbackward) then
      write (*, *) "Use value at t=0 to pad time before t=0"
   else
      write (*, *) "Use zeros to pad time before t=0"
   end if

   if (extendforward) then
      write (*, *) "Use value at t=final to pad time after t=final"
   else
      write (*, *) "Use zeros to pad time after t=final"
   end if
   write (*, *) "Convolving data..."

   !Time the execution
   call cpu_time(start)
   !actual convolution
   call convolve(signal, kernel, nt, nq)
   call cpu_time(finish)
   write (*, *) "Convolving ended successfully"

   print '(" Convolution took ", f9.6, " s")', finish - start

   call cpu_time(start)
   print '(" Writing data to file",A10)', 'conv_pdW'
   ! Write to file
   open (2, file='conv_pdW')
   do i = 1, nt
      write (2, *) signal(i, :)
   end do

   call cpu_time(finish)
   call cpu_time(tfinish)
   print '(" Writing signal took  ", f9.6, " s")', finish - start
   print '(" Program took  ", f9.6, " s")', tfinish - tstart
contains

   subroutine convolve(signal, kernel, nt, nq)

      double precision, dimension(:) :: kernel
      double precision, dimension(:, :) :: signal
      double precision, dimension(:, :), allocatable :: convsignal
      logical :: extendforward = .false., extendbackward = .false.
      integer :: nt, nq, i, j, k, kerneloffset, kernelsize

      kernelsize = size(kernel)
      !Important to note truncate integer division of odd number!
      !i.e. kerneloffset*2 = kernelsize-1
      kerneloffset = kernelsize/2

      allocate (convsignal(nt + 2*kerneloffset, nq))

      ! embeds signal in the middle of the convsignal array
      convsignal(kerneloffset + 1:kerneloffset + nt, :) = signal

      if (extendbackward) then
         do concurrent (i = 1: kerneloffset)
            convsignal(i, :) = signal(1, :)
         end do
      else
      end if

      if (extendforward) then
         do concurrent (i = 1:kerneloffset)
            convsignal(kerneloffset + nt + i, :) = signal(nt, :)
         end do
      else
      end if

      !loops over q and t. overwrites original signal array.
      do concurrent (i = 1: nq)
         do concurrent (j = 1: nt)
            !Dot product approach - marginally faster
            signal(j, i) = DOT_PRODUCT(kernel, convsignal(j:j + kernelsize, i))

            !BLAS approach - marginally faster, but implementation annoying
            !signal(j,i) = DDOT(kernelsize,kernel,1,convsignal(j:j+kernelsize,i),1)

            !Loop approach - marginally slower
            !signal(j, i) = 0.
            !do k = 1, kernelsize
            !   signal(j, i) = signal(j, i) + kernel(k)*convsignal(k + j - 1, i)
            !end do
         end do
      end do
   end subroutine
   subroutine gaussian(kernel, fwhm, cutoff, dt)

      ! Creates Gaussian kernel with all values larger than cutoff
      ! Odd length

      implicit none

      double precision :: fwhm, cutoff, dt
      double precision :: c, a
      double precision, dimension(:), allocatable :: kernel
      double precision, parameter :: pi = 3.1415926536

      integer :: length, i

      c = 1./(2.*sqrt(2.*log(2.)))*fwhm/dt
      a = 1./(c*sqrt(2.*pi))

      !length = 2*ceiling(sqrt(-log(c*sqrt(2.*pi)*cutoff)*(2.*c**2))) + 1
      length = 2*floor(sqrt(-log(c*sqrt(2.*pi)*cutoff)*(2.*c**2))) + 1

      allocate (kernel(length))
      do i = 1, length
         kernel(i) = a*exp(-((i - ((length + 1)*0.5))**2)/(2*(c**2)))
      end do

      return
   end subroutine

   subroutine kernelnorm(kernel)
      double precision, dimension(:) :: kernel

      kernel = kernel/sum(kernel)
      return
   end subroutine

   subroutine lorentzian(kernel, fwhm, cutoff, dt)
 
      ! Creates lorenztian kernel with all values larger than cutoff
      ! Odd length

      implicit none

      double precision :: fwhm, cutoff, dt
      double precision, dimension(:), allocatable :: kernel
      double precision, parameter :: pi = 3.1415926
      integer :: length, i

      length = 2*floor(0.5*(fwhm/dt)*sqrt(1/(PI*0.5*(fwhm/dt)*cutoff)) - 1) + 1
      allocate (kernel(length))

      do i = 1, length
         kernel(i) = 1/(pi*0.5*(fwhm/dt))*(0.25*(fwhm/dt)**2)/((i - (length + 1)*0.5)**2 + 0.25*(fwhm/dt)**2)
      end do
   end subroutine

   subroutine box(kernel, fwhm, dt)
      
      ! Creates box blur style kernel
      ! Kernel of length 2n+1 with values 1/(2n+1), n integer

           
      implicit none

      double precision :: fwhm, dt, val
      integer :: length
      double precision, dimension(:), allocatable :: kernel

      call make_odd(fwhm, length)

      allocate (kernel(length))
      val = 1/DBLE(length)
      do i = 1, length
         kernel(i) = val
      end do
   end subroutine
   subroutine edge_detection(kernel, fwhm, dt)
      implicit none
      !Messing around with this to see if it is useful - apparently not too much

      double precision :: fwhm, dt, val
      integer :: length
      double precision, dimension(:), allocatable :: kernel
      
      call make_odd(fwhm, length)

      allocate (kernel(length))
      kernel(:)=-1
      kernel(length/2+1) = length+1
   end subroutine

   subroutine make_odd(input, output)
      implicit none
      !takes float and makes it a odd integer

      double precision :: input
      integer :: output
      output = floor(input)
      if (mod(output, 2) .EQ. 0) then
         output = output + 1
      end if
   end subroutine
end program
