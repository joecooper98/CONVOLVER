#!/usr/local/bin/julia

using Profile
using MAT
using DelimitedFiles

trajfile = "pdW"

# will work for any local (non-periodic) function

fwhm = 100. # in the same units as the time signal

cutoff = 0.000001 # the smallest value of the kernel

#imptrajfile = matopen(trajfile)
pdW = readdlm(trajfile, ' ', Float64)
#int= read(imptrajfile,"pdW")
#close(imptrajfile)

#pdW=permutedims(int[2:end,:]) # delete first row and take transpose for speed

vectorfile = "Vectors.mat"

impvecfile = matopen(vectorfile)

q = read(impvecfile, "qAng")
t = read(impvecfile, "tt")

# calculate dt, which is contant

const dt = t[2]-t[1]

#function

function create_kernel(fwhm::Float64, dt::Float64, cutoff::Real, type::String)
        type = lowercase(type)
                if type == "lorentzian" 
                         println("Using Lorentzian for convolution...")
                         return lorentzian(fwhm,dt,cutoff)
                 #  elseif typeoffunction  == "step"
                         #  f = step
                         #  println("Using step function")
                 elseif type  == "step"
                         println("Using step function for convolution...")
                         println("The width of the function is $fwhm...")
                         println("The time step is $dt...")
                         return box(fwhm,kernelwidth,dt)

                 else
                         type  = "gaussian"
                         println("Using Gaussian for convolution...")
                         println("The FWHM of the function is $fwhm...")
                         println("The time step is $dt...")
                         return gaussian(fwhm,dt,cutoff)
                 end
        end

function lorentzian(fwhm::Float64, dt::Float64, cutoff::Float64) # define a normalised lorentzian as defined on wikipedia
                         fwhm /= dt
                         gam=0.5*fwhm
                         length = 2*(ceil(gam*√(1/(pi*gam*cutoff)) - 1)) + 1
                         return [1/(pi*gam) * (gam^2)/((x-(length+1)*0.5)^2+gam^2) for x in 1:length]
                 end
        
function gaussian(fwhm::Float64, dt::Float64, cutoff::Float64)::Array{Float64} # return a normalised gaussian vector of a certain length as defined on wikipedia
               fwhm /= dt
               c = 1/(2*√(2log(2)))*fwhm
               a = 1/(c*√(2pi))
               length =2*convert(Int64,ceil(√(-log(c*√(2pi)*cutoff)*(2*c^2)))) + 1
               return  [ a*exp(-(x-(length+1)*0.5)^2/(2*c^2)) for x in 1:length]
               end




function convolve(signal::Array{Float64}, kernel::Array{Float64}, extendforward::Bool=true, extendbackward::Bool=false) # fft convolution for testing
        
        # size of the arrays
        nq = size(signal,2)
        nt = size(signal,1)
        kernelsize = size(kernel,1)

        # find embedding indices
        # we only extend by half minus one the length of the kernel, as we only care about the indices inside the original data region
        cutoffind::Int64= convert(Int64,(kernelsize-1)/2)
        lastindex::Int64= cutoffind + nt
        totalsize::Int64= cutoffind + lastindex
        kep = kernel[end]
        ksp = kernel[begin]
        totalcalcs = nt*nq*kernelsize
        println("Convolving with kernel of size ",kernelsize,"...")
        println("Edge values = ",ksp,", ",kep,"...")
        println("Convolving ",nq," q points and ",nt," time points with ",totalcalcs," total calculations...")

        #pad signal - we only need pad by the total size of the kernel, as we only care about the data in the same region as the original scattering
         
        padsignal = zeros(Float64,totalsize,nq) # initiate empty array

        padsignal[cutoffind+1:end-cutoffind,:] = signal # embed signal in the middle

        if extendforward 
                println("Extending the array forward in time using the last point for each q...")
                for j in 1:nq, i in lastindex:totalsize # copy the last value of the signal forward in time 
                        padsignal[i,j] = padsignal[lastindex-1,j]
                end 
        else
                println("Extending the array forward in time using zeros...")
        end

        if extendbackward
                println("Extending the array backward in time using the last point for each q...")
                for j in 1:nq, i in 1:cutoffind # copy the last value of the signal forward in time 
                        padsignal[i,j] = padsignal[cutoffind+1,j]
                end 
        else
                println("Extending the array backward in time using zeros...")
        end

        # initiate kernel

        # identity kernel for checking 
        #  kernel = zeros(Float64,kernelsize,1)
        #  kernel[cutoffind]=1.
        
       # signal = zeros(Float64,nt,nq) # initiate the convolved signal
        
        # here we only scan the values for which we have the full kernel.
        # The kernel starts with its centre at the first value of the original signal (left side in the padding), and then scans along, performing a sum(hadamard(kernel,signal) for each time
        # This is then repeated for every q. A triple nested loop may seem slow, but testing is actually fairly fast. An alternative method for large series would be using the Convolution Theorem, but this is harder to implement
        print("Total time taken for convolution  =")
        @time for i in 1:nq #nq # for all q
                   for j in 1:nt # for all t
                        signal[j,i] = 0.
                        for k in 1:kernelsize
                              signal[j,i] += kernel[k] * padsignal[k+j-1,i]
                      end
                end
        end
 return signal
end

kernel = create_kernel(fwhm, dt, cutoff, "Gaussian")

convolved_signal = convolve(pdW,kernel)

open("convolved_signal", "w") do io
        writedlm(io,convolved_signal)
end
