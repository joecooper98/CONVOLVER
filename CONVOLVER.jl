#!/usr/local/bin/julia

using Profile
using MAT
using DelimitedFiles

trajfile = "pdW"

fwhm = 10. # in the same units as the time signal

pad = 4 # the number of std. dev. around the gaussian which are calculated. Lowering this is a significant speed up at the cost of less accuracy. 4-6 is a good starting value.

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

function convolve(signal::Array{Float64}, fwhm::Float64, lor::Bool,dt ,cutoff) # fft convolution for testing
        
        # converts from fwhm to sigma for gaussian

        fwhmconv=1/(2*√(2*log(2)))
        

        # size of the arrays
        nq = size(signal,2)
        nt = size(signal,1)


        function gaussian(x::Float64,b::Float64,fwhm::Float64) # define a normalised gaussian as defined on wikipedia
                c = fwhmconv*fwhm
                a = 1/(c*√(2pi))
                return a * exp(-(x-b)^2/(2*c^2))
                end
        
        function lorentzian(x::Float64,b::Float64,fwhm::Float64) # define a normalised lorentzian as defined on wikipedia
                gam=0.5*fwhm
                return 1/(pi*gam) * ((gam^2)/((x-b)^2+gam^2))
                end

        # find embedding indices
        #
        cutoffind = convert(Int64,ceil(cutoff*fwhmconv*fwhm/dt)) # find the length of the total vector - this is a mathematical convenience. By cutting of the vector at (e.g.) 6 std.devs., we get 99% of the kernel at massively reduced comp cost
        lastindex = cutoffind + nt
        totalsize = cutoffind + lastindex        


        #pad signal - we only need pad by the total size of the kernel, as we only care about the data in the same region as the original scattering
         
        padsignal = zeros(Float64,totalsize,nq) # initiate empty array

        padsignal[cutoffind+1:end-cutoffind,:] = signal # embed signal in the middle

        for j in 1:nq, i in lastindex:totalsize # copy the last value of the signal forward in time 
                padsignal[i,j] = padsignal[lastindex-1,j]
        end 
        
        # initiate kernel
        kernel = [gaussian(convert(Float64,x),convert(Float64,cutoffind),fwhm/0.5) for x in 0:2cutoffind] # only gaussian atm. This creates an array of 2 * cuttof * sigma + 1, with gaussian centred . The edge values should be small for good results

        if kernel[begin] > 1e-7 || kernel[end] > 1e-7 # warning function
                println("Warning! Kernel edge values are bigger than 1x10^-7 - this may impact accuracy")
        end

        kernelsize = size(kernel,1) # size of kernel for looping over
       
        # identity kernel for checking 
        #  kernel = zeros(Float64,kernelsize,1)
        #  kernel[cutoffind]=1.
        

        convolved_signal = zeros(Float64,nt,nq) # initiate the convolved signal
        
        # here we only scan the values for which we have the full kernel.
        # The kernel starts with its centre at the first value of the original signal (left side in the padding), and then scans along, performing a sum(hadamard(kernel,signal) for each time
        # This is then repeated for every q. A triple nested loop may seem slow, but testing is actually fairly fast. An alternative method for large series would be using the Convolution Theorem, but this is harder to implement

        for i in 1:nq #nq # for all q
                for j in 1:nt # for all t
                      for k in 1:kernelsize
                              convolved_signal[j,i] += kernel[k] * padsignal[k+j-1,i]
                      end
                end
        end
 return convolved_signal
end

convolved_signal = convolve(pdW, fwhm, false, dt, pad)

open("convolved_signal", "w") do io
        writedlm(io,convolved_signal)
end

