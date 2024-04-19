#=
Tool package for signal processing
    includes window functions, fft processes
=#



using FFTW
using LinearAlgebra



function makeWindow(nt, type)
    # make window function based on time vector

    window = Array{Float64,1}(undef, nt)
    if type == 0
        # no window
        window[:] .= 1
        multiplier = 1
    elseif type == 1
        # Hanning window
        for n in 1:nt
            window[n] = 0.5 * (1 - cos(2 * pi * n / (nt - 1)))
            multiplier = 2
        end
    elseif type == 2
        # Hamming window

    end

    return window, multiplier
end


function overlappedFFT(signal, windowType, num_segments, is, ie)

    # Initialize matrix to store FFT results, STFT = Short Time Fourier Transform
    L_seg = ie[1] - is[1] + 1
    stft_vector = Vector{ComplexF64}(undef, L_seg)
    stft_vector[:] .= 0.0 + 0.0im
    sReal = Vector{Float64}(undef, L_seg)
    sImag = Vector{Float64}(undef, L_seg)


    window, windowMultiplier = makeWindow(L_seg, windowType)    # make window once

    # Apply STFT
    for i in 1:num_segments
        start_idx = is[i]
        end_idx = ie[i]

        # Extract segment and apply window function
        segment = signal[start_idx:end_idx]
        windowed_segment = segment .* window

        # Compute FFT
        fft_result = fft(windowed_segment) ./ L_seg
        sReal = Float64.(sign.(real(fft_result)))
        sImag = Float64.(sign.(imag(fft_result)))
        fft_result = (windowMultiplier^2 / num_segments) .* (real(fft_result) .^ 2 + im .* imag(fft_result) .^ 2)        # segment weighted average 

        # Store FFT result
        stft_vector = stft_vector + fft_result
    end
    stft_vector = sqrt.(real(stft_vector)) + im .* sqrt.(imag(stft_vector))
    stft_vector = sReal .* real(stft_vector) + sImag .* im .* imag(stft_vector)
    return stft_vector
end


function getFrequency(t, is, ie)

    fs = 1 / (t[2] - t[1])      # sample frequency

    L = ie[1] - is[1] + 1
    f = f = (fftfreq(L, fs))        # make frequency vector
    f = f[1:Int(ceil(length(f) / 2))]   # only keep half of frequency we care about

    return f
end


function makeOverlap(inVec, nSeg, pOverlap, rank)
    # overlaps segments of vector based on number of segments, and percent overlap between segments

    L = length(inVec)

    is = Vector{Int64}(undef, nSeg)        # starting indices
    ie = Vector{Int64}(undef, nSeg)        # ending indices

    if nSeg > 1
        Leq = floor(L / nSeg)                     # equal length spacing

        # initialize indices
        for i in 1:nSeg
            #typeof(i)
            #typeof(Leq)
            #typeof([i - 1] * Leq + 1)
            is[i] = (i - 1) * Leq + 1
            ie[i] = is[i] + Leq
        end
        if ie[nSeg] > L
            # no overshooting allowed
            ie[nSeg] = ie[nSeg] - 1
            is[nSeg] = is[nSeg] - 1
        end

        # main loop to grow segments
        counter = 1
        overlapActual = 0
        while overlapActual < pOverlap && counter < L + 1
            # grow ends
            ie[1] = ie[1] + 1
            is[nSeg] = is[nSeg] - 1

            # grow middle segments
            for i in 2:nSeg-1
                if ie[i-1] - is[i] > ie[i] - is[i+1]
                    # add to end
                    ie[i] = ie[i] + 1
                else
                    # add to start
                    is[i] = is[i] - 1
                end
            end

            # check while loop condition
            overlap = is[2] - ie[1]
            Lseg = ie[2] - is[2]
            overlapActual = abs(overlap / Lseg) * 100

            counter = counter + 1
        end
    else
        # only 1 segment
        is[1] = 1
        ie[1] = L
    end


    if rank == 0 # print result
        println("Rank=$rank: Segment Information")

        if nSeg > 1
            println("Total Length: $L, Number of Segments: $nSeg, Percent Overlap: $pOverlap%")
            for i in 1:nSeg
                if i > 1
                    overlap = ie[i-1] - is[i]
                else
                    overlap = 0
                end
                println("Segment $i: is=$(is[i]), ie=$(ie[i]), L=$(ie[i]-is[i]), overlap=$overlap")
            end
        else
            println("Total Length: $L, Number of Segments: $nSeg, Percent Overlap: N/A")
            println("is=$(is[1]), ie=$(ie[1])")
        end
    end
    return is, ie
end