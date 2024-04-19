#=
Tool package for handling the FWH source terms:
    1. get vectors between surface and observers
    2. compute F1 and F2 sources, then FFT them
    3. integrate the sources
=#

using LinearAlgebra
using FFTW
using Statistics
using Memoization

include("signalTools.jl")


function getVector(x, y, z, xOb, yOb, zOb)
    # computes an array of vectors between the observer and surface

    S = size(x)
    iMax = S[1]
    jMax = S[2]
    kMax = S[3]

    rMag = Array{Float64,3}(undef, iMax, jMax, kMax)
    rHatX = Array{Float64,3}(undef, iMax, jMax, kMax)
    rHatY = Array{Float64,3}(undef, iMax, jMax, kMax)
    rHatZ = Array{Float64,3}(undef, iMax, jMax, kMax)

    for i in 1:iMax
        for j in 1:jMax
            for k in 1:kMax
                xTemp = xOb - x[i, j, k]
                yTemp = yOb - y[i, j, k]
                zTemp = zOb - z[i, j, k]

                rMag[i, j, k] = sqrt(xTemp^2 + yTemp^2 + zTemp^2)
                rHatX[i, j, k] = xTemp / rMag[i, j, k]
                rHatY[i, j, k] = yTemp / rMag[i, j, k]
                rHatZ[i, j, k] = zTemp / rMag[i, j, k]
            end
        end
    end

    return rHatX, rHatY, rHatZ, rMag
end

function calcSources(cinf, rx, ry, rz, rmag, ux, uy, uz, u, v, w, p, rho, f, window_type, nSeg, is, ie)
    # computes the source terms, F1 and F2 per Mendez et al. 2013
    # then performs FFT and returns the FFT'd quantities

    # get sizes for the loops
    S = size(u)
    nt = S[1]
    nx = S[2]
    ny = S[3]
    nz = S[4]

    # FFT Start
    # window = makeWindow(nt, window_type)           # window function for later

    L_freq = length(f)
    Fhat1 = Array{ComplexF64}(undef, L_freq, nx, ny, nz)
    Fhat2 = Array{ComplexF64}(undef, L_freq, nx, ny, nz)

    # pre-allocate
    F1 = Array{Float64,1}(undef, nt)    # Un from Lyrintzis
    F2 = Array{Float64,1}(undef, nt)    # Lr from Lyrintzis

    # main loop
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                for m in 1:nt
                    # current variable
                    ui = u[m, i, j, k]
                    vi = v[m, i, j, k]
                    wi = w[m, i, j, k]
                    ppi = p[m, i, j, k]
                    #rhoi = rho[m, i, j, k]
                    rhoi = rhoinf+ppi/(340^2)

                    # establish vectors
                    U_vec = [ui, vi, wi]                                # velocity vector
                    R_vec = [rx[i, j, k], ry[i, j, k], rz[i, j, k]]     # vector to observer
                    N_vec = [ux[i, j, k], uy[i, j, k], uz[i, j, k]]     # unit surface normal

                    # compute the terms
                    un = dot(U_vec, N_vec)
                    dp = dot(ppi .* N_vec, R_vec) + dot(rhoi .* U_vec, un .* R_vec)
                    F1[m] = dp / (cinf * rmag[i, j, k]) + rhoi * un / rmag[i, j, k]
                    F2[m] = dp / (rmag[i, j, k]^2)
                end
                F1 = F1 .- mean(F1)
                F2 = F2 .- mean(F2)

                temp1 = overlappedFFT(F1, window_type, nSeg, is, ie)
                temp2 = overlappedFFT(F2, window_type, nSeg, is, ie)
                #Memoization.empty_cache!(overlappedFFT)

                Fhat1[:, i, j, k] = temp1[1:length(f)]
                Fhat2[:, i, j, k] = temp2[1:length(f)]

            end
        end
    end
    temp1, temp2 = nothing, nothing, nothing, nothing
    GC.gc()

    return Fhat1, Fhat2, f
end


# NOT USING #
#=
function integrateSources(x, y, z, magVec, cinf, f, Fhat1, Fhat2)
    # get sizes for the loops
    S = size(x)
    nx = S[1]
    ny = S[2]
    nz = S[3]
    nf = length(f)

    # pre-allocate solution
    phat_out = Array{ComplexF64}(undef, nf)

    for n in 1:nf
        omega = 2 * pi * f[n]
        # make integrand
        if nx == 1
            # integrate over yz
            integrand = Array{ComplexF64}(undef, ny, nz)
            for j in 1:ny-1
                for k in 1:nz-1
                    integrand[j, k] = (im * omega * Fhat1[n, 1, j, k] + Fhat2[n, 1, j, k]) * exp(-im * omega * magVec[1, j, k] / cinf)
                end
            end
        elseif ny == 1
            # integrate over xz
            integrand = Array{ComplexF64}(undef, nx, nz)
            for i in 1:nx
                for k in 1:nz
                    integrand[i, k] = (im * omega * Fhat1[n, i, 1, k] + Fhat2[n, i, 1, k]) * exp(-im * omega * magVec[i, 1, k] / cinf)
                end
            end
        elseif nz == 1
            # integrate over xy
            integrand = Array{ComplexF64}(undef, nx, ny)
            for i in 1:nx
                for j in 1:ny
                    integrand[i, j] = (im * omega * Fhat1[n, i, j, 1] + Fhat2[n, i, j, 1]) * exp(-im * omega * magVec[i, j, 1] / cinf)
                end
            end
        end

        # integrate
        if nx == 1
            xp = permutedims(y, (2, 3, 1))
            yp = permutedims(z, (2, 3, 1))
        elseif ny == 1
            xp = permutedims(x, (1, 3, 2))
            yp = permutedims(z, (1, 3, 2))
        elseif nz == 1
            xp = x
            yp = y
        end

        S2 = size(xp)
        nx2 = S2[1]
        ny2 = S2[2]

        answer = 0 + 0im
        for i in 1:nx2-1
            for j in 1:ny2-1
                dx = xp[i+1, j, 1] - xp[i, j, 1]
                dy = yp[i, j+1, 1] - yp[i, j, 1]
                integral = dx * dy * (integrand[i, j] + integrand[i+1, j] + integrand[i, j+1] + integrand[i+1, j+1]) / 4
                answer = answer + integral
            end
        end

        phat_out[n] = answer

    end

    return phat_out
end
=#

function integrateSources2(dA, magVec, cinf, rhoinf, f, F1, F2)
    # get sizes for the loops
    S = size(dA)
    nx = S[1]
    ny = S[2]
    nz = S[3]
    nf = length(f)

    # pre-allocate solution
    phat_out = Array{ComplexF64}(undef, nf)
    integrand1 = Array{ComplexF64}(undef, nx, ny, nz)
    integrand2 = Array{ComplexF64}(undef, nx, ny, nz)

    for n in 1:nf
        omega = 2 * pi * f[n]
        for i in 1:nx
            for j in 1:ny
                for k in 1:nz
                    integrand1[i, j, k] = im * omega * F1[n, i, j, k] * exp(-im * omega * magVec[i, j, k] / cinf)
                    integrand2[i, j, k] = F2[n, i, j, k] * exp(-im * omega * magVec[i, j, k] / cinf)
                end
            end
        end
        integrand1 = integrand1 .* dA
        integrand2 = integrand2 .* dA

        I1 = sum(integrand1)
        I2 = sum(integrand2)

        phat_out[n] = (I1 + I2) / (4 * pi)
    end

    return phat_out
end
