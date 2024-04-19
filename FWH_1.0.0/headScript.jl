#=
Head script that does all the work
=#

using MPI
using PencilArrays
using FFTW
using Memoization
include("plot3d_tools.jl")
include("parallel_tools.jl")
include("miscTools.jl")
include("gridMetrics.jl")
include("sources.jl")
include("signalTools.jl")
include("postProcess.jl")
include("fwh.dat")



function computeNoise(rank, comm, sizeMPI, patch, surfDir, gridNames, solNames, onePatch, scales, ambientConds, PSI, overlapInfo, timeInfo, observerInfo, debugInfo)

    # decompose tuples
    lScale, uScale, rhoScale, pScale = scales[1], scales[2], scales[3], scales[4]
    cinf, pinf, rhoinf = ambientConds[1], ambientConds[2], ambientConds[3]
    isOverlap, ieOverlap = overlapInfo[1], overlapInfo[2]
    iterVec, t, L = timeInfo[1], timeInfo[2], timeInfo[3]
    xOb, yOb, zOb = observerInfo[1], observerInfo[2], observerInfo[3]
    writeUnitNormals, writeSingleSol, write_dA, writeObserverVectors, observerIndex = debugInfo[1], debugInfo[2], debugInfo[3], debugInfo[4], debugInfo[5]

    # write grid name and load grid
    if rank == 0
        if onePatch == 1
            gridName = surfDir * gridNames * gridExt
            iLimCurrent = iLim
            jLimCurrent = jLim
            kLimCurrent = kLim
        else
            gridName = surfDir * gridNames[patch] * gridExt
            iLimCurrent = iLim[patch]
            jLimCurrent = jLim[patch]
            kLimCurrent = kLim[patch]
        end
        println("Rank 0: Reading grid: $gridName")
        xGlobe, yGlobe, zGlobe = read_plot3d_grid(gridName, p3dType, [iSkip, jSkip, kSkip], (iLimCurrent, jLimCurrent, kLimCurrent))
        xGlobe, yGlobe, zGlobe = xGlobe .+ gridOffset[1], yGlobe .+ gridOffset[2], zGlobe .+ gridOffset[3]      # apply grid offset
        xGlobe, yGlobe, zGlobe = xGlobe .* lScale, yGlobe .* lScale, zGlobe .* lScale                                 # dimensionalize grid
        nx, ny, nz = size(xGlobe)
        dA_globe = getDA(xGlobe, yGlobe, zGlobe)

        # compute cutoff frequency
        if in(patch, cutoffSurfaces)
            ind_cutoff = findfirst(==(patch), cutoffSurfaces)
            f_cut = computeCutoffFrequency(xGlobe, yGlobe, zGlobe, cinf, xSearchCutoff)
            println("Cutoff Frequency: $f_cut [Hz*PPW]")
            mkpath("output")
        else
            ind_cutoff = 0
        end
        fname = "output/cutoffFrequency"
        if ind_cutoff == 1
            # create new file
            io = open(fname, "w")
            println(io, "Cutoff Frequencies [Hz*PPW]")
            println(io, f_cut)
            close(io)
        elseif ind_cutoff > 1
            # append to file
            io = open(fname, "a")
            println(io, f_cut)
            close(io)
        end



    else
        nx = 0
        ny = 0
        nz = 0
    end
    nx, ny, nz = MPI.bcast((nx, ny, nz), comm, root=0) # Broadcast the values from rank 0 to all processes

    # send the grid to all processors
    if rank != 0
        # initialize x, y, z
        xGlobe = Array{Float64}(undef, nx, ny, nz)
        yGlobe = Array{Float64}(undef, nx, ny, nz)
        zGlobe = Array{Float64}(undef, nx, ny, nz)
        dA_globe = Array{Float64}(undef, nx, ny, nz)
    end
    xGlobe, yGlobe, zGlobe, dA_globe = MPI.bcast((xGlobe, yGlobe, zGlobe, dA_globe), comm, root=0)

    # load data in parallel
    iStart, iEnd = distribute_work(L, sizeMPI, rank)        # distributes for loop indices for each processor
    L_loc = length(iStart:iEnd)

    locSize = iEnd - iStart + 1
    u_loc = Array{Float64}(undef, locSize, nx, ny, nz)
    v_loc = Array{Float64}(undef, locSize, nx, ny, nz)
    w_loc = Array{Float64}(undef, locSize, nx, ny, nz)
    p_loc = Array{Float64}(undef, locSize, nx, ny, nz)
    rho_loc = Array{Float64}(undef, locSize, nx, ny, nz)
    for i in iStart:iEnd
        iLoc = i - iStart + 1   # local i-indice
        if onePatch == 1
            solName = surfDir * solNames * lpad(string(iterVec[i]), numLength, '0') * solExt
            iLimCurrent = iLim
            jLimCurrent = jLim
            kLimCurrent = kLim
        else
            solName = surfDir * solNames[patch] * lpad(string(iterVec[i]), numLength, '0') * solExt
            iLimCurrent = iLim[patch]
            jLimCurrent = jLim[patch]
            kLimCurrent = kLim[patch]
        end
        println("Rank $rank: Reading Solution: $solName")

        u_loc[iLoc, :, :, :], v_loc[iLoc, :, :, :],
        w_loc[iLoc, :, :, :], p_loc[iLoc, :, :, :],
        rho_loc[iLoc, :, :, :] = read_plot3d_sol(solName, p3dType, [iSkip, jSkip, kSkip], (iLimCurrent, jLimCurrent, kLimCurrent))
        Memoization.empty_cache!(read_plot3d_sol)
    end
    u_loc, v_loc, w_loc, rho_loc, p_loc = u_loc .* uScale, v_loc .* uScale, w_loc .* uScale, rho_loc .* rhoScale, p_loc .* pScale .- pinf

    # populate the local range of global variable for rank==0
    if rank == 0
        # initialize global variables
        u0 = Array{Float64}(undef, L, nx, ny, nz)
        v0 = Array{Float64}(undef, L, nx, ny, nz)
        w0 = Array{Float64}(undef, L, nx, ny, nz)
        p0 = Array{Float64}(undef, L, nx, ny, nz)
        rho0 = Array{Float64}(undef, L, nx, ny, nz)

        u0[iStart:iEnd, :, :, :] = u_loc
        v0[iStart:iEnd, :, :, :] = v_loc
        w0[iStart:iEnd, :, :, :] = w_loc
        p0[iStart:iEnd, :, :, :] = p_loc
        rho0[iStart:iEnd, :, :, :] = rho_loc
    end

    # collect all times on zero
    for i in 1:sizeMPI-1    # collect local chunks onto rank=0
        # pass inds
        passToZero(iStart:iEnd, i, comm, rank, 1, "temporal inds")
        inds = zeroRecv(i, comm, rank, 1, "temporal inds")

        # pass u
        passToZero(u_loc, i, comm, rank, 0, "u")
        if rank == 0
            u0[inds, :, :, :] = MPI.recv(comm, source=i, tag=0)    # recieve the message
        end

        # pass v
        passToZero(v_loc, i, comm, rank, 0, "v")
        if rank == 0
            v0[inds, :, :, :] = MPI.recv(comm, source=i, tag=0)    # recieve the message
        end

        # pass w
        passToZero(w_loc, i, comm, rank, 0, "w")
        if rank == 0
            w0[inds, :, :, :] = MPI.recv(comm, source=i, tag=0)    # recieve the message
        end

        # pass p
        passToZero(p_loc, i, comm, rank, 0, "p")
        if rank == 0
            p0[inds, :, :, :] = MPI.recv(comm, source=i, tag=0)    # recieve the message
        end

        # pass rho
        passToZero(rho_loc, i, comm, rank, 0, "rho")
        if rank == 0
            rho0[inds, :, :, :] = MPI.recv(comm, source=i, tag=0)  # recieve the message
        end
    end

    # write solution at first iteration (optional)
    if rank == 0 && writeSingleSol == 1
        println("Writing single grid and solution")
        mkpath("debug")
        if onePatch == 1
            gridSaveName = "debug/" * gridNames * "_grid.x"
            solSaveName = "debug/" * solNames * "_sol.q"
        else
            gridSaveName = "debug/" * gridNames[patch] * "_grid.x"
            solSaveName = "debug/" * solNames[patch] * "_sol.q"
        end
        write_plot3d_grid(xGlobe, yGlobe, zGlobe, gridSaveName, 1)
        write_plot3d_sol(nx, ny, nz, (u0[1, :, :, :], v0[1, :, :, :], w0[1, :, :, :], p0[1, :, :, :], rho0[1, :, :, :]), solSaveName, 1)
    end

    # clear local solutions
    u_loc, v_loc, w_loc, p_loc, rho_loc = nothing, nothing, nothing, nothing, nothing
    GC.gc()

    # get grid normals
    ux, uy, uz = getUnitNormals2(xGlobe, yGlobe, zGlobe, interiorPoint)

    # write unit vectors (optional)
    if rank == 0 && writeUnitNormals == 1
        mkpath("debug")
        if onePatch == 1
            solNameWrite = "debug/" * solNames * "_unitNormals.q"
            gridNameWrite = "debug/" * gridNames * "_grid.x"
        else
            solNameWrite = "debug/" * solNames[patch] * "_unitNormals.q"
            gridNameWrite = "debug/" * gridNames[patch] * "_grid.x"
        end
        dataOut = (ux, uy, uz)
        write_plot3d_sol(nx, ny, nz, dataOut, solNameWrite, 1)
        write_plot3d_grid(xGlobe, yGlobe, zGlobe, gridNameWrite, 1)
    end

    # write dA's (optional)
    if rank == 0 && write_dA == 1
        mkpath("debug")
        if onePatch == 1
            solNameWrite = "debug/" * solNames * "_dA.q"
            gridNameWrite = "debug/" * gridNames * "_grid.x"
        else
            solNameWrite = "debug/" * solNames[patch] * "_dA.q"
            gridNameWrite = "debug/" * gridNames[patch] * "_grid.x"
        end
        dataOut = (dA_globe, dA_globe)
        write_plot3d_sol(nx, ny, nz, dataOut, solNameWrite, 1)
        write_plot3d_grid(xGlobe, yGlobe, zGlobe, gridNameWrite, 1)
    end

    # pencil split
    dimsX, dimsY, dimsZ = pencilSplit(nx, ny, nz, comm)
    println("Rank $rank: xRange: $dimsX, yRange: $dimsY, zRange: $dimsZ")
    x, y, z, dA = xGlobe[dimsX, dimsY, dimsZ], yGlobe[dimsX, dimsY, dimsZ], zGlobe[dimsX, dimsY, dimsZ], dA_globe[dimsX, dimsY, dimsZ] # creates local x,y,z with same spatial split as u,v,w,p,rho

    # pre-allocate bins for local u,v,w,p,rho
    u = Array{Float64}(undef, L, length(dimsX), length(dimsY), length(dimsZ))
    v = Array{Float64}(undef, L, length(dimsX), length(dimsY), length(dimsZ))
    w = Array{Float64}(undef, L, length(dimsX), length(dimsY), length(dimsZ))
    p = Array{Float64}(undef, L, length(dimsX), length(dimsY), length(dimsZ))
    rho = Array{Float64}(undef, L, length(dimsX), length(dimsY), length(dimsZ))

    # distribute grid
    sendDims(dimsX, dimsY, dimsZ, comm, rank)   # rank=0 recieves local dimensions from all other processors
    if rank == 0
        for i in 1:sizeMPI-1
            # recieve local inds from oher rank
            otherInds = MPI.recv(comm, source=i, tag=0)
            #println("Rank $rank: recieved inds: $otherInds")

            # send local chunks out to all other ranks
            MPI.send((u0[:, otherInds[1], otherInds[2], otherInds[3]],
                    v0[:, otherInds[1], otherInds[2], otherInds[3]],
                    w0[:, otherInds[1], otherInds[2], otherInds[3]],
                    p0[:, otherInds[1], otherInds[2], otherInds[3]],
                    rho0[:, otherInds[1], otherInds[2], otherInds[3]]),
                comm, dest=i, tag=0)
            #println("Rank $rank: sent chunk to rank $i")
        end
        u, v, w, p, rho = u0[:, dimsX, dimsY, dimsZ], v0[:, dimsX, dimsY, dimsZ], w0[:, dimsX, dimsY, dimsZ], p0[:, dimsX, dimsY, dimsZ], rho0[:, dimsX, dimsY, dimsZ]
        # clear rank=0's knowledge of full solution to save ram
        u0, v0, w0, p0, rho0 = nothing, nothing, nothing, nothing, nothing
        GC.gc()
    end

    # convert global unit normals to local unit normals
    ux = ux[dimsX, dimsY, dimsZ]
    uy = uy[dimsX, dimsY, dimsZ]
    uz = uz[dimsX, dimsY, dimsZ]

    # each processor recieves its local chunk
    if rank != 0
        # recieve the message
        u, v, w, p, rho = MPI.recv(comm, source=0, tag=0)
        #println("Rank $rank: vars recieved")
    end

    # pre-allocate solution on rank=0
    if rank == 0
        fs = 1 / (t[2] - t[1])              # sampling frequency
        println("Sampling Frequency: $fs Hz")
        f = getFrequency(t, isOverlap, ieOverlap)
        output = Array{ComplexF64}(undef, length(f), length(xOb))

        # write frequency file and observer locations
        if patch == 1

            mkpath("output")
            io1 = open("output/freq", "w")
            for i in eachindex(f)
                println(io1, f[i])
            end
            close(io1)
            io2 = open("output/observers", "w")
            for i in eachindex(xOb)
                println(io2, "$(xOb[i]), $(yOb[i]), $(zOb[i])")
            end
            close(io2)
        end
    else
        # compute frequency vector for everyone else
        f = getFrequency(t, isOverlap, ieOverlap)
    end


    # loop through observers
    for ob in eachindex(xOb)
        if rank == 0
            println("Computing sources for Observer $ob/$(length(xOb))")
        end

        # get vector between observer and patch
        xVec, yVec, zVec, magVec = getVector(x, y, z, xOb[ob], yOb[ob], zOb[ob])

        # write Observer Vectors (debug)
        if rank == 0 && writeObserverVectors == 1 && ob == observerIndex
            mkpath("debug")
            if onePatch == 1
                solNameWrite = "debug/" * solNames * "_obVector.q"
                gridNameWrite = "debug/" * gridNames * "_grid.x"
            else
                solNameWrite = "debug/" * solNames[patch] * "_obVector.q"
                gridNameWrite = "debug/" * gridNames[patch] * "_grid.x"
            end
            xVecTemp, yVecTemp, zVecTemp, magVecTemp = getVector(xGlobe, yGlobe, zGlobe, xOb[ob], yOb[ob], zOb[ob])
            dataOut = (xVecTemp .* magVecTemp, yVecTemp .* magVecTemp, zVecTemp .* magVecTemp)
            write_plot3d_sol(nx, ny, nz, dataOut, solNameWrite, 1)
            write_plot3d_grid(xGlobe, yGlobe, zGlobe, gridNameWrite, 1)
        end

        # calc sources and FFT
        Fhat1, Fhat2, f = calcSources(cinf, xVec, yVec, zVec, magVec, ux, uy, uz, u, v, w, p, rho, f, window_type, nSegments, isOverlap, ieOverlap)
        Fhat1 = Fhat1 .* PSI[patch]
        Fhat2 = Fhat2 .* PSI[patch]

        # double integrate
        integral = integrateSources2(dA, magVec, cinf, rhoinf, f, Fhat1, Fhat2)

        Memoization.empty_cache!(getVector)
        Memoization.empty_cache!(calcSources)
        Memoization.empty_cache!(integrateSources2)


        # sum from all procs
        if rank == 0
            output[:, ob] = integral
        end
        Fhat1, Fhat2 = nothing, nothing
        GC.gc()

        for i in 1:sizeMPI-1
            # pass integral
            passToZero(integral, i, comm, rank, 0, "integral")
            if rank == 0
                output[:, ob] = output[:, ob] + MPI.recv(comm, source=i, tag=0)    # recieve the message
            end
        end
    end

    # write patch solution
    if rank == 0
        if onePatch == 1
            saveName = "output/" * solNames * "_sol"
        else
            saveName = "output/" * solNames[patch] * "_sol"
        end

        writePatchSolution(saveName, output)
    end

    # clear flow solution
    u, v, w, p, rho = nothing, nothing, nothing, nothing, nothing
    GC.gc()

end