#=
Tool package for handling plot3d files
Functions include:
1. read_plot3d_grid
2. read_plot3d_sol
3. write_plot3d_grid
4. write_plot3d_sol
=#
#------------------------------------------------------------------------------------

# 1. read_plot3d_grid
function read_plot3d_grid(filename::AbstractString, type, skip=[1, 1, 1], limits=([1, -1], [1, -1], [1, -1]))
    # Open the Plot3D file
    # type=1, formatted
    # type=2, unformatted


    # parse the skipping here
    iSkip = skip[1]
    jSkip = skip[2]
    kSkip = skip[3]

    #----------------------------------------------------------------------------------------------------
    if type == 1
        # formatted file
        try
            open(filename, "r") do file
                # Read the number of blocks (not used in this function)

                #  Read grid dimensions (nx, ny, nz)
                readline(file)
                dimensions = split(readline(file))
                nx, ny, nz = parse.(Int, dimensions)
                println("nx=", nx, " ny=", ny, " nz=", nz)

                iStart, iEnd, jStart, jEnd, kStart, kEnd = parseLims(limits, nx, ny, nz)

                x = reshape(parse.(Float64, split(readline(file))), (nx, ny, nz))
                y = reshape(parse.(Float64, split(readline(file))), (nx, ny, nz))
                z = reshape(parse.(Float64, split(readline(file))), (nx, ny, nz))

                return x[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       y[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       z[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd]
            end
        catch e
            println("Error reading file: $e")
            return nothing, nothing, nothing
        end
        #------------------------------------------------------------------------------
    elseif type == 2
        # unformatted binary file, single precision
        try
            open(filename, "r") do file
                # Read the number of blocks (not used in this function)
                junk, nx, ny, nz = read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32)
                #junk2, junk3 = read(file,Int32), read(file,Int32)

                #println(junk, junk2, junk3)
                println("nx=", nx, " ny=", ny, " nz=", nz)
                x = Array{Float64,3}(undef, nx, ny, nz)
                y = Array{Float64,3}(undef, nx, ny, nz)
                z = Array{Float64,3}(undef, nx, ny, nz)
                junk, junk = read(file, Float32), read(file, Float32)

                iStart, iEnd, jStart, jEnd, kStart, kEnd = parseLims(limits, nx, ny, nz)

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            x[i, j, k] = Float64(read(file, Float32))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            y[i, j, k] = Float64(read(file, Float32))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            z[i, j, k] = Float64(read(file, Float32))
                        end
                    end
                end

                # Close the file
                close(file)
                iSkip = skip[1]
                jSkip = skip[2]
                kSkip = skip[3]

                return x[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       y[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       z[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd]
            end
        catch e
            println("Error reading file: $e")
            return nothing, nothing, nothing
        end


    elseif type == 3
        # unformatted binary file, double precision
        try
            open(filename, "r") do file
                # Read the number of blocks (not used in this function)
                println("file opened: $filename")
                junk, junk2, junk3, junk4, nx, ny, nz = read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32)

                iStart, iEnd, jStart, jEnd, kStart, kEnd = parseLims(limits, nx, ny, nz)

                println("nx=", nx, " ny=", ny, " nz=", nz)
                x = Array{Float64,3}(undef, nx, ny, nz)
                y = Array{Float64,3}(undef, nx, ny, nz)
                z = Array{Float64,3}(undef, nx, ny, nz)
                junk, junk = read(file, Float32), read(file, Float32)
                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            x[i, j, k] = Float64(read(file, Float64))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            y[i, j, k] = Float64(read(file, Float64))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            z[i, j, k] = Float64(read(file, Float64))
                        end
                    end
                end

                # Close the file
                close(file)

                return x[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd],
                       y[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd],
                       z[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd]
            end
        catch e
            println("Error reading file: $e")
            return nothing, nothing, nothing
        end
    end
end


#-------------------------------------------------------------------------

# 2. read_plot3d_sol
function read_plot3d_sol(filename::AbstractString, type, skip=[1, 1, 1], limits=([1, -1], [1, -1], [1, -1]))
    # Open the Plot3D file
    # type=1, formatted
    # type=2, unformatted

    # parse the skipping here
    iSkip = skip[1]
    jSkip = skip[2]
    kSkip = skip[3]

    #----------------------------------------------------------------------------------------------------
    if type == 1
        # formatted file, doesnt work yet
        try
            open(filename, "r") do file
                # Read the number of blocks (not used in this function)

                #  Read grid dimensions (nx, ny, nz)
                readline(file)
                dimensions = split(readline(file))
                nx, ny, nz = parse.(Int, dimensions)
                println("nx=", nx)
                println("ny=", ny)
                println("nz=", nz)

                x = reshape(parse.(Float64, split(readline(file))), (nx, ny, nz))
                y = reshape(parse.(Float64, split(readline(file))), (nx, ny, nz))
                z = reshape(parse.(Float64, split(readline(file))), (nx, ny, nz))

                return x, y, z
            end
        catch e
            println("Error reading file: $e")
            return nothing, nothing, nothing
        end
        #------------------------------------------------------------------------------
    elseif type == 2
        # unformatted binary file, single precision
        try
            open(filename, "r") do file
                # Read the number of blocks (not used in this function)
                junk, nx, ny, nz, nvar = read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32)
                nStep, tau, junk2 = read(file, Int32), read(file, Float64), read(file, Int64)

                iStart, iEnd, jStart, jEnd, kStart, kEnd = parseLims(limits, nx, ny, nz)

                #println("nx=$nx, ny=$ny, nz=$nz, nvar=$nvar, nStep=$nStep, tau=$tau")
                u = Array{Float64}(undef, nx, ny, nz)
                v = Array{Float64}(undef, nx, ny, nz)
                w = Array{Float64}(undef, nx, ny, nz)
                p = Array{Float64}(undef, nx, ny, nz)
                rho = Array{Float64}(undef, nx, ny, nz)

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            u[i, j, k] = Float64(read(file, Float32))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            v[i, j, k] = Float64(read(file, Float32))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            w[i, j, k] = Float64(read(file, Float32))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            p[i, j, k] = Float64(read(file, Float32))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            rho[i, j, k] = Float64(read(file, Float32))
                        end
                    end
                end

                # Close the file
                close(file)

                iSkip = skip[1]
                jSkip = skip[2]
                kSkip = skip[3]

                return u[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       v[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       w[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       p[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       rho[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd]
            end
        catch e
            println("Error reading file: $e")
            return nothing, nothing, nothing, nothing, nothing
        end

    elseif type == 3
        # unformatted binary file, double precision
        try
            open(filename, "r") do file
                # Read the number of blocks (not used in this function)
                junk, jun2, jun3, junk4, nx, ny, nz, nvar = read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32), read(file, Int32)
                junk20 = read(file, Int64)

                iStart, iEnd, jStart, jEnd, kStart, kEnd = parseLims(limits, nx, ny, nz)

                #println("nx=$nx, ny=$ny, nz=$nz, nvar=$nvar")
                u = Array{Float64}(undef, nx, ny, nz)
                v = Array{Float64}(undef, nx, ny, nz)
                w = Array{Float64}(undef, nx, ny, nz)
                p = Array{Float64}(undef, nx, ny, nz)
                rho = Array{Float64}(undef, nx, ny, nz)

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            u[i, j, k] = Float64(read(file, Float64))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            v[i, j, k] = Float64(read(file, Float64))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            w[i, j, k] = Float64(read(file, Float64))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            p[i, j, k] = Float64(read(file, Float64))
                        end
                    end
                end

                for k in 1:nz
                    for j in 1:ny
                        for i in 1:nx
                            rho[i, j, k] = Float64(read(file, Float64))
                        end
                    end
                end

                # Close the file
                close(file)

                iSkip = skip[1]
                jSkip = skip[2]
                kSkip = skip[3]

                return u[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       v[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       w[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       p[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd], 
                       rho[iStart:iSkip:iEnd,jStart:jSkip:jEnd,kStart:kSkip:kEnd]
            end
        catch e
            println("Error reading file: $e")
            return nothing, nothing, nothing, nothing, nothing
        end
    end
end


#---------------------------------------------------------------------------
function write_plot3d_grid(x, y, z, saveName, type)
    if type == 1
        # formatted plot3d file
        println("Writing grid file: $saveName")
        nx, ny, nz = size(x)
        io = open(saveName, "w")
        println(io, 1)
        println(io, "$nx $ny $nz")

        # write x
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    print(io, x[i, j, k], " ")
                end
            end
        end
        print(io, "\n")
        # write y
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    print(io, y[i, j, k], " ")
                end
            end
        end
        print(io, "\n")
        # write z
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    print(io, z[i, j, k], " ")
                end
            end
        end
        close(io)

    elseif type == 2
        # unformatted plot3d file (doesn't work)


    end
end

#---------------------------------------------------------------------------
function write_plot3d_sol(nx, ny, nz, data, saveName, type)
    if type == 1
        # formatted p3d file
        println("Writing solution file: $saveName")
        nData = length(data)
        dataSave = Array{Float64}(undef, nData, nx, ny, nz)
        io = open(saveName, "w")
        println(io, 1)
        println(io, "$nx $ny $nz $nData")
        for N in 1:nData
            q = data[N]
            #qflat = permutedims(q, [1,2,3])[:]
            #println(io,qflat[:])
            for k in 1:nz
                for j in 1:ny
                    for i in 1:nx
                        print(io, q[i, j, k], " ")
                    end
                end
            end
            if N != nData
                print(io, "\n")
            end
        end

        close(io)


    elseif type == 2
        # unformatted p3d file
        nData = length(data)
        println("Writing solution file: $saveName")
        dataSave = Array{Float64}(undef, nData, nx, ny, nz)
        for N in 1:nData
            dataSave[N, :, :, :] = data[N]
        end
        flattened_q = permutedims(dataSave, [2, 3, 4, 1])[:]
        io = open(saveName, "w")
        write(io, Int32(0), Int32(nx), Int32(ny), Int32(nz), Int32(nData), Int32(50), Float64(1.0), Int64(1))
        for N in 1:nData
            q = data[N]
            println(size(q))
            qflat = permutedims(q, [1, 2, 3])[:]
            println(size(qflat))
            write(io, convert(Array{Float32,1}, qflat))
        end
        #write(io, flattened_q)
        close(io)
    end
end

#---------------------------------------------------------------------------
function parseLims(lims, nx, ny, nz)
    # parses the limits to clean up the script in the read functions

    # i
    iLims = lims[1]
    iStart = iLims[1]
    iEnd = iLims[2]
    if iLims[2] == -1
        # overwrite iEnd
        iEnd = nx
    end

    # j
    jLims = lims[2]
    jStart = jLims[1]
    jEnd = jLims[2]
    if jLims[2] == -1
        # overwrite iEnd
        jEnd = ny
    end

    # k
    kLims = lims[3]
    kStart = kLims[1]
    kEnd = kLims[2]
    if kLims[2] == -1
        # overwrite iEnd
        kEnd = nz
    end

    return iStart, iEnd, jStart, jEnd, kStart, kEnd
end