#=
Function package used to post process FWH solution
=#

using DelimitedFiles


function writePatchSolution(saveName, data)
    nf, nOb = size(data)

    println("Writing  patch output: $saveName")
    io = open(saveName, "w")
    for i in 1:nf
        for j in 1:nOb
            if j == nOb
                if imag(data[i, j]) < 0
                    println(io, "$(real(data[i,j]))$(imag(data[i,j]))j")
                elseif imag(data[i, j]) > 0
                    println(io, "$(real(data[i,j]))+$(imag(data[i,j]))j")
                else
                    println(io, "$(real(data[i,j]))+0.0j")
                end
            else
                if imag(data[i, j]) < 0
                    print(io, "$(real(data[i,j]))$(imag(data[i,j]))j, ")
                elseif imag(data[i, j]) > 0
                    print(io, "$(real(data[i,j]))+$(imag(data[i,j]))j, ")
                else
                    print(io, "$(real(data[i,j]))+0.0j, ")
                end
            end
        end
    end
    close(io)
end




function sumSurfacesAndWriteAbs(folder, solNames, solExt, onePatch, writeName)
    # adds real and imaginary components from each surface and writes a similarly formatted file that has the absolute value

    if onePatch == 1
        readName = folder * solNames * solExt
        println("Reading File: $readName")
        data = readdlm(readName)

        nf, nOb = size(data)
        data = parseComplex.(data)
        absVal = abs.(data)
    else
        nPatches = length(solNames)

        # read first patch solution to get size
        readName = folder * solNames[1] * solExt
        println("Reading File: $readName")
        data = readdlm(readName)
        nf, nOb = size(data)
        dataOut = parseComplex.(data)

        # loop through surfaces
        for i in 2:nPatches
            readName = folder * solNames[i] * solExt
            println("Reading File: $readName")
            data = readdlm(readName)
            dataOut = dataOut + parseComplex.(data)
        end

        # convert to absolute value
        absVal = abs.(dataOut)
    end


    # write final solution
    println("Writing Solution: $writeName")
    io = open(writeName, "w")
    for i in 1:nf
        for j in 1:nOb
            if j == nOb
                # go to next line
                println(io, "$(absVal[i,j])")
            else
                # add comma
                print(io, "$(absVal[i,j]), ")
            end
        end
    end

end

function parseComplex(s)
    return parse(ComplexF64, replace(s, "j" => "im"))
end