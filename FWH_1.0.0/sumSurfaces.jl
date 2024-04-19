#=
Script to sum the surfaces after a run
=#

include("postProcess.jl")
include("miscTools.jl")
include("fwh.dat")

# print header
printCleanupHeader()

# determine if there's one patch or not
if length(PSI)==1
    onePatch = 1
else
    onePatch = 0
end

# read all patch solutions, sum them, and write a file with the absolute value
println("Writing absolute value of summed surfaces")
writeName = "output/" * summedName
sumSurfacesAndWriteAbs("output/", solNames, "_sol", onePatch, writeName)
