#=
Ffowcs Williams-Hawking code
Written by: Benjamin Malczewski & Sam Salehian, Ph.D.
=#

# Load Libraries ---------------------------------------------
using MPI
using Memoization
include("headScript.jl")
include("signalTools.jl")
include("fwh.dat")
include("whatPatch.dat")

# start MPI
MPI.Init()
comm = MPI.COMM_WORLD       # MPI communicator
rank = MPI.Comm_rank(comm)  # rank of local process
sizeMPI = MPI.Comm_size(comm)

# print output header
if rank == 0 && currentPatch==1
    printHeader()
end

tScale = lScale / uScale

# Iteration Vector
iterVec = nStart:nSkip:nEnd
t = iterVec * dt
t = t * tScale
L = length(iterVec)
timeInfo = (iterVec, t, L)

# segment time vector
isOverlap, ieOverlap = makeOverlap(t, nSegments, overlap_percentage, rank)
overlapInfo = (isOverlap, ieOverlap)

if length(PSI) == 1
    onePatch = 1
    nPatches = 1
else
    onePatch = 0
    nPatches = length(solNames)
end

# some tuples for convenience
scales = (lScale, uScale, rhoScale, pScale)
ambientConds = (cinf, pinf, rhoinf)
observerInfo = (xOb, yOb, zOb)
debugInfo = (writeUnitNormals, writeSingleSol, write_dA, writeObserverVectors, observerIndex)

# MAIN LOOP
@memoize computeNoise(rank, comm, sizeMPI, currentPatch, surfDir, gridNames, solNames, onePatch, scales, ambientConds, PSI, overlapInfo, timeInfo, observerInfo, debugInfo)
Memoization.empty_all_caches!()



