#=
Some tools that make parallelization simpler
=#

using MPI
using PencilArrays


function distribute_work(total_iterations, num_processes, rank)
    # distributes a load in single dimension
    # handy for parallel for loops
    chunk_size, remainder = divrem(total_iterations, num_processes)
    start_idx = rank * chunk_size + min(rank, remainder) + 1
    end_idx = start_idx + chunk_size - 1 + (rank < remainder)
    return start_idx, end_idx
end


function pencilSplit(nx, ny, nz, comm)
    if nz == 1
        # permute to put nz as first dim
        dims_global = (nz, ny, nx)  # we know the temporal dimension so no need to include it here
        pen = Pencil(dims_global, comm)
        dummy = PencilArray{Float64}(undef, pen)
        dims_local_t = range_local(dummy)                                   # temporary local dims in wrong order
        dims_local = (dims_local_t[3], dims_local_t[2], dims_local_t[1])    # put dims back in right order
    elseif ny == 1
        # permute to put ny as first dim
        dims_global = (ny, nx, nz)  # we know the temporal dimension so no need to include it here
        pen = Pencil(dims_global, comm)
        dummy = PencilArray{Float64}(undef, pen)
        dims_local_t = range_local(dummy)                                   # temporary local dims in wrong order
        dims_local = (dims_local_t[2], dims_local_t[1], dims_local_t[3])    # put dims back in right order
    else
        # proceed as normal
        dims_global = (nx, ny, nz)  # we know the temporal dimension so no need to include it here
        pen = Pencil(dims_global, comm)
        dummy = PencilArray{Float64}(undef, pen)
        dims_local = range_local(dummy)    # local dimensions for splitting the global variables
    end
    dimsX = dims_local[1]
    dimsY = dims_local[2]
    dimsZ = dims_local[3]

    return dimsX, dimsY, dimsZ
end


function sendDims(dimsX, dimsY, dimsZ, comm, rank)
    if rank != 0
        # send local indices to rank=0
        MPI.send((dimsX, dimsY, dimsZ), comm, dest=0, tag=0)
    end
end

function passToZero(message, i, comm, rank, t, memo)
    # function intended to be run in a for loop where all rank!=0 processors send to rank=0
    if i == rank
        # send message
        MPI.send(message, comm, dest=0, tag=t)
        #println("Rank $rank: Message Sent, $memo")
    end
end

function zeroRecv(sender, comm, rank, t, memo)
    if rank == 0
        # recieve the message
        message = MPI.recv(comm, source=sender, tag=t)
        #println("Rank $rank: Message recieved from rank $sender, $memo")
        return message
    else
        return nothing
    end
end