#=
Tool package for getting grid metrics for a 2D grid in 3d space.  Specifically:
    1. dA's
    2. unit normal vectors
=#

using LinearAlgebra
using Statistics

function getDA(x, y, z)
    # computes differential areas at each node
    # areas are based on area of triangle between node and midpoint of two other nodes
    # there definitely is a more efficient way to write this, but since its only done once, it's okay
    # one of the dimensions in x,y,z needs to be 1

    nx, ny, nz = size(x)
    nx_orig, ny_orig, nz_orig = nx, ny, nz

    if nx == 1
        # permute to be yz
        x = permutedims(x, (2, 3, 1))
        y = permutedims(y, (2, 3, 1))
        z = permutedims(z, (2, 3, 1))

        xt = x
        yt = y
        zt = z

        x = yt
        y = zt
        z = xt
    elseif ny == 1
        # permute to be xz
        x = permutedims(x, (1, 3, 2))
        y = permutedims(y, (1, 3, 2))
        z = permutedims(z, (1, 3, 2))

        xt = x
        yt = y
        zt = z

        x = xt
        y = zt
        z = yt
    end

    # collapse z dimension
    x = x[:, :, 1]
    y = y[:, :, 1]
    z = z[:, :, 1]

    nx, ny = size(x)

    dA = Array{Float64,2}(undef, nx, ny)
    for i in 1:nx
        for j in 1:ny
            A = 0
            p1 = [x[i, j], y[i, j], z[i, j]]
            if i == 1 && j == 1
                # bottom left corner
                p2t = [x[i+1, j], y[i+1, j], z[i+1, j]]
                p3t = [x[i, j+1], y[i, j+1], z[i, j+1]]
                p2m = getMidpoint(p1, p2t)
                p3m = getMidpoint(p1, p3t)
                v1 = p2m - p1
                v2 = p3m - p1

                cp = cross(v1, v2)
                A = norm(cp)
            elseif i == nx && j == 1
                # bottom right corner
                p2t = [x[i-1, j], y[i-1, j], z[i-1, j]]
                p3t = [x[i, j+1], y[i, j+1], z[i, j+1]]
                p2m = getMidpoint(p1, p2t)
                p3m = getMidpoint(p1, p3t)
                v1 = p2m - p1
                v2 = p3m - p1

                cp = cross(v1, v2)
                A = norm(cp)
            elseif i == 1 && j == ny
                # top left corner
                p2t = [x[i+1, j], y[i+1, j], z[i+1, j]]
                p3t = [x[i, j-1], y[i, j-1], z[i, j-1]]
                p2m = getMidpoint(p1, p2t)
                p3m = getMidpoint(p1, p3t)
                v1 = p2m - p1
                v2 = p3m - p1

                cp = cross(v1, v2)
                A = norm(cp)
            elseif i == nx && j == ny
                # top right corner
                p2t = [x[i-1, j], y[i-1, j], z[i-1, j]]
                p3t = [x[i, j-1], y[i, j-1], z[i, j-1]]
                p2m = getMidpoint(p1, p2t)
                p3m = getMidpoint(p1, p3t)
                v1 = p2m - p1
                v2 = p3m - p1

                cp = cross(v1, v2)
                A = norm(cp)
            elseif i == 1 && j != 1 && j != ny
                # left wall
                p2t = [x[i+1, j], y[i+1, j], z[i+1, j]]
                p3t = [x[i, j-1], y[i, j-1], z[i, j-1]]
                p4t = [x[i, j+1], y[i, j+1], z[i, j+1]]
                p2m = getMidpoint(p1, p2t)
                p3m = getMidpoint(p1, p3t)
                p4m = getMidpoint(p1, p4t)
                v2 = p2m - p1
                v3 = p3m - p1
                v4 = p4m - p1

                cp1 = cross(v2, v4)
                cp2 = cross(v2, v3)
                A = norm(cp1) + norm(cp2)
            elseif i == nx && j != 1 && j != ny
                # right boundary
                p2t = [x[i-1, j], y[i-1, j], z[i-1, j]]
                p3t = [x[i, j-1], y[i, j-1], z[i, j-1]]
                p4t = [x[i, j+1], y[i, j+1], z[i, j+1]]
                p2m = getMidpoint(p1, p2t)
                p3m = getMidpoint(p1, p3t)
                p4m = getMidpoint(p1, p4t)
                v2 = p2m - p1
                v3 = p3m - p1
                v4 = p4m - p1

                cp1 = cross(v2, v3)
                cp2 = cross(v2, v4)
                A = norm(cp1) + norm(cp2)
            elseif j == 1 && i != 1 && i != nx
                # bottom boundary
                p2t = [x[i-1, j], y[i-1, j], z[i-1, j]]
                p3t = [x[i+1, j], y[i+1, j], z[i+1, j]]
                p4t = [x[i, j+1], y[i, j+1], z[i, j+1]]
                p2m = getMidpoint(p1, p2t)
                p3m = getMidpoint(p1, p3t)
                p4m = getMidpoint(p1, p4t)
                v2 = p2m - p1
                v3 = p3m - p1
                v4 = p4m - p1

                cp1 = cross(v4, v2)
                cp2 = cross(v4, v3)
                A = norm(cp1) + norm(cp2)
            elseif j == ny && i != 1 && i != nx
                # top boundary
                p2t = [x[i-1, j], y[i-1, j], z[i-1, j]]
                p3t = [x[i+1, j], y[i+1, j], z[i+1, j]]
                p4t = [x[i, j-1], y[i, j-1], z[i, j-1]]
                p2m = getMidpoint(p1, p2t)
                p3m = getMidpoint(p1, p3t)
                p4m = getMidpoint(p1, p4t)
                v2 = p2m - p1
                v3 = p3m - p1
                v4 = p4m - p1

                cp1 = cross(v4, v2)
                cp2 = cross(v4, v3)
                A = norm(cp1) + norm(cp2)
            else
                # interior point
                p2t = [x[i-1, j], y[i-1, j], z[i-1, j]]
                p3t = [x[i+1, j], y[i+1, j], z[i+1, j]]
                p4t = [x[i, j-1], y[i, j-1], z[i, j-1]]
                p5t = [x[i, j+1], y[i, j+1], z[i, j+1]]
                p2m = getMidpoint(p1, p2t)
                p3m = getMidpoint(p1, p3t)
                p4m = getMidpoint(p1, p4t)
                p5m = getMidpoint(p1, p5t)
                v2 = p2m - p1
                v3 = p3m - p1
                v4 = p4m - p1
                v5 = p5m - p1

                cp1 = cross(v2, v4)
                cp2 = cross(v2, v5)
                cp3 = cross(v3, v4)
                cp4 = cross(v3, v5)
                A = norm(cp1) + norm(cp2) + norm(cp3) + norm(cp4)
            end
            dA[i, j] = A
        end
    end

    # put dA back into 3D form
    dA_out = Array{Float64,3}(undef, nx_orig, ny_orig, nz_orig)
    if nx_orig == 1
        dA_out[1, :, :] = dA
    elseif ny_orig == 1
        dA_out[:, 1, :] = dA
    elseif nz_orig == 1
        dA_out[:, :, 1] = dA
    end
    println("Checksum: $(sum(dA_out))")
    return dA_out
end


function getUnitNormals(x, y, z, interiorPoint)
    # computes unit normal vector for a z-surface defined by x and y
    # to do side walls or end-caps, will need to switch order of arguments
    nx, ny, nz = size(x)
    ux = Array{Float64,3}(undef, nx, ny, nz)
    uy = Array{Float64,3}(undef, nx, ny, nz)
    uz = Array{Float64,3}(undef, nx, ny, nz)

    if nz == 1
        for i in 1:nx
            for j in 1:ny
                for k in 1:nz
                    # dz/dx
                    if i == 1
                        # calc dz/dx with forward difference
                        dzdx = (z[i+1, j, k] - z[i, j, k]) / (x[i+1, j, k] - x[i, j, k])
                    elseif i == nx
                        # calc dz/dx with backward difference
                        dzdx = (z[i, j, k] - z[i-1, j, k]) / (x[i, j, k] - x[i-1, j, k])
                    else
                        # calc dz/dx with central difference
                        dzdx = (z[i+1, j, k] - z[i-1, j, k]) / (x[i+1, j, k] - x[i-1, j, k])
                    end

                    # dz/dy
                    if j == 1
                        # calc dz/dy with forward difference
                        dzdy = (z[i, j+1, k] - z[i, j, k]) / (y[i, j+1, k] - y[i, j, k])
                    elseif j == ny
                        # calc dz/dy with backward difference
                        dzdy = (z[i, j, k] - z[i, j-1, k]) / (y[i, j, k] - y[i, j-1, k])
                    else
                        # calc dz/dy with central difference
                        dzdy = (z[i, j+1, k] - z[i, j-1, k]) / (y[i, j+1, k] - y[i, j-1, k])
                    end

                    # unit vector
                    denom = sqrt(dzdx^2 + dzdy^2 + 1)
                    ux[i, j, k] = -dzdx / denom
                    uy[i, j, k] = -dzdy / denom
                    uz[i, j, k] = 1 / denom
                end
            end
        end

    elseif ny == 1
        for i in 1:nx
            for j in 1:ny
                for k in 1:nz
                    # dy/dx
                    if i == 1
                        # calc dy/dx with forward difference
                        dydx = (y[i+1, j, k] - y[i, j, k]) / (x[i+1, j, k] - x[i, j, k])
                    elseif i == nx
                        # calc dy/dx with backward difference
                        dydx = (y[i, j, k] - y[i-1, j, k]) / (x[i, j, k] - x[i-1, j, k])
                    else
                        # calc dy/dx with central difference
                        dydx = (y[i+1, j, k] - y[i-1, j, k]) / (x[i+1, j, k] - x[i-1, j, k])
                    end

                    # dy/dz
                    if k == 1
                        # compute dy/dz with forward difference
                        dydz = (y[i, j, k+1] - y[i, j, k]) / (z[i, j, k+1] - z[i, j, k])
                    elseif k == nz
                        # compute dy/dz with backward difference
                        dydz = (y[i, j, k] - y[i, j, k-1]) / (z[i, j, k] - z[i, j, k-1])
                    else
                        # use central difference
                        dydz = (y[i, j, k+1] - y[i, j, k-1]) / (z[i, j, k+1] - z[i, j, k-1])
                    end

                    # unit vector
                    denom = sqrt(dydx^2 + dydz^2 + 1)
                    ux[i, j, k] = -dydx / denom
                    uy[i, j, k] = 1 / denom
                    uz[i, j, k] = -dydz / denom
                end
            end
        end


    elseif nx == 1
        for i in 1:nx
            for j in 1:ny
                for k in 1:nz
                    # dx/dy
                    if j == 1
                        # compute dx/dy with forward difference
                        dxdy = (x[i, j+1, k] - x[i, j, k]) / (y[i, j+1, k] - y[i, j, k])
                    elseif j == ny
                        # compute dx/dy with backward difference
                        dxdy = (x[i, j, k] - x[i, j-1, k]) / (y[i, j, k] - y[i, j-1, k])
                    else
                        # use central difference
                        dxdy = (x[i, j+1, k] - x[i, j-1, k]) / (y[i, j+1, k] - y[i, j-1, k])
                    end

                    # dx/dz
                    if k == 1
                        # compute dx/dz with forward difference
                        dxdz = (x[i, j, k+1] - x[i, j, k]) / (z[i, j, k+1] - z[i, j, k])
                    elseif k == nz
                        # compute dx/dz with backward difference
                        dxdz = (x[i, j, k] - x[i, k, k-1]) / (z[i, j, k] - z[i, j, k-1])
                    else
                        # use central difference
                        dxdz = (x[i, j, k+1] - x[i, j, k-1]) / (z[i, j, k+1] - z[i, j, k-1])
                    end

                    # unit vector
                    denom = sqrt(dxdy^2 + dxdz^2 + 1)
                    ux[i, j, k] = 1 / denom
                    uy[i, j, k] = -dxdy / denom
                    uz[i, j, k] = -dxdz / denom
                end
            end
        end


    end

    # flip unit normals if need be
    point = [interiorPoint[1], interiorPoint[2], interiorPoint[3]]
    #point = [0.0, 0.0, 0.0]
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                N = [x[i, j, k], y[i, j, k], z[i, j, k]]
                inVec = N .- point
                dotProd = ux[i, j, k] * inVec[1] + uy[i, j, k] * inVec[2] + uz[i, j, k] * inVec[3]
                if dotProd < 0.0
                    # flip unit normal
                    #println("we doin the flip")
                    ux[i, j, k], uy[i, j, k], uz[i, j, k] = -ux[i, j, k], -uy[i, j, k], -uz[i, j, k]
                end
            end
        end
    end

    return ux, uy, uz
end

function getUnitNormals2(x, y, z, interiorPoint)
    # computes unit normals using cross product

    nx, ny, nz = size(x)
    ux = Array{Float64,3}(undef, nx, ny, nz)
    uy = Array{Float64,3}(undef, nx, ny, nz)
    uz = Array{Float64,3}(undef, nx, ny, nz)

    if nx == 1
        for i in 1:nx
            for j in 1:ny
                for k in 1:nz
                    p1 = [x[i, j, k], y[i, j, k], z[i, j, k]]
                    if j == ny
                        # pull from behind
                        p2 = [x[i, j-1, k], y[i, j-1, k], z[i, j-1, k]]
                    else
                        # pull from front
                        p2 = [x[i, j+1, k], y[i, j+1, k], z[i, j+1, k]]
                    end
                    if k == nz
                        # pull from behind
                        p3 = [x[i, j, k-1], y[i, j, k-1], z[i, j, k-1]]
                    else
                        # pull from front
                        p3 = [x[i, j, k+1], y[i, j, k+1], z[i, j, k+1]]
                    end

                    v1 = p2 - p1
                    v2 = p3 - p1

                    crossOut = cross(v1, v2)
                    mag = norm(crossOut)
                    ux[i, j, k] = crossOut[1] / mag
                    uy[i, j, k] = crossOut[2] / mag
                    uz[i, j, k] = crossOut[3] / mag
                end
            end
        end

    elseif ny == 1
        for i in 1:nx
            for j in 1:ny
                for k in 1:nz
                    p1 = [x[i, j, k], y[i, j, k], z[i, j, k]]
                    if i == nx
                        # pull from behind
                        p2 = [x[i-1, j, k], y[i-1, j, k], z[i-1, j, k]]
                    else
                        # pull from front
                        p2 = [x[i+1, j, k], y[i+1, j, k], z[i+1, j, k]]
                    end
                    if k == nz
                        # pull from behind
                        p3 = [x[i, j, k-1], y[i, j, k-1], z[i, j, k-1]]
                    else
                        # pull from front
                        p3 = [x[i, j, k+1], y[i, j, k+1], z[i, j, k+1]]
                    end

                    v1 = p2 - p1
                    v2 = p3 - p1

                    crossOut = cross(v1, v2)
                    mag = norm(crossOut)
                    ux[i, j, k] = crossOut[1] / mag
                    uy[i, j, k] = crossOut[2] / mag
                    uz[i, j, k] = crossOut[3] / mag

                end
            end
        end

    elseif nz == 1
        for i in 1:nx
            for j in 1:ny
                for k in 1:nz
                    p1 = [x[i, j, k], y[i, j, k], z[i, j, k]]
                    if i == nx
                        # pull from behind
                        p2 = [x[i-1, j, k], y[i-1, j, k], z[i-1, j, k]]
                    else
                        # pull from front
                        p2 = [x[i+1, j, k], y[i+1, j, k], z[i+1, j, k]]
                    end
                    if j == ny
                        # pull from behind
                        p3 = [x[i, j-1, k], y[i, j-1, k], z[i, j-1, k]]
                    else
                        # pull from front
                        p3 = [x[i, j+1, k], y[i, j+1, k], z[i, j+1, k]]
                    end

                    v1 = p2 - p1
                    v2 = p3 - p1

                    crossOut = cross(v1, v2)
                    mag = norm(crossOut)
                    ux[i, j, k] = crossOut[1] / mag
                    uy[i, j, k] = crossOut[2] / mag
                    uz[i, j, k] = crossOut[3] / mag

                end
            end
        end

    end



    # flip unit normals if need be
    point = [interiorPoint[1], interiorPoint[2], interiorPoint[3]]
    #point = [0.0, 0.0, 0.0]
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                N = [x[i, j, k], y[i, j, k], z[i, j, k]]
                inVec = N .- point
                dotProd = ux[i, j, k] * inVec[1] + uy[i, j, k] * inVec[2] + uz[i, j, k] * inVec[3]
                if dotProd < 0.0
                    # flip unit normal
                    ux[i, j, k], uy[i, j, k], uz[i, j, k] = -ux[i, j, k], -uy[i, j, k], -uz[i, j, k]
                end
            end
        end
    end

    # get rid of any NaN's
    ux = replace(ux, NaN => 0.0)
    uy = replace(uy, NaN => 0.0)
    uz = replace(uz, NaN => 0.0)

    return ux, uy, uz
end




function getTriangleArea(p1, p2, p3)
    # computes area of a triangle from 3 points
    AB = p2 - p1
    AC = p3 - p1
    A = 0.5 * norm(cross(AB, AC))
    return A
end

function getMidpoint(p1, p2)
    # gets midpoint of 2 points
    x = (p1[1] + p2[1]) / 2
    y = (p1[2] + p2[2]) / 2
    z = (p1[3] + p2[3]) / 2
    p3 = Vector{Float64}([x, y, z])
    return p3
end


function computeCutoffFrequency(x, y, z, c, xFind)
    # computes cutoff frequency of grid
    nx, ny, nz = size(x)

    pFind = [xFind, 0, 0]                                                               # point we will look for
    d = sqrt.((x .- pFind[1]) .^ 2 + (y .- pFind[2]) .^ 2 + (z .- pFind[3]) .^ 2)       # distance from point we are looking for

    ind_min = argmin(d)

    if nx == 1
        # negect dim=1
        d_eta = sqrt.(diff(x, dims=2) .^ 2 + diff(y, dims=2) .^ 2 + diff(z, dims=2) .^ 2)
        d_zeta = sqrt.(diff(x, dims=3) .^ 2 + diff(y, dims=3) .^ 2 + diff(z, dims=3) .^ 2)

        ds = max(d_eta[ind_min[1], ind_min[2], ind_min[3]], d_zeta[ind_min[1], ind_min[2], ind_min[3]])
    elseif ny == 1
        # neglect dim=2
        d_xi = sqrt.(diff(x, dims=1) .^ 2 + diff(y, dims=1) .^ 2 + diff(z, dims=1) .^ 2)
        d_zeta = sqrt.(diff(x, dims=3) .^ 2 + diff(y, dims=3) .^ 2 + diff(z, dims=3) .^ 2)

        ds = max(d_xi[ind_min[1], ind_min[2], ind_min[3]], d_zeta[ind_min[1], ind_min[2], ind_min[3]])
    elseif nz == 1
        # neglect dim=3
        d_xi = sqrt.(diff(x, dims=1) .^ 2 + diff(y, dims=1) .^ 2 + diff(z, dims=1) .^ 2)
        d_eta = sqrt.(diff(x, dims=2) .^ 2 + diff(y, dims=2) .^ 2 + diff(z, dims=2) .^ 2)

        ds = max(d_xi[ind_min[1], ind_min[2], ind_min[3]], d_eta[ind_min[1], ind_min[2], ind_min[3]])
    end

    println(ind_min)
    f_cut = c / ds


    return f_cut
end