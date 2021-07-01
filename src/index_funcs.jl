
# index functions
Base.getindex(cube::AbstractSpectralCube, i::Integer, j::Integer, k::Integer) = cube.data[i, j, k]
function Base.getindex(cube::AbstractSpectralCube, i::AbstractFloat, j::AbstractFloat, k::Colon)
    ii = MapAndCube.closest_index(cube.xaxis, i)
    jj = MapAndCube.closest_index(cube.yaxis, j)
    return cube.data[ii, jj, :]
end