
# index functions
Base.getindex(cube::AbstractSpectralCube, i::Union{Colon, Integer}, j::Union{Colon, Integer}, k::Union{Colon, Integer}) = cube.data[i, j, k]
function Base.getindex(cube::AbstractSpectralCube, i::AbstractFloat, j::AbstractFloat, k::Colon)
    ii = closest_index(cube.xaxis, i)
    jj = closest_index(cube.yaxis, j)
    return cube.data[ii, jj, :]
end
