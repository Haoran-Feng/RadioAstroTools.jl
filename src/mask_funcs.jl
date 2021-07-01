# --------------------------------------------------
# mask functions 
                    
function apply_mask!(datacube::AbstractSpectralCube, maskcube::AbstractSpectralCube)
    datacube.data[isnan.(maskcube.data)] .= NaN;
    datacube.data[iszero.(maskcube.data)] .= NaN;
    return datacube
end

function apply_mask(datacube::AbstractSpectralCube, maskcube::AbstractSpectralCube)
    datacube2 = deepcopy(datacube)
    datacube2.data[isnan.(maskcube.data)] .= NaN;
    datacube2.data[iszero.(maskcube.data)] .= NaN;
    return datacube2
end

function apply_mask!(datacube::AbstractSpectralCube, maskmap::AbstractMap)
    mask2d = maskmap.data
    s = size(maskmap.data)
    mask3d = zeros(eltype(maskmap.data), (s[1], s[2], size(datacube.data, 3)))
    mask3d .= mask2d
    datacube.data[isnan.(mask3d)] .= NaN;
    datacube.data[iszero.(mask3d)] .= NaN;
    return datacube
end

function apply_mask!(datacube::AbstractSpectralCube, maskarr::Union{Array{T, 2} where T <: Real, BitArray})
    mask2d = maskarr
    s = size(mask2d)
    T = eltype(maskarr)
    mask3d = zeros(T, (s[1], s[2], size(datacube.data, 3)))
    mask3d .= mask2d
    datacube.data[isnan.(mask3d)] .= NaN;
    datacube.data[iszero.(mask3d)] .= NaN;
    return datacube
end