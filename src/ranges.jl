
# --------------------------------------------------
# Ranges

mutable struct MapRange
    xmin::Real
    xmax::Real
    ymin::Real
    ymax::Real
end

mutable struct CubeRange
    xmin::Real
    xmax::Real
    ymin::Real
    ymax::Real
    zmin::Real
    zmax::Real
end

function xy_range(map_::AbstractMap)
    MapRange(map_.xaxis[begin], map_.xaxis[end], map_.yaxis[begin], map_.yaxis[end])
end
        
function xy_range(cr::CubeRange)
    MapRange(cr.xmin, cr.xmax, cr.ymin, cr.ymax)
end

function xy_range(cube::AbstractSpectralCube)
    cr = xyz_range(cube)
    return MapRange(cr.xmin, cr.xmax, cr.ymin, cr.ymax)
end

function xyz_range(cube::AbstractSpectralCube)
    CubeRange(cube.xaxis[begin], cube.xaxis[end], cube.yaxis[begin], cube.yaxis[end], cube.zaxis[begin], cube.zaxis[end])
end

