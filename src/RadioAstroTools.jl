module RadioAstroTools

using FITSIO
using WCS
using Unitful, UnitfulAstro
using StatsBase

export AbstractSpectralCube, AbstractMap
export Map, SpectralCube
export create_map, spectralcube
export MapRange, CubeRange
export xy_range, xyz_range
export read_map_from_file, read_cube_from_file
export save
export submap, subcube
export apply_mask, apply_mask!
export copy


# abstract types
abstract type AbstractMap end
abstract type AbstractSpectralCube end

include("helpers.jl")
include("ranges.jl")
include("io.jl")
include("WCShelpers.jl")
include("maps.jl")
include("cubes.jl")
include("index_funcs.jl")
include("mask_funcs.jl")

end
