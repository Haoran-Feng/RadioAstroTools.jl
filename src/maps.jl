
mutable struct Map <: AbstractMap
    data::Array{T, 2} where T <: Real
    value_unit::Unitful.FreeUnits
    header::FITSHeader
    wcs::WCSTransform
    xaxis::Array{T, 1} where T <: Real
    yaxis::Array{T, 1} where T <: Real
    xy_unit::Unitful.FreeUnits
    xy_resolution::Real
end

# --------------------------------------------------
# map functions 
        
function create_map(data::Array{T, 2} where T <: Real, header::FITSHeader)
    wcs = wcs2d_from_header(header)

    idx = zeros((2, size(data, 2)))
    idx[2, :] = 1:size(data, 2)
    y_axis = pix_to_world(wcs, idx)[2, :];

    idx = zeros((2, size(data, 1)))
    idx[1, :] = 1:size(data, 1)
    x_axis = pix_to_world(wcs, idx)[1, :];

    xy_unit = u"°";
    
    val_unit = u"K";
    
    xy_resolution = abs(header["CDELT1"])
    
    return Map(data, val_unit, header, wcs, x_axis, y_axis, xy_unit, xy_resolution);
end

function submap(map_::Map, xy_range::MapRange)
    xmin_incube = closest_index(map_.xaxis, xy_range.xmin)
    xmax_incube = closest_index(map_.xaxis, xy_range.xmax)
    ymin_incube = closest_index(map_.yaxis, xy_range.ymin)
    ymax_incube = closest_index(map_.yaxis, xy_range.ymax)

    # make sure the indeces are ascending
    xmin_incube, xmax_incube = sort([xmin_incube, xmax_incube])
    ymin_incube, ymax_incube = sort([ymin_incube, ymax_incube])
                
    data = map_.data[xmin_incube:xmax_incube, ymin_incube:ymax_incube]
                
    header = deepcopy(map_.header)
    header["CRPIX1"] = header["CRPIX1"] - xmin_incube + 1
    header["CRPIX2"] = header["CRPIX2"] - ymin_incube + 1
    header["NAXIS"] = 2
    header["NAXIS1"] = size(data, 1)
    header["NAXIS2"] = size(data, 2)
                
    return create_map(data, header)
end

function Base.string(map_::Map)
    # TODO
    return "Map info"
end