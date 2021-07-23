
mutable struct SpectralCube <: AbstractSpectralCube
    data::Array{T, 3} where T <: Real
    value_unit::Unitful.FreeUnits
    header::FITSHeader
    wcs::WCSTransform
    xaxis::Array{T, 1} where T <: Real
    yaxis::Array{T, 1} where T <: Real
    zaxis::Array{T, 1} where T <: Real
    xy_unit::Unitful.FreeUnits
    z_unit::Unitful.FreeUnits
    xy_resolution::Real
    z_resolution::Real
end


# --------------------------------------------------
# Cube functions 
function spectralcube(data::Array{T, 3} where T <: Real, header::FITSHeader)
    wcs = wcs3d_from_header(header)

    idx = zeros((3, size(data, 3)))
    idx[3, :] = 1:size(data, 3)
    z_axis = pix_to_world(wcs, idx)[3, :];
    if ("CUNIT3" ∈ keys(header) && occursin("km", header["CUNIT3"]))
        z_unit = u"km / s";
    end
    z_unit = u"m / s"
#     z_axis = to_kps_val.(z_axis * z_unit);
#     z_unit = u"km / s";

    idx = zeros((3, size(data, 2)))
    idx[2, :] = 1:size(data, 2)
    y_axis = pix_to_world(wcs, idx)[2, :];

    idx = zeros((3, size(data, 1)))
    idx[1, :] = 1:size(data, 1)
    x_axis = pix_to_world(wcs, idx)[1, :];

    xy_unit = u"°";
    
    val_unit = u"K";
    
    xy_resolution = abs(header["CDELT1"])
    z_resolution = abs(z_axis[1] - z_axis[2])
    
    return SpectralCube(data, val_unit, header, wcs, x_axis, y_axis, z_axis, xy_unit, z_unit, xy_resolution, z_resolution);
end

            
function subcube(cube::SpectralCube, xyz_range::CubeRange)
    xmin_incube = closest_index(cube.xaxis, xyz_range.xmin)
    xmax_incube = closest_index(cube.xaxis, xyz_range.xmax)
    ymin_incube = closest_index(cube.yaxis, xyz_range.ymin)
    ymax_incube = closest_index(cube.yaxis, xyz_range.ymax)
    zmin_incube = closest_index(cube.zaxis, xyz_range.zmin)
    zmax_incube = closest_index(cube.zaxis, xyz_range.zmax)

    # make sure the indeces are ascending
    xmin_incube, xmax_incube = sort([xmin_incube, xmax_incube])
    ymin_incube, ymax_incube = sort([ymin_incube, ymax_incube])
    zmin_incube, zmax_incube = sort([zmin_incube, zmax_incube])
                
    if (cube.data |> size |> length) == 4
        data = cube.data[xmin_incube:xmax_incube, ymin_incube:ymax_incube, zmin_incube: zmax_incube, 1]
    else
        data = cube.data[xmin_incube:xmax_incube, ymin_incube:ymax_incube, zmin_incube: zmax_incube]
    end 
    
                
    header = deepcopy(cube.header)
    header["CRPIX1"] = header["CRPIX1"] - xmin_incube + 1
    header["CRPIX2"] = header["CRPIX2"] - ymin_incube + 1
    header["CRPIX3"] = header["CRPIX3"] - zmin_incube + 1
    header["NAXIS"] = 3
    header["NAXIS1"] = size(data, 1)
    header["NAXIS2"] = size(data, 2)
    header["NAXIS3"] = size(data, 3)
                
    return spectralcube(data, header)
end
                    
                    
function map_to_new_zaxis!(cube1::AbstractSpectralCube, cube2::AbstractSpectralCube)
    # this method will modify cube2
    # from cube1 to cube2
    new_data = zeros(eltype(cube1.data), size(cube2.data))
    for i = 1:size(cube1.data, 3)
        v = cube1.zaxis[i]
        new_ch_idx = closest_index(cube2.zaxis, v)
        new_data[:, :, new_ch_idx] = cube1.data[:, :, i]
    end
    new_data[new_data .== 0] .= NaN;
    cube2.data = new_data;
    return cube2
end

function map_to_new_zaxis!(cube::AbstractSpectralCube, new_zaxis::AbstractArray) 
    # this method will modify cube
    new_data = zeros(eltype(cube.data), tuple([size(cube.data)[1:2]...,  length(new_zaxis)]...) )
    for (i, v) in enumerate(new_zaxis)
        old_i = closest_index(cube.zaxis, v)
        new_data[:, :, i] = cube.data[:, :, old_i]
    end
    # new_data[new_data .== 0] .= NaN;
    cube.data = new_data;
    cube.zaxis = new_zaxis
    return cube
end

function voxel_num(maskcube::AbstractSpectralCube)
    maskcube.data .|> isnan .|> (!) |> sum
end

function pixel_num(maskcube::AbstractSpectralCube)
    maskcube.data .|> isnan .|> (!) |> (x->sum(x, dims=3) .!= 0) |> sum
end

function avg_channel_num(maskcube::AbstractSpectralCube)
    voxel_num(maskcube) / pixel_num(maskcube)
end
                    
function valid_map(masked_datacube::SpectralCube, rmsmap::Map; n_sigma=3)
    cover_map = masked_datacube.data .|> isnan .|> (!) |> x->sum(x, dims=3)[:, :, 1];
    data = masked_datacube.data[:, :, :] # make a copy
    data[isnan.(data)] .= 0.;
    sum_map = sum(data, dims=3)[:, :, 1];
    vm = sum_map .> (n_sigma .* rmsmap.data .* sqrt.(cover_map))
    return vm 
end


function three_channels_13(spec::Vector, rms::Real)
    Nσ = 2.0
    max_idx = argmax(spec)
    if max_idx >= 2 && max_idx <= (length(spec) - 1)
        chs = [max_idx - 1, max_idx, max_idx + 1]
    else
        for i=1:length(spec) - 2
            chs = [i, i + 1, i + 2]
            if all(spec[chs] .>= Nσ * rms)
                return true
            end
        end
        return false
    end
    return all(spec[chs] .>= Nσ * rms)
end

function three_channels_18(spec::Vector, rms::Real)
    Nσ = 1.5
    max_idx = argmax(spec)
    if max_idx >= 2 && max_idx <= (length(spec) - 1)
        chs = [max_idx - 1, max_idx, max_idx + 1]
    else
        for i=1:length(spec) - 2
            chs = [i, i + 1, i + 2]
            if all(spec[chs] .>= Nσ * rms)
                return true
            end
        end
        return false
    end
    return all(spec[chs] .>= Nσ * rms)
end


function Base.string(cube::SpectralCube)
    # TODO
    return "SpectralCube info"
end