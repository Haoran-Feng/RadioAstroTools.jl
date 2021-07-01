# --------------------------------------------------
# IO functions 

function read_map_from_file(fitsfilename::String; hdu_id::Int=1)
    hdu = FITS(fitsfilename)[hdu_id]
    data = read(hdu);
    header = read_header(hdu);
    wcs = wcs2d_from_header(header)

    idx = zeros((2, size(data, 2)))
    idx[2, :] = 1:size(data, 2)
    y_axis = pix_to_world(wcs, idx)[2, :];

    idx = zeros((2, size(data, 1)))
    idx[1, :] = 1:size(data, 1)
    x_axis = pix_to_world(wcs, idx)[1, :];

    xy_unit = u"°";
    
    val_unit = u"K";
    
    xy_resolution = abs(x_axis[1] - x_axis[2])
    
    return Map(data, val_unit, header, wcs, x_axis, y_axis, xy_unit, xy_resolution);
end

function read_cube_from_file(fitsfilename::String; hdu_id::Int=1)
    hdu = FITS(fitsfilename)[hdu_id]
    if ndims(hdu) == 4
        data = read(hdu, :, :, :, 1);
    else
        data = read(hdu)
    end
    header = read_header(hdu);
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
            

function save(filename::AbstractString, map_::AbstractMap)
    f = FITS(filename, "w")
    write(f, map_.data, header=map_.header)
    close(f)
    return nothing
end

function save(filename::String, cube::AbstractSpectralCube)
    FITS(filename, "w") do f
        write(f, cube.data, header=cube.header)
    end
end
    
        
function header2d_of(cube::AbstractSpectralCube, line::AbstractString, content_description::AbstractString)
    hkeys, hvalues, comments = MapAndCube.wcs2d_from_header(cube.header) |> WCS.to_header |> MapAndCube.parse_fitsheader
    hkeys = ["SIMPLE" ; hkeys; "DATA"; "LINE"]
    hvalues = [true; hvalues; content_description; line];
    comments = [cube.header.comments[1]; comments; ""; ""]
    header = FITSHeader(hkeys, hvalues, comments);
    return header
end