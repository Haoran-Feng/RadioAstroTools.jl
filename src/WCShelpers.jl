
# --------------------------------------------------
# WCS functions 

function wcs2d_from_header(header::FITSHeader)
    cdelt = [header["CDELT" * string(i)] for i in [1, 2]]
    ctype = [header["CTYPE" * string(i)] for i in [1, 2]]
    crpix = convert(Array{Float64, 1}, [header["CRPIX" * string(i)] for i in [1, 2]])
    crval = [header["CRVAL" * string(i)] for i in [1, 2]];
    wcs = WCSTransform(2;
        cdelt = cdelt,
        ctype = ctype,
        crpix = crpix,
        crval = crval
    );
    return wcs
end

function wcs3d_from_header(header::FITSHeader)
    cdelt = [header["CDELT" * string(i)] for i in [1, 2, 3]]
    ctype = [header["CTYPE" * string(i)] for i in [1, 2, 3]]
    crpix = convert(Array{Float64, 1}, [header["CRPIX" * string(i)] for i in [1, 2, 3]])
    crval = [header["CRVAL" * string(i)] for i in [1, 2, 3]];
    wcs = WCSTransform(3;
        cdelt = cdelt,
        ctype = ctype,
        crpix = crpix,
        crval = crval
    );
    return wcs
end

function parse_fitsheader(header::String)
    n = length(header)
    lines = [header[(i - 1)*80 + 1: i * 80] for i in 1:(length(header) ÷ 80)]
    keys = []
    values = []
    comments = []
    for i in lines
        if (startswith(i, "COMMENT")) || length(strip(i)) == 0
            continue
        end
        k, v = split(i, "=")
        push!(keys, string(strip(k)))
        val, comment = split(v, "/")
        val = string(strip(val))
        val_types = [tryparse(T, val) for T in [Int, Float64]]
        if val_types[1] != nothing
            val = val_types[1]
        elseif val_types[2] != nothing
            val = val_types[2]
        else
            val = replace(val, "'" => "")
        end
        
        push!(values, val)
        push!(comments, string(strip(comment)))
    end
    return convert(Array{String, 1}, keys), values, convert(Array{String, 1}, comments)
end

function header2d_of(cube::AbstractSpectralCube, line::AbstractString, content_description::AbstractString)
    hkeys, hvalues, comments = wcs2d_from_header(cube.header) |> WCS.to_header |> parse_fitsheader
    hkeys = ["SIMPLE" ; hkeys; "DATA"; "LINE"]
    hvalues = [true; hvalues; content_description; line];
    comments = [cube.header.comments[1]; comments; ""; ""]
    header = FITSHeader(hkeys, hvalues, comments);
    return header
end
                    