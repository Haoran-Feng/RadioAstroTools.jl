function moment0(spectrum::Vector{<:Real}, channel_width::Real)
    return sum(spectrum) * channel_width;
end


function moment0(cube::AbstractSpectralCube)
    m0 = mapslices(i->moment0(i, cube.z_resolution), cube.data; dims=[3])[:, :, 1]
    header = header2d_of(cube, cube.header["LINE"], "Moment0")
    return create_map(m0, header)
end