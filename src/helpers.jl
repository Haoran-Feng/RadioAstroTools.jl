# --------------------------------------------------
# Helper functions 

function closest_index(arr::AbstractArray, val)
    arr .- val .|> abs |> argmin
end

function cart_index_to_array(idx::Array{CartesianIndex{3}, 1})
    a = zeros(Float32, (3, length(idx)));
    Threads.@threads for i=1:length(idx)
        a[:, i] .= idx[i].I
    end
    return a
end

function cart_index_to_array(idx::Array{CartesianIndex{2}, 1})
    a = zeros(Float32, (2, length(idx)));
    Threads.@threads for i=1:length(idx)
        a[:, i] .= idx[i].I
    end
    return a
end

function closest_index(arr::AbstractArray, val)
    arr .- val .|> abs |> argmin
end

                        
function to_kps_val(x)
    uconvert(u"km / s", x).val
end
            


