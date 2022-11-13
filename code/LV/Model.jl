include("./temp_params.jl")

# the ODE function
function GLV_model!(dx, x, p, t)
    dx .= x.*(p.r+p.α*x)
end
