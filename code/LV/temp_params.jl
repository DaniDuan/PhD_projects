using Random, Distributions

function temp_func(T, Tr, B0, Ea, Ed, Th)
    k =  0.0000862 # Boltzman constant
    B = B0 .* exp.(-Ea/k .* (1/T.-1/Tr))# ./(1 .+ (Ea./(Ed - Ea)) .* exp.((Ed/k) .* (1 ./Th .- 1/T)))
    return B
end

function Params(N, T)
    Tr = 10+273.15
    r0 = rand(Normal(0.21, 0.01), N) # from data
    Ea_r = rand(Normal(0.95, 0.16), N) # from data
    Ed_r = fill(4, N)
    Th_r = rand(Normal(25+273.15, 5), N) # from data
    r = temp_func(T, Tr, r0, Ea_r, Ed_r, Th_r) 
    # assuming every α has a seperated TPC
    α0 =  rand(Normal(-1.18, 0.5), N,N) # mean value from data, test variation
    ran_Eα = rand(Normal(2.2, 0.5), N)
    Ea_α = reshape(repeat(ran_Eα, N),(N,N)) # mean value from data, test variation
    Ed_α = fill(4, N, N)
    Th_α = rand(Normal(24+273.15, 1), N, N) # test variation
    α = temp_func(T, Tr, α0, Ea_α, Ed_α, Th_α)
    return (N = N, α = α, r = r)
end

