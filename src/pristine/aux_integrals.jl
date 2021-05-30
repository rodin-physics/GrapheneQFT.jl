# Nearest-neighbor hopping parameter
const NN_hopping = 2.8
# Integrals used in computing the propagator
@inline function Ω_Integrand(z, u, v, x::Float64)
    t = NN_hopping
    W = ((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x)
    return (
        exp(1.0im * (u - v) * x) / cos(x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v)) / (√(W - 1) * √(W + 1))
    )
end

@inline function Ω(z, u, v)
    t = NN_hopping
    return ((quadgk(
        x -> Ω_Integrand(z, u, v, x) / (8.0 * π * t^2),
        0.0,
        2.0 * π,
    ))[1])
end

@inline function Ωp_Integrand(z, u, v, x::Float64)
    t = NN_hopping
    W = ((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x)
    return (
        2 *
        exp(1.0im * (u - v) * x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v + 1)) / (√(W - 1) * √(W + 1))
    )
end

@inline function Ωp(z, u, v)
    t = NN_hopping
    return ((quadgk(
        x -> Ωp_Integrand(z, u, v, x) / (8.0 * π * t^2),
        0.0,
        2.0 * π,
    ))[1])
end

@inline function Ωn_Integrand(z, u, v, x::Float64)
    t = NN_hopping
    W = ((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x)
    return (
        2 *
        exp(1.0im * (u - v) * x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v - 1)) / (√(W - 1) * √(W + 1))
    )
end

@inline function Ωn(z, u, v)
    t = NN_hopping
    return ((quadgk(
        x -> Ωn_Integrand(z, u, v, x) / (8.0 * π * t^2),
        0.0,
        2.0 * π,
    ))[1])
end
