using Magnetostatics: MagneticField
using StaticArrays

"""
A cusp field configuration where ``B_x = B_0 * (z/L) / (1 + a * (z/L)^2)^2``
"""
@kwdef struct BConfig{𝐵, 𝐿} <: MagneticField{1}
    B0::𝐵
    By::𝐵 = 0.0
    Bz::𝐵 = 1.0
    L::𝐿 = 1.0
    a::Float64 = 0.0
end

function (c::BConfig)(z::Number)
    z̃ = z / c.L
    Bx = c.B0 * z̃ / (1 + c.a * z̃^2)^2
    return SA[Bx, c.By, c.Bz]
end

"A cusp field configuration where ``B_x = B_0 * tanh(z/L) / (1 + a * (z/L)^2)^2``"
@kwdef struct BConfig2{𝐵, 𝐿} <: MagneticField{1}
    B0::𝐵
    By::𝐵 = 0.0
    Bz::𝐵 = 1.0
    L::𝐿 = 1.0
    a::Float64 = 0.0
end

function (c::BConfig2)(z::Number)
    z̃ = z / c.L
    Bx = c.B0 * tanh(z̃) / (1 + c.a * z̃^2)^2
    return SA[Bx, c.By, c.Bz]
end
