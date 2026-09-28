# ---------------------------------------------------------------------------
# Material data sets and kinetic parameters
# ---------------------------------------------------------------------------

"""
Kinetic and geometric parameters of the model (all dimensionless).

`kappa` conductivity, `beta` Faraday factor, `A2`/`A1` positive/negative
exponential slopes, `B` radius of the insulated container `D = B_B(0)`.
"""
Base.@kwdef struct Params
    kappa::Float64 = 0.5
    beta::Float64 = 0.12
    A2::Float64 = 1.0
    A1::Float64 = 1.0
    B::Float64 = 2.0
end

withkappa(P::Params, kappa) = Params(kappa=kappa, beta=P.beta, A2=P.A2, A1=P.A1, B=P.B)

"""
A material data set: exchange-current scale `i0(x1,x2)`, its gradient, the
equilibrium (mixed) potential `phieq(x1,x2)`, and material bounds
over the closed container `B_B(0)`:
`imin = min i0`, `a1 = ||grad i0||_inf`, `a2 = ||D^2 i0||_inf`.
"""
struct Material{F1,F2,F3}
    name::String
    i0::F1
    gradi0::F2
    phieq::F3
    imin::Float64
    a1::Float64
    a2::Float64
end

Base.show(io::IO, m::Material) = print(io, "Material(", m.name, ")")

# M0: homogeneous data -> exact parallel-curve solutions (the parallel-curve lemma).
const M0 = Material("M0", (x, y) -> 1.0, (x, y) -> (0.0, 0.0), (x, y) -> 0.15, 1.0, 0.0, 0.0)

# M1: frozen single Fourier mode on the interface; phi_eq enters the discrete
# problem only through its trace, and ||phi_eq||_inf <= delta on the closed disk.
M1(n::Int, delta::Real) = Material("M1(n=$n,delta=$delta)",
    (x, y) -> 1.0, (x, y) -> (0.0, 0.0),
    (x, y) -> delta * cos(n * atan(y, x)), 1.0, 0.0, 0.0)

# M2 family with contrast c: i0 = 1 + c(0.25 x1 + 0.10 x2^2).
function M2family(c::Real)
    c = float(c)
    Material("M2(c=$c)",
        (x, y) -> 1.0 + c * (0.25 * x + 0.10 * y^2),
        (x, y) -> (0.25 * c, 0.20 * c * y),
        (x, y) -> 0.30 * (x^2 - y^2) + 0.15 * x * y,
        1.0 - 0.5 * c, c * hypot(0.25, 0.40), 0.20 * c)
end
const M2 = M2family(1.0)

"Dissolution-excess constant c_{A2,A1} of the dissolution-excess lemma."
excess_constant(P::Params) = P.A2 * P.A1 * min(P.A2, P.A1) / (2 * (P.A2 + P.A1))

"""
Dirichlet-to-Neumann eigenvalue lambda_n(R) of the annulus R < |x| < B with
insulated outer circle, the sharpness lemma of the paper; lambda_0 = 0.
"""
function dtn_eigenvalue(n::Integer, R::Real, B::Real)
    n == 0 && return 0.0
    n = abs(n)
    q = (R / B)^(2n)
    return n / R * (1 - q) / (1 + q)
end
