# ---------------------------------------------------------------------------
# P1 finite elements on the fitted annular mesh
# ---------------------------------------------------------------------------

const GXI = ((1 - sqrt(3 / 5)) / 2, 0.5, (1 + sqrt(3 / 5)) / 2)   # 3-point Gauss on [0,1]
const GW = (5 / 18, 4 / 9, 5 / 18)

"""
Fitted annular mesh: `nt` equally spaced angular nodes, `nr` radial intervals
between the interface r = R(theta) (level 0) and the wall r = B (level nr).
Node (level k, angle j) has index k*nt + j. Each quadrilateral is split into
the two triangles (a,c,b), (a,d,c). The sparsity pattern of the stiffness
matrix and the positions of every local entry are precomputed once; the mesh
is rebuilt from new radii without changing connectivity.
"""
struct AnnulusMesh
    nt::Int
    nr::Int
    B::Float64
    grading::Float64
    theta::Vector{Float64}
    cs::Vector{Float64}
    sn::Vector{Float64}
    eta::Vector{Float64}
    tri::Matrix{Int}
    pattern::SparseMatrixCSC{Float64,Int}
    pos::Matrix{Int}
    epos::Matrix{Int}
    boundary::Vector{Int}
    interior::Vector{Int}
end

nnodes(m::AnnulusMesh) = (m.nr + 1) * m.nt

function _pattern_position(S::SparseMatrixCSC, i::Int, j::Int)
    r = S.colptr[j]:S.colptr[j+1]-1
    k = searchsortedfirst(view(S.rowval, r), i)
    @assert k <= length(r) && S.rowval[r[k]] == i "entry ($i,$j) not in pattern"
    return r[k]
end

function AnnulusMesh(nt::Int, nr::Int; B::Real=2.0, grading::Real=1.0)
    theta = theta_grid(nt)
    eta = (collect(0:nr) ./ nr) .^ grading
    ntri = 2 * nr * nt
    tri = Matrix{Int}(undef, 3, ntri)
    e = 0
    for k in 0:nr-1, j in 1:nt
        jp = j == nt ? 1 : j + 1
        a = k * nt + j
        b = k * nt + jp
        c = (k + 1) * nt + jp
        d = (k + 1) * nt + j
        e += 1; tri[:, e] .= (a, c, b)
        e += 1; tri[:, e] .= (a, d, c)
    end
    n = (nr + 1) * nt
    I = Vector{Int}(undef, 9ntri); J = similar(I)
    idx = 0
    for e in 1:ntri, β in 1:3, α in 1:3
        idx += 1
        I[idx] = tri[α, e]; J[idx] = tri[β, e]
    end
    pattern = sparse(I, J, ones(length(I)), n, n)
    fill!(pattern.nzval, 0.0)
    pos = Matrix{Int}(undef, 9, ntri)
    idx = 0
    for e in 1:ntri, β in 1:3, α in 1:3
        idx += 1
        pos[α+3(β-1), e] = _pattern_position(pattern, I[idx], J[idx])
    end
    epos = Matrix{Int}(undef, 4, nt)
    for j in 1:nt
        jp = j == nt ? 1 : j + 1
        epos[1, j] = _pattern_position(pattern, j, j)
        epos[2, j] = _pattern_position(pattern, jp, j)
        epos[3, j] = _pattern_position(pattern, j, jp)
        epos[4, j] = _pattern_position(pattern, jp, jp)
    end
    return AnnulusMesh(nt, nr, float(B), float(grading), theta, cos.(theta), sin.(theta), eta,
        tri, pattern, pos, epos, collect(1:nt), collect(nt+1:n))
end

"Node coordinates of the fitted mesh for the radii R (in place)."
function node_coordinates!(x::Vector{Float64}, y::Vector{Float64}, m::AnnulusMesh, R::AbstractVector)
    nt = m.nt
    @inbounds for k in 0:m.nr, j in 1:nt
        r = R[j] + m.eta[k+1] * (m.B - R[j])
        idx = k * nt + j
        x[idx] = r * m.cs[j]
        y[idx] = r * m.sn[j]
    end
    return nothing
end

function node_coordinates(m::AnnulusMesh, R::AbstractVector)
    n = nnodes(m)
    x = Vector{Float64}(undef, n); y = similar(x)
    node_coordinates!(x, y, m, R)
    return x, y
end

"Unit-conductivity stiffness values into `vals` (pattern order); returns the minimal triangle determinant."
function assemble_stiffness!(vals::Vector{Float64}, m::AnnulusMesh, x::Vector{Float64}, y::Vector{Float64})
    fill!(vals, 0.0)
    tri, pos = m.tri, m.pos
    mindet = Inf
    @inbounds for e in 1:size(tri, 2)
        a, b, c = tri[1, e], tri[2, e], tri[3, e]
        x1, y1, x2, y2, x3, y3 = x[a], y[a], x[b], y[b], x[c], y[c]
        det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        mindet = min(mindet, det)
        b1, b2, b3 = y2 - y3, y3 - y1, y1 - y2
        c1, c2, c3 = x3 - x2, x1 - x3, x2 - x1
        f = 1 / (2 * abs(det))
        bb = (b1, b2, b3); cc = (c1, c2, c3)
        for β in 1:3, α in 1:3
            vals[pos[α+3(β-1), e]] += f * (bb[α] * bb[β] + cc[α] * cc[β])
        end
    end
    return mindet
end

"P1 mass matrix (for L^2(Omega) norms of bulk fields)."
function mass_matrix(m::AnnulusMesh, x::Vector{Float64}, y::Vector{Float64})
    tri = m.tri
    ntri = size(tri, 2)
    I = Vector{Int}(undef, 9ntri); J = similar(I); V = Vector{Float64}(undef, 9ntri)
    idx = 0
    @inbounds for e in 1:ntri
        a, b, c = tri[1, e], tri[2, e], tri[3, e]
        det = abs((x[b] - x[a]) * (y[c] - y[a]) - (x[c] - x[a]) * (y[b] - y[a]))
        nodes = (a, b, c)
        for β in 1:3, α in 1:3
            idx += 1
            I[idx] = nodes[α]; J[idx] = nodes[β]
            V[idx] = det / 24 * (α == β ? 2.0 : 1.0)
        end
    end
    return sparse(I, J, V, nnodes(m), nnodes(m))
end

"Lengths of the reactive edges (j, j+1) of the polygonal interface."
function edge_lengths!(len::Vector{Float64}, m::AnnulusMesh, x, y)
    nt = m.nt
    @inbounds for j in 1:nt
        jp = j == nt ? 1 : j + 1
        len[j] = hypot(x[jp] - x[j], y[jp] - y[j])
    end
    return len
end

"""
Boundary mass matrix M_Gamma of the polygonal interface (consistent P1) and the
reaction-weighted mass matrix M_a with weight a = (A2+A1) i0 (nodal interpolant
of i0 at the Gauss points), both cyclic tridiagonal.
"""
function boundary_mass(len::Vector{Float64}, weight::Vector{Float64})
    nt = length(len)
    I = Int[]; J = Int[]; V = Float64[]
    for j in 1:nt
        jp = j == nt ? 1 : j + 1
        maa = mab = mbb = 0.0
        for q in 1:3
            ξ = GXI[q]; w = GW[q]
            wq = (1 - ξ) * weight[j] + ξ * weight[jp]
            maa += len[j] * w * wq * (1 - ξ)^2
            mab += len[j] * w * wq * (1 - ξ) * ξ
            mbb += len[j] * w * wq * ξ^2
        end
        push!(I, j); push!(J, j); push!(V, maa)
        push!(I, j); push!(J, jp); push!(V, mab)
        push!(I, jp); push!(J, j); push!(V, mab)
        push!(I, jp); push!(J, jp); push!(V, mbb)
    end
    return sparse(I, J, V, nt, nt)
end

# ---------------------------------------------------------------------------
# Nonlinear Robin problem: kappa K phi + Q_Gamma[i(x, phi) psi_i] = 0
# ---------------------------------------------------------------------------

"""
Solver state for the electrical problem on a mesh with fixed connectivity.
The Cholesky factorization keeps its symbolic analysis between calls.
"""
mutable struct RobinSolver
    mesh::AnnulusMesh
    P::Params
    x::Vector{Float64}
    y::Vector{Float64}
    Kvals::Vector{Float64}
    A::SparseMatrixCSC{Float64,Int}
    factor::Any
    phi::Vector{Float64}
    res::Vector{Float64}
    len::Vector{Float64}
    eqb::Vector{Float64}
    i0b::Vector{Float64}
    tol::Float64
    maxit::Int
end

function RobinSolver(mesh::AnnulusMesh, P::Params; tol=2e-12, maxit=40)
    n = nnodes(mesh)
    A = SparseMatrixCSC(n, n, copy(mesh.pattern.colptr), copy(mesh.pattern.rowval), zeros(nnz(mesh.pattern)))
    RobinSolver(mesh, P, zeros(n), zeros(n), zeros(nnz(mesh.pattern)), A, nothing,
        zeros(n), zeros(n), zeros(mesh.nt), zeros(mesh.nt), zeros(mesh.nt), tol, maxit)
end

"Set the geometry: coordinates, unit stiffness, edge lengths and nodal boundary data."
function set_geometry!(S::RobinSolver, mat::Material, R::AbstractVector)
    m = S.mesh
    node_coordinates!(S.x, S.y, m, R)
    mindet = assemble_stiffness!(S.Kvals, m, S.x, S.y)
    mindet > 1e-14 || error("degenerate fitted element (min det = $mindet)")
    edge_lengths!(S.len, m, S.x, S.y)
    @inbounds for j in 1:m.nt
        S.eqb[j] = mat.phieq(S.x[j], S.y[j])
        S.i0b[j] = mat.i0(S.x[j], S.y[j])
    end
    return mindet
end

"""
Residual r = kappa K phi + reaction(phi) into `r`; if `jac`, also the Jacobian
into S.A. Returns `false` if the exponential overflow guard triggers.
"""
function residual!(r::Vector{Float64}, S::RobinSolver, phi::Vector{Float64}, jac::Bool)
    m = S.mesh; P = S.P
    nt = m.nt
    κ, A2, A1 = P.kappa, P.A2, P.A1
    # bulk part
    Av = S.A.nzval
    if jac
        @. Av = κ * S.Kvals
    end
    Kfull = SparseMatrixCSC(size(S.A, 1), size(S.A, 2), S.A.colptr, S.A.rowval, S.Kvals)
    mul!(r, Kfull, phi)
    r .*= κ
    # reactive edges
    @inbounds for j in 1:nt
        jp = j == nt ? 1 : j + 1
        L = S.len[j]
        ra = rb = 0.0
        jaa = jab = jbb = 0.0
        for q in 1:3
            ξ = GXI[q]; w = GW[q] * L
            φq = (1 - ξ) * phi[j] + ξ * phi[jp]
            eq = (1 - ξ) * S.eqb[j] + ξ * S.eqb[jp]
            i0 = (1 - ξ) * S.i0b[j] + ξ * S.i0b[jp]
            δ = φq - eq
            abs(δ) > 60 && return false
            d = i0 * exp(A2 * δ)
            rr = i0 * exp(-A1 * δ)
            ia = d - rr
            ra += w * ia * (1 - ξ)
            rb += w * ia * ξ
            if jac
                ip = A2 * d + A1 * rr
                jaa += w * ip * (1 - ξ)^2
                jab += w * ip * (1 - ξ) * ξ
                jbb += w * ip * ξ^2
            end
        end
        r[j] += ra; r[jp] += rb
        if jac
            Av[m.epos[1, j]] += jaa
            Av[m.epos[2, j]] += jab
            Av[m.epos[3, j]] += jab
            Av[m.epos[4, j]] += jbb
        end
    end
    return true
end

"""
Newton's method with residual backtracking for the electrical problem at the
current geometry. `phi0` is the initial guess (default: nodal phi_eq extended
by its value at the same angle). Returns (iterations, final residual norm).
"""
function newton!(S::RobinSolver; phi0::Union{Nothing,Vector{Float64}}=nothing)
    m = S.mesh; n = nnodes(m); nt = m.nt
    phi = S.phi
    if phi0 === nothing
        @inbounds for k in 0:m.nr, j in 1:nt
            phi[k*nt+j] = S.eqb[j]
        end
    else
        copyto!(phi, phi0)
    end
    r = S.res
    cand = similar(phi); rc = similar(phi)
    iters = 0
    norm_r = Inf
    for it in 0:S.maxit
        ok = residual!(r, S, phi, true)
        ok || error("overflow in the reaction term at the Newton iterate")
        norm_r = maximum(abs, r)
        if norm_r < S.tol
            iters = it
            break
        end
        it == S.maxit && error("Newton failed to converge: residual $norm_r")
        Asym = Symmetric(S.A, :U)
        if S.factor === nothing
            S.factor = cholesky(Asym)
        else
            cholesky!(S.factor, Asym)
        end
        δ = S.factor \ (-r)
        damping = 1.0
        accepted = false
        while damping >= 2.0^-16
            @. cand = phi + damping * δ
            if residual!(rc, S, cand, false) && maximum(abs, rc) <= (1 - 1e-4 * damping) * norm_r
                copyto!(phi, cand)
                accepted = true
                break
            end
            damping *= 0.5
        end
        accepted || error("Newton backtracking failed at residual $norm_r")
        iters = it + 1
    end
    return iters, norm_r
end

"""
Interface diagnostics of the current solution: edge-quadrature integrals of the
net current, the dissolution current, i0, the the dissolution-excess lemma lower bound and the
nodal trace defect.
"""
function interface_diagnostics(S::RobinSolver)
    m = S.mesh; P = S.P; nt = m.nt
    A2, A1 = P.A2, P.A1
    phi = S.phi
    charge = diss = i0int = bound = 0.0
    cA = excess_constant(P)
    @inbounds for j in 1:nt
        jp = j == nt ? 1 : j + 1
        L = S.len[j]
        for q in 1:3
            ξ = GXI[q]; w = GW[q] * L
            φq = (1 - ξ) * phi[j] + ξ * phi[jp]
            eq = (1 - ξ) * S.eqb[j] + ξ * S.eqb[jp]
            i0 = (1 - ξ) * S.i0b[j] + ξ * S.i0b[jp]
            δ = φq - eq
            d = i0 * exp(A2 * δ); rr = i0 * exp(-A1 * δ)
            charge += w * (d - rr)
            diss += w * d
            i0int += w * i0
            bound += w * cA * i0 * δ^2
        end
    end
    return (charge=charge, dissolution=diss, i0_integral=i0int, excess=diss - i0int, excess_bound=bound)
end

"Nodal normal speed of the radial graph: -beta d(x_j, phi_j) sqrt(1 + (R_theta/R)^2)."
function nodal_speed(S::RobinSolver, R::AbstractVector)
    m = S.mesh; P = S.P
    N = metric_factor(R)
    nt = m.nt
    sp = Vector{Float64}(undef, nt)
    @inbounds for j in 1:nt
        sp[j] = -P.beta * S.i0b[j] * exp(P.A2 * (S.phi[j] - S.eqb[j])) * N[j]
    end
    return sp
end

"Discrete energy J_{kappa,h}(phi) = 1/2 phi'K phi + kappa^{-1} Q_Gamma[I(x,phi)]."
function discrete_energy(S::RobinSolver)
    m = S.mesh; P = S.P; nt = m.nt
    A2, A1 = P.A2, P.A1
    phi = S.phi
    Kfull = SparseMatrixCSC(size(S.A, 1), size(S.A, 2), S.A.colptr, S.A.rowval, S.Kvals)
    bulk = 0.5 * dot(phi, Kfull * phi)
    surf = 0.0
    @inbounds for j in 1:nt
        jp = j == nt ? 1 : j + 1
        L = S.len[j]
        for q in 1:3
            ξ = GXI[q]; w = GW[q] * L
            φq = (1 - ξ) * phi[j] + ξ * phi[jp]
            eq = (1 - ξ) * S.eqb[j] + ξ * S.eqb[jp]
            i0 = (1 - ξ) * S.i0b[j] + ξ * S.i0b[jp]
            δ = φq - eq
            surf += w * i0 * ((exp(A2 * δ) - 1) / A2 + (exp(-A1 * δ) - 1) / A1)
        end
    end
    return bulk + surf / P.kappa
end

"Solve the electrical problem at radii R; returns a NamedTuple with speed, diagnostics and Newton data."
function solve_electrical!(S::RobinSolver, mat::Material, R::AbstractVector; phi0=nothing)
    mindet = set_geometry!(S, mat, R)
    iters, resid = newton!(S; phi0=phi0)
    diag = interface_diagnostics(S)
    speed = nodal_speed(S, R)
    return (speed=speed, iterations=iters, residual=resid, min_triangle_area=mindet / 2, diag...)
end

# ---------------------------------------------------------------------------
# Linear Dirichlet problem and variational flux recovery
# ---------------------------------------------------------------------------

"""
Discrete harmonic extension with Dirichlet data on the interface and natural
(zero-flux) condition on the wall, with variational recovery of the flux
F_i = (K u_h)_i on interface nodes. The interior Cholesky factorization keeps
its symbolic analysis between calls (the interior pattern is fixed).
"""
mutable struct DirichletSolver
    mesh::AnnulusMesh
    x::Vector{Float64}
    y::Vector{Float64}
    Kvals::Vector{Float64}
    K::SparseMatrixCSC{Float64,Int}
    factor::Any
    len::Vector{Float64}
end

function DirichletSolver(mesh::AnnulusMesh)
    n = nnodes(mesh)
    K = SparseMatrixCSC(n, n, copy(mesh.pattern.colptr), copy(mesh.pattern.rowval), zeros(nnz(mesh.pattern)))
    DirichletSolver(mesh, zeros(n), zeros(n), K.nzval, K, nothing, zeros(mesh.nt))
end

"Assemble the unit stiffness at radii R."
function set_geometry!(D::DirichletSolver, R::AbstractVector)
    node_coordinates!(D.x, D.y, D.mesh, R)
    mindet = assemble_stiffness!(D.Kvals, D.mesh, D.x, D.y)
    mindet > 1e-14 || error("degenerate fitted element (min det = $mindet)")
    edge_lengths!(D.len, D.mesh, D.x, D.y)
    return mindet
end

"""
Harmonic extension of the nodal interface data `g` on the current geometry.
Returns (u, F) with u the full nodal vector and F = (K u)[interface] the
variational flux functional, F_i = int_Omega grad u_h . grad psi_i.
"""
function harmonic_extension(D::DirichletSolver, g::AbstractVector)
    m = D.mesh
    bnd, int = m.boundary, m.interior
    KII = D.K[int, int]
    KIB = D.K[int, bnd]
    rhs = -(KIB * g)
    Asym = Symmetric(KII, :U)
    if D.factor === nothing
        D.factor = cholesky(Asym)
    else
        cholesky!(D.factor, Asym)
    end
    uI = D.factor \ rhs
    u = Vector{Float64}(undef, nnodes(m))
    u[bnd] .= g
    u[int] .= uI
    F = (D.K * u)[bnd]
    return u, F
end

"Discrete Dirichlet-to-Neumann map Lambda_h g = M_Gamma^{-1} (K u_h)|_Gamma (electrolyte-exterior normal)."
function dtn_apply(D::DirichletSolver, g::AbstractVector)
    u, F = harmonic_extension(D, g)
    MΓ = boundary_mass(D.len, ones(D.mesh.nt))
    return MΓ \ F, u, F
end

"""
First-order corrector data of the discrete problem at the geometry R:
`gamma1 = -M_a^{-1} F` (trace corrector), `u0` (discrete harmonic extension of
the nodal phi_eq), `F` (flux functional), `J0 = 1/2 u0'K u0`, `J1 = -1/2 F'M_a^{-1}F`.
These are exactly the kappa-derivatives at kappa = 0 of the discrete Robin problem.
"""
function corrector_data(D::DirichletSolver, P::Params, mat::Material, R::AbstractVector)
    set_geometry!(D, R)
    nt = D.mesh.nt
    eqb = [mat.phieq(D.x[j], D.y[j]) for j in 1:nt]
    i0b = [mat.i0(D.x[j], D.y[j]) for j in 1:nt]
    u0, F = harmonic_extension(D, eqb)
    Ma = boundary_mass(D.len, (P.A2 + P.A1) .* i0b)
    gamma1 = -(Ma \ F)
    J0 = 0.5 * dot(u0, D.K * u0)
    J1 = 0.5 * dot(F, gamma1)          # = -1/2 F' M_a^{-1} F
    return (gamma1=gamma1, u0=u0, F=F, J0=J0, J1=J1, eqb=eqb, i0b=i0b, Ma=Ma)
end

"Naive one-sided nodal differencing of the normal derivative at the interface (for comparison only)."
function dtn_nodal_differencing(D::DirichletSolver, u::AbstractVector)
    m = D.mesh; nt = m.nt
    out = Vector{Float64}(undef, nt)
    @inbounds for j in 1:nt
        r0 = hypot(D.x[j], D.y[j]); r1 = hypot(D.x[nt+j], D.y[nt+j])
        out[j] = -(u[nt+j] - u[j]) / (r1 - r0)     # nu = -e_r on the inner circle
    end
    return out
end
