# ---------------------------------------------------------------------------
# Geometry of radial graphs r = R(theta)
# ---------------------------------------------------------------------------

"Curvature K = (R^2 + 2R_t^2 - R R_tt)/(R^2 + R_t^2)^{3/2}; equals 1/R on a circle."
function radial_curvature(R::AbstractVector{<:Real})
    R1 = fourier_derivative(R, 1)
    R2 = fourier_derivative(R, 2)
    return @. (R^2 + 2R1^2 - R * R2) / (R^2 + R1^2)^1.5
end

"Metric factor N(R) = sqrt(1 + (R_theta/R)^2) of the radial graph."
metric_factor(R::AbstractVector{<:Real}) = (R1 = fourier_derivative(R, 1); @. sqrt(1 + (R1 / R)^2))

"Maximal |R_theta/R| (the radial-graph condition degenerates when it blows up)."
graph_stretch(R::AbstractVector{<:Real}) = maximum(abs, fourier_derivative(R, 1) ./ R)

"Area of the polygon with vertices R_j e^{i theta_j}."
polygon_area(R::AbstractVector{<:Real}) = 0.5 * sin(2π / length(R)) * dot(R, circshift(R, -1))

"Area of the smooth radial graph, (1/2) int R^2 dtheta by the trapezoidal (spectral) rule."
spectral_area(R::AbstractVector{<:Real}) = π / length(R) * sum(abs2, R)

"Length of the polygon with vertices R_j e^{i theta_j}."
function polygon_perimeter(R::AbstractVector{<:Real})
    N = length(R)
    c = cos(2π / N)
    s = 0.0
    for j in 1:N
        k = j == N ? 1 : j + 1
        s += sqrt(R[j]^2 + R[k]^2 - 2 * R[j] * R[k] * c)
    end
    return s
end

"Length of the smooth radial graph, int sqrt(R^2 + R_theta^2) dtheta (trapezoidal rule)."
function spectral_perimeter(R::AbstractVector{<:Real})
    R1 = fourier_derivative(R, 1)
    return 2π / length(R) * sum(sqrt.(R .^ 2 .+ R1 .^ 2))
end

"Radial graph of the ellipse with semi-axes A (along x1) and B (along x2)."
ellipse_radius(theta, A::Real, B::Real) = A * B / hypot(B * cos(theta), A * sin(theta))

"Radial graph of the flower R = 1 + eps cos(m theta)."
flower_radius(theta, eps::Real, m::Integer) = 1 + eps * cos(m * theta)

"""
Exact solution of the parallel-curve lemma for the ellipse under homogeneous data, as a radial
graph on the grid `theta`: the inward parallel curve X_0 + c t nu_0 of the
ellipse (A cos phi, B sin phi), re-expressed as r = R(theta, t).
Valid while 1 - c t K_0 > 0, i.e. t < B^2/(c A).
"""
function parallel_ellipse_radius(theta::AbstractVector{<:Real}, t::Real, A::Real, B::Real, c::Real)
    R = similar(collect(float.(theta)))
    for (j, th) in enumerate(theta)
        # point of the parallel curve at ellipse parameter phi
        P = function (phi)
            sp, cp = sincos(phi)
            nrm = hypot(B * cp, A * sp)
            nx, ny = -B * cp / nrm, -A * sp / nrm          # nu_0 = Rot(tau), points into the solid
            return (A * cp + c * t * nx, B * sp + c * t * ny)
        end
        e1, e2 = cos(th), sin(th)
        # cross(e, P(phi)) vanishes when P(phi) lies on the ray of angle theta
        f = function (phi)
            px, py = P(phi)
            return e1 * py - e2 * px
        end
        # initial guess: ellipse parameter of the point of the ellipse on this ray
        phi0 = atan(A * sin(th), B * cos(th))
        lo, hi = phi0 - 1.0, phi0 + 1.0
        flo, fhi = f(lo), f(hi)
        # widen the bracket if necessary
        iter = 0
        while flo * fhi > 0 && iter < 10
            lo -= 0.5; hi += 0.5; flo, fhi = f(lo), f(hi); iter += 1
        end
        for _ in 1:100
            mid = 0.5 * (lo + hi)
            fm = f(mid)
            if fm * flo <= 0
                hi, fhi = mid, fm
            else
                lo, flo = mid, fm
            end
            hi - lo < 1e-15 && break
        end
        phi = 0.5 * (lo + hi)
        px, py = P(phi)
        R[j] = hypot(px, py)
    end
    return R
end
