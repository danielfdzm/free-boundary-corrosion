# ---------------------------------------------------------------------------
# Periodic spectral tools on the uniform theta-grid
# ---------------------------------------------------------------------------

theta_grid(N::Integer) = collect(0:N-1) .* (2π / N)

"Periodic Fourier derivative of order 1 or 2 (Nyquist mode of the odd derivative set to zero)."
function fourier_derivative(u::AbstractVector{<:Real}, order::Int=1)
    N = length(u)
    U = rfft(u)
    k = 0:(N ÷ 2)
    if order == 1
        @. U *= im * k
        iseven(N) && (U[end] = 0)
    elseif order == 2
        @. U *= -(k^2)
    else
        error("order must be 1 or 2")
    end
    return irfft(U, N)
end

"""
Sobolev norm of a periodic grid function with the convention of the paper,
||u||_{H^m}^2 = sum_{k=0}^m ||d^k u/dtheta^k||_{L^2(0,2pi)}^2, evaluated spectrally.
For m = 0 this is (int_0^{2pi} |u|^2 dtheta)^{1/2}.
"""
function sobolev_norm(u::AbstractVector{<:Real}, m::Integer)
    N = length(u)
    U = rfft(u) ./ N
    s = 0.0
    for (idx, k) in enumerate(0:(N ÷ 2))
        mult = (k == 0 || (iseven(N) && k == N ÷ 2)) ? 1.0 : 2.0
        w = sum(float(k)^(2j) for j in 0:m)
        s += mult * w * abs2(U[idx])
    end
    return sqrt(2π * s)
end

linf(u) = maximum(abs, u)

"Amplitude 2|c_n| of the mode n (0 < n < N/2) of a real grid function; |c_0| for n = 0."
function mode_amplitude(u::AbstractVector{<:Real}, n::Integer)
    N = length(u)
    U = rfft(u) ./ N
    n == 0 && return abs(U[1])
    return (iseven(N) && n == N ÷ 2) ? abs(U[n+1]) : 2 * abs(U[n+1])
end

"Cosine and sine coefficients (a_n, b_n) with u = a_0 + sum a_n cos + b_n sin."
function mode_coefficients(u::AbstractVector{<:Real}, n::Integer)
    N = length(u)
    U = rfft(u) ./ N
    c = U[n+1]
    return (2 * real(c), -2 * imag(c))
end

"Fraction of the L^2 energy of a mean-free grid function carried by the modes |n| = 1."
function translation_fraction(u::AbstractVector{<:Real})
    N = length(u)
    U = rfft(u) ./ N
    tot = sum(abs2(U[k]) * ((k == 1 || (iseven(N) && k == N ÷ 2 + 1)) ? 1 : 2) for k in 2:length(U))
    return tot == 0 ? 0.0 : 2 * abs2(U[2]) / tot
end

"Spectral interpolation of a periodic grid function from N to M points."
function resample(u::AbstractVector{<:Real}, M::Integer)
    N = length(u)
    N == M && return collect(float.(u))
    U = fft(u) ./ N
    V = zeros(ComplexF64, M)
    if M > N
        h = N ÷ 2
        V[1:h] .= U[1:h]                       # k = 0 .. N/2-1
        V[end-h+2:end] .= U[end-h+2:end]       # k = -N/2+1 .. -1
        if iseven(N)                           # split the Nyquist mode
            V[h+1] += U[h+1] / 2
            V[M-h+1] += U[h+1] / 2
        else
            V[h+1] = U[h+1]; V[M-h] = U[end-h]
        end
    else
        h = M ÷ 2
        V[1:h] .= U[1:h]
        V[end-h+2:end] .= U[end-h+2:end]
        if iseven(M)
            V[h+1] = real(U[h+1] + U[N-h+1])
        else
            V[h+1] = U[h+1]; V[M-h] = U[N-h]
        end
    end
    return real.(ifft(V)) .* M
end

"Exponential spectral filter sigma(k) = exp(-alpha (k/k_max)^p), applied in place."
function exponential_filter!(u::AbstractVector{Float64}; p::Int=36, alpha::Float64=36.0)
    N = length(u)
    U = rfft(u)
    kmax = N / 2
    for (idx, k) in enumerate(0:(N ÷ 2))
        U[idx] *= exp(-alpha * (k / kmax)^p)
    end
    u .= irfft(U, N)
    return u
end

"Observed order between consecutive parameter values: log(e_i/e_{i+1}) / log(p_i/p_{i+1})."
function observed_orders(p::AbstractVector, e::AbstractVector)
    n = length(e)
    ord = fill(NaN, n)
    for i in 1:n-1
        if e[i] > 0 && e[i+1] > 0
            ord[i] = log(e[i] / e[i+1]) / log(p[i] / p[i+1])
        end
    end
    return ord
end
