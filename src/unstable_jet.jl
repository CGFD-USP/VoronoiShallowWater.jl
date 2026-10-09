# Two zonal jets of opposite sign, u = u0*sin(2π*y/y_period)^81 and v = 0, centred at y_period/4
# and 3*y_period/4, with the depth deviation in geostrophic balance, g*∂h/∂y = -f*u, from h = 0 at
# y = 0, plus two Gaussian bumps A*exp(-1000*d^2) that trigger the barotropic instability, d the
# periodic distance, in domain units, to (0.85, 0.75) and to (0.15, 0.25). Planar analogue of
# Galewsky et al. (2004), as in Peixoto & Schreiber (2019). The balanced depth is integrated with
# the trapezoidal rule on `profile_points` intervals of y and interpolated linearly, as the
# `unstable_jet` benchmark of vulpes, so that both codes start from the same state.
# Defaults: u0 = 50, f = 2*Ωₑ, H = 1e4, A = 0.01*H and profile_points = 1000.

const UNSTABLE_JET_POWER = 81
const UNSTABLE_JET_BUMP_SHARPNESS = 1000.0
const UNSTABLE_JET_BUMP_CENTRES = ((0.85, 0.75), (0.15, 0.25))

function _unstable_jet_parameters(d::Dict)

    di = d["initial_condition"]

    p = if di isa String
        @assert lowercase(di::String) == "unstable jet"
        Dict{String, Any}()
    else#if isa Dict
        @assert lowercase(di["name"]::String) == "unstable jet"
        di
    end

    _g = haskey(p, "g") ? Float64(p["g"]::Number)::Float64 : g
    _H = haskey(p, "H") ? Float64(p["H"]::Number)::Float64 : 1e4
    _f = haskey(p, "f") ? Float64(p["f"]::Number)::Float64 : 2*Ωₑ
    _u0 = haskey(p, "u0") ? Float64(p["u0"]::Number)::Float64 : 50.0
    _A = haskey(p, "amplitude") ? Float64(p["amplitude"]::Number)::Float64 : 0.01 * _H
    _n = haskey(p, "profile_points") ? (p["profile_points"]::Int) : 1000

    return (g = _g, H = _H, f = _f, u0 = _u0, A = _A, n = _n)
end

function _get_constants(d::Dict,
                        m::AbstractVoronoiMesh{false},
                        ::Val{Symbol("unstable jet")})

    p = _unstable_jet_parameters(d)

    return (H = p.H, g = p.g, b = ConstArray(Zero()), f = ConstArray(p.f))
end

function _create_initial_condition(d::Dict,
                                   mesh::AbstractVoronoiMesh{false},
                                   ::Val{Symbol("unstable jet")})

    p = _unstable_jet_parameters(d)

    return init_cond_unstable_jet(mesh; u0 = p.u0, A = p.A, n = p.n, g = p.g, f = p.f)
end

unstable_jet_velocity(yn, u0) = u0 * sin(2π * yn)^UNSTABLE_JET_POWER

# Balanced depth deviation at yn = i/n, i = 0, ..., n, by the trapezoidal rule from h = 0 at yn = 0
function unstable_jet_depth_profile(u0, n, g, f, y_period)
    h = zeros(n + 1)
    dy = 1 / n
    for i in 1:n
        u_mean = 0.5 * (unstable_jet_velocity((i - 1) * dy, u0) + unstable_jet_velocity(i * dy, u0))
        h[i + 1] = h[i] - (f / g) * u_mean * dy * y_period
    end
    return h
end

function unstable_jet_profile_value(h, yn)
    n = length(h) - 1
    t = yn * n
    i = min(floor(Int, t), n - 1)
    w = t - i
    return h[i + 1] * (1 - w) + h[i + 2] * w
end

unstable_jet_periodic_distance(a, c) = (d = abs(a - c); min(d, 1 - d))

function unstable_jet_bumps(xn, yn)
    s = zero(xn)
    for (xb, yb) in UNSTABLE_JET_BUMP_CENTRES
        d2 = unstable_jet_periodic_distance(xn, xb)^2 + unstable_jet_periodic_distance(yn, yb)^2
        s += exp(-UNSTABLE_JET_BUMP_SHARPNESS * d2)
    end
    return s
end

function init_cond_unstable_jet(m::AbstractVoronoiMesh{false}; A, u0 = 50.0, n = 1000, g = g, f = 2*Ωₑ)
    h_c = zeros(m.cells.n)
    u_e = zeros(m.edges.n)
    return init_cond_unstable_jet!(h_c, u_e, m; u0 = u0, A = A, n = n, g = g, f = f)
end

function init_cond_unstable_jet!(h_c, u_e, m::AbstractVoronoiMesh{false}; A, u0 = 50.0, n = 1000, g = g, f = 2*Ωₑ)

    xp = m.x_period
    yp = m.y_period

    h_profile = unstable_jet_depth_profile(u0, n, g, f, yp)

    c_x = m.cells.position.x
    c_y = m.cells.position.y

    @batch for i in eachindex(h_c)
        @inbounds begin
            xn = mod(c_x[i] / xp, 1)
            yn = mod(c_y[i] / yp, 1)
            h_c[i] = unstable_jet_profile_value(h_profile, yn) + A * unstable_jet_bumps(xn, yn)
        end
    end

    e_y = m.edges.position.y
    nₑₓ = m.edges.normal.x

    @batch for i in eachindex(u_e)
        @inbounds begin
            u_e[i] = unstable_jet_velocity(mod(e_y[i] / yp, 1), u0) * nₑₓ[i]
        end
    end

    return (h_c, u_e)
end
