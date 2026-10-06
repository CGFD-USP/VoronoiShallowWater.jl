# Cyclone in gradient-wind balance centred in the domain: depth deviation
# h = -D*cos(π*r/(2*r0))^(2n) for r < r0 and 0 beyond, with D and r0 set so that the
# azimuthal speed peaks at u0 at radius R. Exact steady state of the nonlinear equations.
# Defaults: f = 2*sqrt(g*H)/x_period (domain two deformation radii wide), R = 0.1*x_period,
# u0 = 0.5*f*R (Rossby number 0.5) and n = 4.

function _stationary_vortex_parameters(d::Dict, m::AbstractVoronoiMesh{false})

    di = d["initial_condition"]

    p = if di isa String
        @assert lowercase(di::String) == "stationary vortex"
        Dict{String, Any}()
    else#if isa Dict
        @assert lowercase(di["name"]::String) == "stationary vortex"
        di
    end

    _g = haskey(p, "g") ? Float64(p["g"]::Number)::Float64 : g
    _H = haskey(p, "H") ? Float64(p["H"]::Number)::Float64 : 1e4
    _f = haskey(p, "f") ? Float64(p["f"]::Number)::Float64 : 2 * sqrt(_g * _H) / m.x_period
    _R = (haskey(p, "radius_fraction") ? Float64(p["radius_fraction"]::Number)::Float64 : 0.1) * m.x_period
    _u0 = haskey(p, "u0") ? Float64(p["u0"]::Number)::Float64 : 0.5 * _f * _R
    _n = haskey(p, "bump_power") ? (p["bump_power"]::Int) : 4

    return (g = _g, H = _H, f = _f, R = _R, u0 = _u0, n = _n)
end

function _get_constants(d::Dict,
                        m::AbstractVoronoiMesh{false},
                        ::Val{Symbol("stationary vortex")})

    p = _stationary_vortex_parameters(d, m)

    return (H = p.H, g = p.g, b = ConstArray(Zero()), f = ConstArray(p.f))
end

function _create_initial_condition(d::Dict,
                                   mesh::AbstractVoronoiMesh{false},
                                   ::Val{Symbol("stationary vortex")})

    p = _stationary_vortex_parameters(d, mesh)

    return init_cond_stationary_vortex(mesh; u0 = p.u0, R = p.R, n = p.n, g = p.g, f = p.f)
end

# h'(r)/r of the unit-depth cos^(2n) bump at x = k*r, k = π/(2*r0)
function cos_bump_slope_over_radius(x, k, n)
    sinc_x = x > 0 ? sin(x) / x : one(x)
    return 2n * k * k * cos(x)^(2n - 1) * sinc_x
end

# D and r0 such that the azimuthal speed peaks at u0 at R. With w = u0/R and
# x = π*R/(2*r0), a zero derivative at R reads
# (2n-1)*x*tan(x) - x*cot(x) + 1 = (2w+f)/(w+f), increasing in x.
function cos_bump_vortex(u0, R, n, g, f)
    w = u0 / R
    q = (2w + f) / (w + f)
    lo, hi = 0.0, π / 2
    x = (lo + hi) / 2
    while lo < x < hi
        if (2n - 1) * x * tan(x) - x / tan(x) + 1 < q
            lo = x
        else
            hi = x
        end
        x = (lo + hi) / 2
    end
    k = x / R
    return (D = w * (w + f) / (g * cos_bump_slope_over_radius(x, k, n)), r0 = π / (2k))
end

cos_bump_depth(r, D, r0, n) = r < r0 ? -D * cos(π * r / (2r0))^(2n) : zero(r)

# u_θ/r, root of (u_θ/r)^2 + f*u_θ/r = g*h'/r, rationalized so it does not cancel where h' -> 0
function cos_bump_angular_velocity(r, D, r0, n, g, f)
    r < r0 || return zero(r)
    k = π / (2r0)
    ga = g * D * cos_bump_slope_over_radius(k * r, k, n)
    return 2ga / (f + sqrt(f * f + 4ga))
end

function init_cond_stationary_vortex(m::AbstractVoronoiMesh{false}; u0, R, n = 4, g, f)
    h_c = zeros(m.cells.n)
    u_e = zeros(m.edges.n)
    return init_cond_stationary_vortex!(h_c, u_e, m; u0 = u0, R = R, n = n, g = g, f = f)
end

function init_cond_stationary_vortex!(h_c, u_e, m::AbstractVoronoiMesh{false}; u0, R, n = 4, g, f)

    D, r0 = cos_bump_vortex(u0, R, n, g, f)

    @assert r0 < min(m.x_period, m.y_period) / 2 "vortex support must fit in the periodic domain"

    xc = m.x_period / 2
    yc = m.y_period / 2

    c_x = m.cells.position.x
    c_y = m.cells.position.y

    @batch for i in eachindex(h_c)
        @inbounds begin
            h_c[i] = cos_bump_depth(hypot(c_x[i] - xc, c_y[i] - yc), D, r0, n)
        end
    end

    e_x = m.edges.position.x
    e_y = m.edges.position.y
    nₑₓ = m.edges.normal.x
    nₑᵧ = m.edges.normal.y

    @batch for i in eachindex(u_e)
        @inbounds begin
            dx = e_x[i] - xc
            dy = e_y[i] - yc
            ω = cos_bump_angular_velocity(hypot(dx, dy), D, r0, n, g, f)
            u_e[i] = ω * (dx * nₑᵧ[i] - dy * nₑₓ[i])
        end
    end

    return (h_c, u_e)
end
