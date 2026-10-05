function _get_constants(d::Dict,
                        m::AbstractVoronoiMesh,
                        ::Val{Symbol("zonal geostrophic balance")})

    local _g, _H, _b, _f

    _b = ConstArray(Zero())

    di = d["initial_condition"]

    if di isa String

        @assert lowercase(di::String) == "zonal geostrophic balance"

        _g = g
        _H = 1e4
        _f = ConstArray(2*Ωₑ)

    else#if isa Dict

        @assert lowercase(di["name"]::String) == "zonal geostrophic balance"

        _g = haskey(di, "g") ? Float64(di["g"]::Number)::Float64 : g
        _H = haskey(di, "H") ? Float64(di["H"]::Number)::Float64 : 1e4
        _f = haskey(di, "f") ? ConstArray(Float64(di["f"]::Number)::Float64) : ConstArray(2*Ωₑ)
    end

    return (H = _H, g = _g, b = _b, f = _f)
end

function _create_initial_condition(d::Dict,
                                   mesh::AbstractVoronoiMesh{false},
                                   ::Val{Symbol("zonal geostrophic balance")})

    local _u0, _g, _f

    di = d["initial_condition"]

    if di isa String # Using default values

        @assert lowercase(di::String) == "zonal geostrophic balance"

        _u0 = 50.0
        _g = g
        _f = 2*Ωₑ

    else#if isa Dict

        @assert lowercase(di["name"]::String) == "zonal geostrophic balance"

        _u0 = haskey(di, "u0") ? Float64(di["u0"]::Number)::Float64 : 50.0
        _g = haskey(di, "g") ? Float64(di["g"]::Number)::Float64 : g
        _f = haskey(di, "f") ? Float64(di["f"]::Number)::Float64 : 2*Ωₑ
    end

    return init_cond_zonal_geostrophic_balance(mesh; u0 = _u0, g = _g, f = _f)
end

function init_cond_zonal_geostrophic_balance(m::AbstractVoronoiMesh{false}; u0 = 50.0, g = g, f = 2*Ωₑ)
    h_c = zeros(m.cells.n)
    u_e = zeros(m.edges.n)
    return init_cond_zonal_geostrophic_balance!(h_c, u_e, m, u0=u0, g=g, f=f)
end

function init_cond_zonal_geostrophic_balance!(h_c, u_e, m::AbstractVoronoiMesh{false}; u0 = 50.0, g = g, f = 2*Ωₑ)

    yp = m.y_period

    ω = 2*pi/yp

    h_c_parameter = (f * u0 / (ω*g))

    c_y = m.cells.position.y

    @batch for i in eachindex(h_c)
        @inbounds begin
            h_c[i] = h_c_parameter * cos(ω*c_y[i])
        end
    end

    nₑₓ = m.edges.normal.x
    e_y = m.edges.position.y

    @batch for i in eachindex(u_e)
        @inbounds begin
            u_e[i] = (u0 * sin(ω*e_y[i])) * nₑₓ[i]
        end
    end

    return (h_c, u_e)
end

