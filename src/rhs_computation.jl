@inline negate_second_arg(_, x) = -x

function shallow_water_trisk_rhs!(rhs_h::AbstractArray, rhs_u::AbstractArray,
                                  h_c, u_e,
                                  mesh::AbstractVoronoiMesh,
                                  constants, diag, operators)

    #compute nedded diagnostics
    shallow_water_trisk_diagnostics!(diag, h_c, u_e, constants, operators)

    uh_e = diag.uh_e

    #Compute -(∇⋅(u*(h+H))) and store in `rhs_h`
    operators.div_c(rhs_h, negate_second_arg, uh_e)

    pv_e = diag.pv_e

    k_c = diag.k_c

    shallow_water_trisk_rhs_u!(rhs_u, mesh.edges.lengthDual, mesh.edges.cells, operators.tangent_reconstruction, uh_e, pv_e, k_c, h_c, constants.b, constants.g)

    return (rhs_h, rhs_u)
end

struct uh_e_helper{T} <: Function
    H::T
end

@inline (f::uh_e_helper)(u, h) = u * (f.H + h)

function shallow_water_trisk_diagnostics!(diag, h_c, u_e, constants, operators)

    uh_e = diag.uh_e

    #Initializes uh_e with u
    mytmap!(uh_e, identity, u_e)

    #Compute u * (H + h) at edge
    operators.cell_to_edge(uh_e, uh_e_helper(constants.H), h_c)

    k_c = diag.k_c

    #compute K at cell
    operators.k_reconstruction(k_c, u_e)

    compute_potential_vorticity!(diag, h_c, u_e, constants, operators)

    return diag
end

function compute_potential_vorticity!(diag, h_c, u_e, constants, operators)
    if haskey(operators, :curl_e)
        compute_potential_vorticity_from_edge!(diag, h_c, u_e, constants, operators)
    elseif haskey(operators, :curl_v)
        compute_potential_vorticity_from_vertex!(diag, h_c, u_e, constants, operators)
    end
end

struct pv_helper{T} <: Function
    H::T
end

@inline (f::pv_helper)(w, h) = w / (f.H + h)

function compute_potential_vorticity_from_vertex!(diag, h_c, u_e, constants, operators)
    #Initializes pv_v array with coriolis parameter f
    mytmap!(diag.pv_v, identity, constants.f)

    #Computes ωₐ = curl_v + f
    operators.curl_v(diag.pv_v, +, u_e)

    #Computes pv_v = ωₐ / (H + h_v)
    operators.cell_to_vertex(diag.pv_v, pv_helper(constants.H), h_c)

    #compute pv_v at edge
    operators.vertex_to_edge(diag.pv_e, diag.pv_v)
end

function compute_potential_vorticity_from_edge!(diag, h_c, u_e, constants, operators)
    #Initializes pv_v array with coriolis parameter f
    mytmap!(diag.pv_e, identity, constants.f)

    #Computes ωₐ = curl_e + f
    operators.curl_e(diag.pv_e, +, u_e)

    #Computes pv_v = ωₐ / (H + h_v)
    operators.cell_to_edge(diag.pv_e, pv_helper(constants.H), h_c)
end

function shallow_water_trisk_rhs_u!(rhs_u::AbstractArray, dc_e, cellsOnEdge, tangent_e, uh_e, pv_e, k_c, h_c, b, g::Real)

    w = tangent_e.weights
    edgesOnEdge = tangent_e.indices

    @batch for e in eachindex(rhs_u)
        @inbounds begin
            c1, c2 = cellsOnEdge[e]

            Gₑ = (g * ((h_c[c2] - h_c[c1]) + (b[c2] - b[c1])) + (k_c[c2] - k_c[c1])) / dc_e[e]

            Qₑ = local_trisk_coriolis_term((e,), w[e], edgesOnEdge[e], uh_e, pv_e) 

            rhs_u[e] = -(Gₑ + Qₑ)
        end
    end

    return rhs_u
end

