struct MPASSWModel{S, NE, TI, TF, TZ, TB, TFV, NER} <: SWModel
    mesh::VoronoiMesh{S, NE, TI, TF, TZ}
    constants::@NamedTuple{
        H::TF,
        g::TF,
        b::TB,
        f::TFV
    }
    diagnostics::@NamedTuple{
        vort_v::Vector{TF},
        uh_e::Vector{TF},
        vort_e::Vector{TF},
        k_c::Vector{TF},
        k_v::Vector{TF}
    }
    operators::@NamedTuple{
        div_c::DivAtCell{NE, TI, TF},
        tangent_reconstruction::TangentialVelocityReconstructionThuburn{NER, TI, TF},
        curl_v::CurlAtVertex{TI, TF},
        cell_to_edge::CellToEdgeMean{TI},
        vertex_to_edge::VertexToEdgeMean{TI},
        k_reconstruction::CellKineticEnergyMPAS{NE, TI, TF}
    }
    configuration::Dict{String, Any}
end

function MPASSWModel(m::AbstractVoronoiMesh, constants, conf)
    diag = create_diag_fields_mpas(m)
    op = create_operators_mpas(m)
    return MPASSWModel(m, constants, diag, op, conf)
end

function (model::MPASSWModel)((rhs_h, rhs_u), (h_c, u_e))
    return shallow_water_mpas_rhs!(rhs_h, rhs_u, h_c, u_e,
                                    model.mesh,
                                    model.constants, model.diagnostics, model.operators)
end

function create_diag_fields_mpas(m::VoronoiMesh{S, NE, TI, TF}) where {S, NE, TI, TF}
    vort_v = zeros(TF, m.vertices.n)
    uh_e = zeros(TF, m.edges.n)
    vort_e = zeros(TF, m.edges.n)
    k_c = zeros(TF, m.cells.n)
    k_v = zeros(TF, m.vertices.n)

    return (
        vort_v = vort_v,
        uh_e = uh_e,
        vort_e = vort_e,
        k_c = k_c,
        k_v = k_v,
    )
end

function create_operators_mpas(m::VoronoiMesh{S, NE, TI, TF}) where {S, NE, TI, TF}
    div_c = DivAtCell(m)
    tangent_reconstruction = TangentialVelocityReconstructionThuburn(m)
    curl_v = CurlAtVertex(m)
    cell_to_edge = CellToEdgeMean(m)
    vertex_to_edge = VertexToEdgeMean(m)
    k_reconstruction = CellKineticEnergyMPAS(m)

    return (
        div_c = div_c,
        tangent_reconstruction = tangent_reconstruction,
        curl_v = curl_v,
        cell_to_edge = cell_to_edge,
        vertex_to_edge = vertex_to_edge,
        k_reconstruction = k_reconstruction,
    )
end

