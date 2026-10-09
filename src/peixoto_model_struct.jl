
struct PeixotoSWModel{S, NE, TI, TF, TZ, TB, TFV, NER} <: SWModel
    mesh::VoronoiMesh{S, NE, TI, TF, TZ}
    constants::@NamedTuple{
        H::TF,
        g::TF,
        b::TB,
        f::TFV
    }
    diagnostics::@NamedTuple{
        pv_v::Vector{TF},
        uh_e::Vector{TF},
        pv_e::Vector{TF},
        k_c::Vector{TF}
    }
    operators::@NamedTuple{
        div_c::DivAtCell{NE, TI, TF},
        tangent_reconstruction::TangentialVelocityReconstructionPeixoto{NER, TI, TF},
        curl_v::CurlAtVertex{TI, TF},
        cell_to_vertex::CellToVertexBaricentric{TI, TF},
        cell_to_edge::CellToEdgeMean{TI},
        vertex_to_edge::VertexToEdgeMean{TI},
        k_reconstruction::CellKineticEnergyPerot{NE, TI, TF, TZ}
    }
    configuration::Dict{String, Any}
end

function PeixotoSWModel(m::AbstractVoronoiMesh, constants, conf)
    diag = create_diag_fields_trisk(m)
    op = create_operators_peixoto(m)
    return PeixotoSWModel(m, constants, diag, op, conf)
end

function create_operators_peixoto(m::VoronoiMesh{S, NE, TI, TF}) where {S, NE, TI, TF}
    div_c = DivAtCell(m)
    tangent_reconstruction = TangentialVelocityReconstructionPeixoto(m)
    curl_v = CurlAtVertex(m)
    cell_to_vertex = CellToVertexBaricentric(m)
    cell_to_edge = CellToEdgeMean(m)
    vertex_to_edge = VertexToEdgeMean(m)
    k_reconstruction = CellKineticEnergyPerot(m)

    return (
        div_c = div_c,
        tangent_reconstruction = tangent_reconstruction,
        curl_v = curl_v,
        cell_to_vertex = cell_to_vertex,
        cell_to_edge = cell_to_edge,
        vertex_to_edge = vertex_to_edge,
        k_reconstruction = k_reconstruction,
    )
end

function (model::PeixotoSWModel)((rhs_h, rhs_u), (h_c, u_e))
    return shallow_water_trisk_rhs!(rhs_h, rhs_u, h_c, u_e,
                                    model.mesh,
                                    model.constants, model.diagnostics, model.operators)
end

