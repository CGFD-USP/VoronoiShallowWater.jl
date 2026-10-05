abstract type SWModel<:Function end

function SWModel(d::Dict)

    !haskey(d, "grid_file") && @error "No `grid_file` provided in the simulation configuration. Aborting simulation."

    m = VoronoiMesh(d["grid_file"]::String)

    if haskey(d, "grid_scale_factor")
        factor = Float64(d["grid_scale_factor"]::Number)::Float64
        scale!(m, factor)
    end

    model = if haskey(d, "model")
        lowercase(d["model"]::String)
    else
        @warn "No `model` specified in the simulation configuration. Using the \"TRiSK\" model."
        "trisk"
    end

    constants = get_constants(d, m) # constants given by the initial condition.

    if model == "trisk"
        return TriskSWModel(m, constants, d)
    elseif model == "peixoto"
        return PeixotoSWModel(m, constants, d)
    end
end

struct TriskSWModel{S, NE, TI, TF, TZ, TB, TFV, NER} <: SWModel
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
        tangent_reconstruction::TangentialVelocityReconstructionThuburn{NER, TI, TF},
        curl_v::CurlAtVertex{TI, TF},
        cell_to_vertex::CellToVertexArea{TI, TF},
        cell_to_edge::CellToEdgeMean{TI},
        vertex_to_edge::VertexToEdgeMean{TI},
        k_reconstruction::CellKineticEnergyRingler{NE, TI, TF}
    }
    configuration::Dict{String, Any}
end

function TriskSWModel(m::AbstractVoronoiMesh, constants, conf)
    diag = create_diag_fields_trisk(m)
    op = create_operators_trisk(m)
    return TriskSWModel(m, constants, diag, op, conf)
end

function (model::TriskSWModel)((rhs_h, rhs_u), (h_c, u_e))
    return shallow_water_trisk_rhs!(rhs_h, rhs_u, h_c, u_e,
                                    model.mesh,
                                    model.constants, model.diagnostics, model.operators)
end

function create_diag_fields_trisk(m::VoronoiMesh{S, NE, TI, TF}) where {S, NE, TI, TF}
    pv_v = zeros(TF, m.vertices.n)
    uh_e = zeros(TF, m.edges.n)
    pv_e = zeros(TF, m.edges.n)
    k_c = zeros(TF, m.cells.n)

    return (
        pv_v = pv_v,
        uh_e = uh_e,
        pv_e = pv_e,
        k_c = k_c
    )
end

function create_operators_trisk(m::VoronoiMesh{S, NE, TI, TF}) where {S, NE, TI, TF}
    div_c = DivAtCell(m)
    tangent_reconstruction = TangentialVelocityReconstructionThuburn(m)
    curl_v = CurlAtVertex(m)
    cell_to_vertex = CellToVertexArea(m)
    cell_to_edge = CellToEdgeMean(m)
    vertex_to_edge = VertexToEdgeMean(m)
    k_reconstruction = CellKineticEnergyRingler(m)

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
