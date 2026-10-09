using Test

using NCDatasets
using Zeros
using TensorsLite
using VoronoiMeshes
using VoronoiOperators
using VoronoiShallowWater
using VoronoiShallowWater: parse_seconds
using TOML


@testset "Time duration parsing" begin
    @test parse_seconds("17") === 17.0
    @test parse_seconds("13:17") == 13*60 + 17.0
    @test parse_seconds("12:13:17") == 12*60*60 + 13*60 + 17.0

    @test parse_seconds("1_12:13:17") == 1*24*60*60 + 12*60*60 + 13*60 + 17.0
    @test parse_seconds("01-02_12:13:17") == 1*30*24*60*60 +  2*24*60*60 + 12*60*60 + 13*60 + 17.0

    @test parse_seconds("02-00-00_00:00:00") == 2.0 * (365*24*60*60)

    @test parse_seconds("2-03_12:30") == 2.0*(30*24*60*60) + 3*24*60*60 + 12*60 + 30
end

@inline square(x) = x*x

const config_trisk = TOML.parsefile("config-mesh.toml")

@testset "TRiSK Model" begin

    model = SWModel(config_trisk)

    @test typeof(model) === TriskSWModel{false, 6, Int32, Float64, Zero, VoronoiShallowWater.ConstArray{Zero}, VoronoiShallowWater.ConstArray{Float64}, 10}

    m = model.mesh

    h, u = create_initial_condition(config_trisk, m)

    @test length(h) == m.cells.n
    @test length(u) == m.edges.n

    rhs_h, rhs_u = VoronoiShallowWater.similar_variables((h,u))

    @test length(rhs_h) === length(h)
    @test length(rhs_u) === length(u)

    model((rhs_h, rhs_u), (h,u))

    At = m.x_period*m.y_period

    err_h = sqrt(surface_integral(m, square, rhs_h) / At)
    err_u = sqrt(surface_integral(m, square, rhs_u) / At)

    @test (err_h / model.constants.H) <= eps()

    @test (err_u / maximum(abs, u)) <= 1e-6

end

const config_peixoto = TOML.parsefile("config-mesh-peixoto.toml")

@testset "Peixoto Model" begin

    model = SWModel(config_peixoto)

    @test typeof(model) === PeixotoSWModel{false, 6, Int32, Float64, Zero, VoronoiShallowWater.ConstArray{Zero}, VoronoiShallowWater.ConstArray{Float64}, 10}

    m = model.mesh

    h, u = create_initial_condition(config_peixoto, m)

    @test length(h) == m.cells.n
    @test length(u) == m.edges.n

    rhs_h, rhs_u = VoronoiShallowWater.similar_variables((h,u))

    @test length(rhs_h) === length(h)
    @test length(rhs_u) === length(u)

    model((rhs_h, rhs_u), (h,u))

    At = m.x_period*m.y_period

    err_h = sqrt(surface_integral(m, square, rhs_h) / At)
    err_u = sqrt(surface_integral(m, square, rhs_u) / At)

    @test (err_h / model.constants.H) <= eps()

    @test (err_u / maximum(abs, u)) <= 1.5e-6

end

@testset "Stationary vortex" begin

    conf = merge(config_trisk, Dict{String, Any}("initial_condition" => "Stationary vortex"))

    model = SWModel(conf)

    m = model.mesh

    h, u = create_initial_condition(conf, m)

    @test length(h) == m.cells.n
    @test length(u) == m.edges.n

    g, f = model.constants.g, model.constants.f[1]

    R = 0.1 * m.x_period
    u0 = 0.5 * f * R

    D, r0 = VoronoiShallowWater.cos_bump_vortex(u0, R, 4, g, f)

    uθ(r) = r * VoronoiShallowWater.cos_bump_angular_velocity(r, D, r0, 4, g, f)

    @test uθ(R) ≈ u0
    @test uθ(0.99R) < u0 && uθ(1.01R) < u0

    @test r0 < m.x_period / 2
    @test -D <= minimum(h) < 0
    @test maximum(abs, u) <= u0

end

@testset "Unstable jet" begin

    conf = merge(config_trisk, Dict{String, Any}("initial_condition" => "Unstable jet"))

    model = SWModel(conf)

    m = model.mesh

    h, u = create_initial_condition(conf, m)

    @test length(h) == m.cells.n
    @test length(u) == m.edges.n

    g, f = model.constants.g, model.constants.f[1]

    profile = VoronoiShallowWater.unstable_jet_depth_profile(50.0, 1000, g, f, m.y_period)

    # the two jets carry opposite transports, so the balanced depth is periodic in y
    @test abs(profile[end] - profile[1]) <= 1e-9 * maximum(abs, profile)

    @test VoronoiShallowWater.unstable_jet_bumps(0.85, 0.75) ≈ 1
    @test VoronoiShallowWater.unstable_jet_bumps(0.15, 0.25) ≈ 1
    @test maximum(abs, u) <= 50.0

end

function rk4_rhs_test((o1, o2), (i1, i2))
    o1 .= 2 # dy1 = 2; y1 = 2t
    o2 .=  2 .* sqrt.(i2) # dy2 = 2√y2 ; y2 = t²
    return (o1, o2)
end

@testset "RK4" begin
    y0 = (rand(100), rand(150))
    ti = VoronoiShallowWater.RK4(y0)

    @test size.(ti.rhs) === size.(y0)
    @test size.(ti.stage) === size.(y0)

    dt = 1e-3

    yⁿ⁺¹ = VoronoiShallowWater.similar_variables(y0)

    ti(yⁿ⁺¹, y0, rk4_rhs_test, dt)

    @test mapreduce(isapprox, &, yⁿ⁺¹[1], (y0[1] .+ 2*dt))

    @test mapreduce(isapprox, &, yⁿ⁺¹[2], (y0[2] .+ (dt*dt) .+ ((2*dt) .* sqrt.(y0[2]))))

end

@testset "WriteOutputCallback" begin

    model = SWModel(config_trisk)

    m = model.mesh

    output_callback = VoronoiShallowWater.WriteOutputCallBack(m, config_trisk["output"])

    @test typeof(output_callback) === VoronoiShallowWater.WriteOutputCallBack{
        Float32,
        Tuple{
            VoronoiShallowWater.HOutput,
            VoronoiShallowWater.UOutput,
            VoronoiShallowWater.VelocityCellOutput{Vec2DxyArray{Float64,1}, CellVelocityReconstructionLSq2{6, Int32, Float64, Zero}},
            VoronoiShallowWater.PVOutput,
            VoronoiShallowWater.KineticEnergyCellOutput
        }
    }


    y = VoronoiShallowWater.create_initial_condition(config_trisk, m)

    file_out = (config_trisk["output"]["file_name"])::String

    isfile(file_out) && Base.Filesystem.rm(file_out)

    output_callback(y, y, model, 0.0, 1)

    @test isfile(file_out)

    ncfile = NCDataset(file_out, "r")

    @test ncfile.dim["nCells"] == m.cells.n
    @test ncfile.dim["nVertices"] == m.vertices.n
    @test ncfile.dim["nEdges"] == m.edges.n
    @test ncfile.dim["Time"] == 1

    @test ncfile["time"][1] == 0.0

    @test ncfile["H"][] === model.constants.H
    @test ncfile["g"][] === model.constants.g
    @test ncfile["f"][] === model.constants.f[1]

    @test eltype(ncfile["h"]) === output_callback.precision

    @test ncfile["h"][:,1] ≈ y[1]
    @test ncfile["u"][:,1] ≈ y[2]

    close(ncfile)

    ti = VoronoiShallowWater.RK4(y)
    yp1 = VoronoiShallowWater.similar_variables(y)

    current_dt = Ref(output_callback.next_time[])
    @show current_dt
    current_time = Ref(0.0)
    current_iteration = Ref(1)

    VoronoiShallowWater.advance_in_time!(yp1, y, model, ti, current_dt, current_time, current_iteration, current_dt[], output_callback)

    ncfile = NCDataset(file_out, "r")

    @test ncfile.dim["Time"] == 2

    close(ncfile)
    

    # Base.Filesystem.rm(file_out)

end
