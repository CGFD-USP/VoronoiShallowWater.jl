function time_step_RK4!(yⁿ⁺¹, yⁿ, compute_rhs!::Function, dt::Real, rhs, yˢ)

    copy_variables!(yⁿ⁺¹, yⁿ)

    compute_rhs!(rhs, yⁿ) # k₁

    # yˢ = (dt/2)*k + yⁿ # yˢ₂
    # yⁿ⁺¹ = (dt/6)*k + yⁿ⁺¹
    stage_computation_diagonally_explicit_RK!(yˢ, yⁿ⁺¹, rhs, yⁿ, dt, 1/2, 1/6)

    compute_rhs!(rhs, yˢ) # k₂

    # yˢ = (dt/2)*k + yⁿ # yˢ₃
    # yⁿ⁺¹ = (dt/3)*k + yⁿ⁺¹
    stage_computation_diagonally_explicit_RK!(yˢ, yⁿ⁺¹, rhs, yⁿ, dt, 1/2, 1/3)

    compute_rhs!(rhs, yˢ) # k₃

    # yˢ = dt*k + yⁿ # yˢ₄
    # yⁿ⁺¹ = (dt/3)*k + yⁿ⁺¹
    stage_computation_diagonally_explicit_RK!(yˢ, yⁿ⁺¹, rhs, yⁿ, dt, 1.0, 1/3)

    compute_rhs!(rhs, yˢ) # k₄

    # yⁿ⁺¹ = (dt/6)*k + yⁿ⁺¹ 
    result_computation_diagonally_explicit_RK!(yⁿ⁺¹, rhs, dt, 1/6)

    return yⁿ⁺¹
end

function stage_computation_diagonally_explicit_RK!(y_s::AbstractArray, accumulator::AbstractArray, f::AbstractArray, y0, dt::Real, a::Real, b::Real)
    cstage = dt*a
    cacc = dt*b
    @batch for i in eachindex(accumulator)
        @inbounds begin
            fi = f[i]
            y_s[i] = muladd(cstage, fi, y0[i])
            accumulator[i] = muladd(cacc, fi, accumulator[i])
        end
    end
    return (y_s, accumulator)
end

stage_computation_diagonally_explicit_RK!(y_s::NTuple{1}, acc::NTuple{1}, f::NTuple{1}, y0::NTuple{1}, dt::Real, a::Real, b::Real) =
   stage_computation_diagonally_explicit_RK!(y_s[1], acc[1], f[1], y0[1], dt, a, b)

function stage_computation_diagonally_explicit_RK!(y_s::Tuple, acc::Tuple, f::Tuple, y0::Tuple, dt::Real, a::Real, b::Real)
   stage_computation_diagonally_explicit_RK!(y_s[1], acc[1], f[1], y0[1], dt, a, b)
   stage_computation_diagonally_explicit_RK!(Base.tail(y_s), Base.tail(acc), Base.tail(f), Base.tail(y0), dt, a, b)
   return (y_s, acc)
end

function result_computation_diagonally_explicit_RK!(accumulator::AbstractArray, f::AbstractArray, dt::Real, b::Real)
    cacc = dt*b
    @batch for i in eachindex(accumulator)
        @inbounds begin
            accumulator[i] = muladd(cacc, f[i], accumulator[i])
        end
    end
    return accumulator
end

result_computation_diagonally_explicit_RK!(acc::NTuple{1}, f::NTuple{1}, dt::Real, b::Real) =
   result_computation_diagonally_explicit_RK!(acc[1], f[1], dt, b)

function result_computation_diagonally_explicit_RK!(acc::Tuple, f::Tuple, dt::Real, b::Real)
   result_computation_diagonally_explicit_RK!(acc[1], f[1], dt, b)
   result_computation_diagonally_explicit_RK!(Base.tail(acc), Base.tail(f), dt, b)
   return acc
end

struct RK4{T}
    rhs::T
    stage::T
end

RK4(y0) = RK4(similar_variables(y0), similar_variables(y0))

function (ti::RK4{T})(yⁿ⁺¹::T, yⁿ::T, rhs_function::Function, dt::Real) where {T}
    return time_step_RK4!(yⁿ⁺¹, yⁿ, rhs_function, dt, ti.rhs, ti.stage)
end
