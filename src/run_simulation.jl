function run_simulation(sim::SWSimulation)

    sim.current_time[] = 0.0
    sim.current_iteration[] = 1



    sim.current_dt[] = sim.dt

    yⁿ = sim.initial_condition
    yⁿ⁺¹ = similar_variables(yⁿ)

    sim.output_callback(yⁿ, yⁿ⁺¹, sim.model, sim.current_time[], sim.current_iteration[])

    tuple_call(sim.callbacks, yⁿ, yⁿ⁺¹, sim.model, sim.current_time[], sim.current_iteration[])

    while (sim.current_time[] <= sim.run_duration[])

        advance_in_time!(yⁿ⁺¹, yⁿ, sim)
        
        yⁿ, yⁿ⁺¹ = (yⁿ⁺¹, yⁿ)
    end

    return nothing
end

function advance_in_time!(yp1, y0, sim::SWSimulation)
    advance_in_time!(
        yp1, y0, sim.model, sim.time_integration, sim.current_dt, sim.current_time,
        sim.current_iteration, sim.dt, sim.output_callback, sim.callbacks
    )
end
        
function advance_in_time!(yⁿ⁺¹, yⁿ, model, time_integration_method, current_dt::Base.RefValue{Float64}, current_time::Base.RefValue{Float64}, current_iteration::Base.RefValue{Int}, dt::Real, write_output_callback, other_callbacks=())

    time_integration_method(yⁿ⁺¹, yⁿ, model, current_dt[])

    current_time[] += current_dt[]

    current_iteration[] += 1

    #Inside this call it will be checked whether it is really time to write the output
    write_output_callback(yⁿ⁺¹, yⁿ, model, current_time[], current_iteration[])

    tuple_call(other_callbacks, yⁿ⁺¹, yⁿ, model, current_time[], current_iteration[])

    # Adjust current_dt such that output is written at the desired time
    if (current_time[] + current_dt[]) > write_output_callback.next_time[]
        current_dt[] = write_output_callback.next_time[] - current_time[]
    else
        current_dt[] = dt
    end

    return nothing
end

tuple_call(::Tuple{}, ::Vararg) = nothing

tuple_call(cb, args::Vararg) = cb(args...)

tuple_call(cbs::NTuple{1}, args::Vararg) = tuple_call(cbs[1], args...)

function tuple_call(cbs::Tuple, args::Vararg)
    tuple_call(cbs[1], args...)
    tuple_call(Base.tail(cbs), args...)
end
