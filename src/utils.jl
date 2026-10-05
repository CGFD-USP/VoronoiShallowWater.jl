"""
    ConstArray{T}

An utility type with the property that `Base.getindex(ConstArray(value), I::Vararg{Integer}) always return `value`.
"""
struct ConstArray{T}
    value::T
end

Base.@constprop :aggressive Base.getindex(a::ConstArray, _::Vararg{Integer}) = a.value

"""
    similar_variables(y::T) where T<:Union{<:Tuple, <:AbstractArray} -> T

If `y` is a tuple returns the equivalent of `similar.(y)`. If it is an `AbstractArray` returns `similar(y)`.
"""
similar_variables(a) = similar(a)
similar_variables(a::Tuple) = similar_variables.(a)

"""
    copy_variables!(x::T, y::T) where T<:Union{<:Tuple, <:AbstractArray} -> T

If `x` and `y` are tuples returns the equivalent of `copy!.(x, y)`. If it is an `AbstractArray` returns `copy!(x, y)`.
"""
function copy_variables!(a, b)
    @assert size(a) == size(b)
    mytmap!(a, identity, b)
end

function copy_variables!(a::Tuple, b::Tuple)
    @assert length(a) == length(b)
    copy_variables!.(a, b)
    return a
end

"""
    parse_seconds(s::AbstractString) -> Float64

Parse a duration string and return the equivalent duration in seconds.
The full format is `YYYY-MM-DD_hh:mm:ss` where `YYYY`, `MM`, `DD`, `hh`, `mm`, and `ss` denote
years, months, days, hours, minutes, and seconds, respectively. Larger units may be omitted,
so the following forms are also accepted:

- `DD_hh:mm:ss`
- `hh:mm:ss`
- `mm:ss`
- `ss`

Years and months are interpreted as fixed-duration units, with
- 1 year  = 365 days
- 1 month = 30 days

All components are parsed as floating-point values, so fractional seconds are supported.

The substrings on both side of a `_` character are parsed independently and, if any components
of the substring is omitted, the values are interpreted as those of the smaller units.
For example, `01_30:15` is parsed as 1 day, 30 minutes and 15 seconds.

# Examples

```julia-repl
julia> parse_seconds("1:30:00")
5400.0

julia> parse_seconds("30:00")
1800.0

julia> parse_seconds("34")
34.0

julia> parse_seconds("2-03_12:30") # 2 months + 3 days + 12 minutes + 30 seconds
5.44395e6
```
"""
function parse_seconds(s::AbstractString)
    parts = split(s, '_')

    if length(parts) == 1
        date_parts = String[]
        time_parts = split(parts[1], ':')
    elseif length(parts) == 2
        date_parts = split(parts[1], '-')
        time_parts = split(parts[2], ':')
    else
        throw(ArgumentError("invalid time format: $s"))
    end

    length(date_parts) ≤ 3 || throw(ArgumentError("invalid date format: $s"))
    length(time_parts) ≤ 3 || throw(ArgumentError("invalid time format: $s"))

    while length(date_parts) < 3
        pushfirst!(date_parts, "0")
    end

    while length(time_parts) < 3
        pushfirst!(time_parts, "0")
    end

    # Seconds corresponding to [years, months, days, hours, minutes, seconds]
    factors = (365 * 86400, 30 * 86400, 86400,
               3600, 60, 1)

    values = (
        parse(Float64, date_parts[1]),
        parse(Float64, date_parts[2]),
        parse(Float64, date_parts[3]),
        parse(Float64, time_parts[1]),
        parse(Float64, time_parts[2]),
        parse(Float64, time_parts[3]),
    )

    return sum(values .* factors)
end

"""
    format_seconds(s::Real) -> String

Format a duration in seconds as `DD_hh:mm:ss`.

The duration is decomposed into days, hours, minutes, and seconds, using 1 day = 86400 seconds.
Fractional seconds are preserved.

# Examples

```julia-repl
julia> format_seconds(5400.0)
"00_01:30:00"

julia> format_seconds(1800.0)
"00_00:30:00"

julia> format_seconds(5.44395e6)
"63_00:19:15"

julia> format_seconds(34.5)
"00_00:00:34.5"

"""
function format_seconds(s::Real)

    s ≥ 0 || throw(ArgumentError("duration must be non-negative"))

    days, remainder = divrem(s, 86400)
    hours, remainder = divrem(remainder, 3600)
    minutes, seconds = divrem(remainder, 60)

    return "$(Int(days))_$(lpad(Int(hours), 2, '0')):" * "$(lpad(Int(minutes), 2, '0')):$(seconds)"
end

function toml_string(d::Dict{String,Any})
    io = IOBuffer()
    TOML.print(io, d)
    return String(take!(io))
end

function tminmax(r::AbstractArray)
    rmin = r[1]
    rmax = r[1]
    @batch reduction=((min, rmin), (max, rmax)) for i in eachindex(r)
        @inbounds begin
            ri = r[i]
            rmin = min(rmin, ri)
            rmax = max(rmax, ri)
        end
    end
    return (rmin, rmax)
end

@inline square(x) = x*x

@inline surface_area(m::AbstractVoronoiMesh{true}) = 4*pi*m.sphere_radius^2
@inline surface_area(m::AbstractVoronoiMesh{false}) = m.x_period*m.y_period
