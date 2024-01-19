struct Pairs
    u₋::Float64
    f₋::Float64
    u₀::Float64
    f₀::Float64
    u₊::Float64
    f₊::Float64
end

struct Parabola
    a    ::Float64
    slope::Float64
end

function parabola_parameters(pairs::Pairs)
    Δf = pairs.f₊ - pairs.f₋
    Δu = pairs.u₊ - pairs.u₋
    f̄₀ = pairs.f₋ + Δf/Δu * (pairs.u₀ - pairs.u₋)
    a  = (pairs.f₀-f̄₀) / ((pairs.u₀ - pairs.u₋)*(pairs.u₊ - pairs.u₀))

    return Parabola(a, Δf/Δu)
end
    
function parabola_from_pairs(u, pairs::Pairs)
    parab = parabola_parameters(pairs)

    # This is the parabola in dependence of u
    return pairs.f₋ + parab.slope * (u - pairs.u₋) + parab.a * ((u - pairs.u₋)*(pairs.u₊ - u))
end

function parabola_max_interval(pairs::Pairs)
    parab = parabola_parameters(pairs)

    if parab.a <= 0 # in that case, convex parabola, maximum is at boundaries
        uₘₐₓ = pairs.f₊ >= pairs.f₋ ? pairs.u₊ : pairs.u₋
    else                        # concave parabola
        # first compute parabola maximum (which exists for a > 0)
        uₘₐₓ = 1/(2 * parab.a)*parab.slope + 1/2*(pairs.u₊ + pairs.u₋)
        # clip to interval:
        uₘₐₓ = min(max(pairs.u₋, uₘₐₓ), pairs.u₊)
    end

    return uₘₐₓ
end
