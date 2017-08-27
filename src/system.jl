zero_like(x) = x - x

type Beam
    """Complex beam parameter."""
    q

    """Wavelength."""
    λ
end

beam(waist_radius, offset, λ) = Beam(-offset + 1im * rayleigh_range(waist_radius, λ), λ)

rayleigh_range(b::Beam) = imag(b.q)
rayleigh_range(waist_radius, λ) = π * waist_radius^2 / λ

waist_offset(b::Beam) = -real(b.q)
waist_radius(b::Beam) = sqrt(rayleigh_range(b) * b.λ / π)

radius(b::Beam) = sqrt(-b.λ / (π * imag(1 / b.q)))
curvature(b::Beam) = 1 / real(1 / b.q)

abstract type Element end

struct Lens <: Element
    focal_length
end

struct FreeSpace <: Element
    distance
end

function expand_system(positions_and_elems, up_to_pos=nothing)
    last_pos = zero_like(up_to_pos == nothing ? positions_and_elems[1][1] : up_to_pos)
    system = Element[]
    for (p, e) in sort(positions_and_elems, by=a -> a[1])
        if up_to_pos != nothing && p > up_to_pos
            break
        end
        delta_pos = p - last_pos
        if delta_pos != zero_like(p)
            push!(system, FreeSpace(delta_pos))
        end
        last_pos = p

        push!(system, e)
    end
    if up_to_pos != nothing && last_pos < up_to_pos
        push!(system, FreeSpace(up_to_pos - last_pos))
    end
    system
end

function apply_abcd(beam::Beam, matrix)
    Beam((matrix[1, 1] * beam.q + matrix[1, 2]) /
        (matrix[2, 1] * beam.q + matrix[2, 2]), beam.λ)
end

function apply(beam::Beam, elem::Lens)
    apply_abcd(beam, [1 zero_like(elem.focal_length); (-1 / elem.focal_length) 1])
end

function apply(beam::Beam, elem::FreeSpace)
    apply_abcd(beam, [1 elem.distance; zero_like(1 / elem.distance) 1])
end

function propagate{E<:Element}(elements::Array{E, 1}, beam::Beam)
    reduce(apply, beam, elements)
end

export Beam, beam, Element, Lens, FreeSpace, expand_system, propagate
