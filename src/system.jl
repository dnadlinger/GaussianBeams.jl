zero_like(x) = x - x

struct Beam
    """Complex beam parameter."""
    q

    """Wavelength."""
    λ

    """Current index of refraction."""
    index
end

beam(waist_radius, offset, λ, index=1.0) =
    Beam(-offset + 1im * rayleigh_range(waist_radius, λ, index), λ, index)

rayleigh_range(b::Beam) = imag(b.q)
rayleigh_range(waist_radius, λ, index) = π * index * waist_radius^2 / λ

waist_offset(b::Beam) = -real(b.q)
waist_radius(b::Beam) = sqrt(rayleigh_range(b) * b.λ / (π * b.index))

radius(b::Beam) = sqrt(-b.λ / (π * b.index * imag(1 / b.q)))
curvature(b::Beam) = 1 / real(1 / b.q)

abstract type Element end

struct Lens <: Element
    focal_length
end

struct FreeSpace <: Element
    distance
end

struct Interface <: Element
    new_index  # New index of refraction.
end

new_index(old_index, elem::Element) = old_index
new_index(old_index, elem::Interface) = elem.new_index

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
        (matrix[2, 1] * beam.q + matrix[2, 2]), beam.λ, beam.index)
end

function apply(beam::Beam, elem::Lens)
    apply_abcd(beam, [1 zero_like(elem.focal_length); (-1 / elem.focal_length) 1])
end

function apply(beam::Beam, elem::Interface)
    new = apply_abcd(beam, [1 zero_like(beam.λ); zero_like(1 / beam.λ) (beam.index / elem.new_index)])
    Beam(new.q, new.λ, elem.new_index)
end

function apply(beam::Beam, elem::FreeSpace)
    apply_abcd(beam, [1 elem.distance; zero_like(1 / elem.distance) 1])
end

function propagate(elements::Array{E, 1}, beam::Beam) where E<:Element
    reduce(apply, elements, init=beam)
end

export Beam, beam, Element, Lens, FreeSpace, Interface, expand_system, propagate, waist_radius, waist_offset
