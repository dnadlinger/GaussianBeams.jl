using SpecialFunctions

clipped_intensity(edge_pos, beam_radius) = 1 / 2 * (1 - erf(sqrt(2) * edge_pos / beam_radius))

export clipped_intensity
