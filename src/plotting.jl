using PyPlot
using Unitful.DefaultSymbols

function plot_system(positions_and_elems, initial_beam, min_pos, max_pos, steps=500)
    positions = linspace(min_pos, max_pos, steps) |> collect
    radii = map(p -> radius(propagate(expand_system(positions_and_elems, p), initial_beam)), positions)
    plot(map(p -> p / mm |> Float64, positions), map(r -> r / mm |> Float64, radii))
end

export plot_system
