using PyPlot
using Unitful
using Unitful.DefaultSymbols

function plot_system(positions_and_elems, initial_beam, min_pos, max_pos, steps=500)
    positions = range(min_pos, max_pos, length=steps) |> collect
    radii = map(p -> radius(propagate(expand_system(positions_and_elems, p), initial_beam)), positions)
    plot(map(p -> uconvert(NoUnits, p / mm), positions), map(r -> uconvert(NoUnits, r / mm), radii))
end

export plot_system
