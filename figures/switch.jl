using DifferentialEquations
using ForwardDiff,LinearAlgebra
using CairoMakie,Colors,InvertedIndices

CairoMakie.activate!()

set_theme!(Theme(
    Axis = (
        xgridvisible  = false,
		ygridvisible  = false,
        # xticksvisible  = false,
		# yticksvisible  = false,
        # xticklabelsvisible  = false,
        # yticklabelsvisible  = false,
    )
))

figure = Figure(resolution = 72 .* (12,6) )
ax = Axis(figure[1,1],
    ylabel = L"$\beta$-galactosidase, Gal",
    xlabel = L"(allo)lactose, $p$"
)

streamplot!( ax, x-> Point2( -x[1]*x[2], asinh(3x[1]-6)+3-x[2] ),
    -0.1..4, -0.1..6, density=0.5, linewidth=1, colormap=[:gray], stepsize=0.1,
    arrow_size=12)

lines!( ax, -0.1..4, x->asinh(3x-6)+3, linewidth=3, color=:darkblue)
vlines!( ax, [0], linewidth=3, color=:darkblue)
scatter!( ax, [0], [1/2], markersize=10, color=:darkblue)


ax = Axis(figure[1,2],
    ylabel = "",
    xlabel = L"(allo)lactose, $p$"
)

streamplot!( ax, x-> Point2( (-1/2-1.3x[1]+x[2]), asinh(3x[1]-6)+3-x[2] ),
    -0.1..4, -0.1..6, density=0.5, linewidth=1, colormap=[:gray], stepsize=0.1,
    arrow_size=12)

lines!( ax, -0.1..4, x->asinh(3x-6)+3, linewidth=3, color=:darkblue)
lines!( ax, -0.1..4, x->1/2+1.3x, linewidth=3, color=:darkblue)

scatter!( ax, [0], [1/2], markersize=10, color=:darkblue)
scatter!( ax, [2.05], [3.15], markersize=10, color=:lightblue)
scatter!( ax, [3.7], [5.3], markersize=10, color=:darkblue)

figure
save("figures/switch.pdf",figure)