using CairoMakie,Colors

#########################################################
#########################################################
#########################################################

set_theme!(Theme(
    Axis = (
        xgridvisible  = false,
		ygridvisible  = false,
        xticksvisible  = false,
		yticksvisible  = false,
        xticklabelsvisible  = false,
        yticklabelsvisible  = false,
        ylabelpadding = 10,
    )
))

figure = Figure(resolution = 72 .* (9,7) )
ax = Axis(figure[1:3,1:3],
    ylabel = L"argon concentration $\theta$",
    xlabel = L"voltage $p$"
)

lines!(ax, fill(0,2), range(0,1,length=2), color=:gold)
lines!(ax, range(-1,1,length=2), fill(-0.1,2), color=:gray, linestyle=:dash)
lines!(ax, range(-1,1,length=2), fill(-0.6,2), color=:gray, linestyle=:dash)
lines!(ax, range(-1,1,length=2), fill(0.4,2), color=:gray, linestyle=:dash)

scatter!(ax,[0],[0.4], color=:gold, markersize=10)

xlims!(-1,1); ylims!(-1,0.8)


for (i,f) ∈ enumerate([x->x,x->tanh(2*x),x->tanh(x/2)])
    ax = Axis(figure[i,4],
        ylabel = i == 2 ? L"brightness $u$" : "",
        xlabel = i == 3 ? L"voltage $p$" : ""
    )

    if i ≠ 1
        parameter = -3:0.1:3
        lines!( ax, parameter, map( f , parameter ), color=:darkblue )
    else
        lines!( ax, -3:0.1:0, x->tanh(x)-1, color=:darkblue )
        lines!( ax,  0:0.1:3, x->tanh(x)+1, color=:darkblue )
    end

    if i == 1
        scatter!(ax,[0],[0], color=:gold, markersize=10)
    end
end

figure