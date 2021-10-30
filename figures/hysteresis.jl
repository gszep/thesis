using Pkg,Distributed,DifferentialEquations,Parameters,Polyhedra,StatsBase
addprocs(4;exeflags="--project=$(Pkg.project().path)")

@everywhere using DifferentialEquations
using CairoMakie,Colors

######################################################## model
@everywhere function f(du,u,parameters,time)
	@unpack t,t′,τ = parameters
	p,p′ = time < τ ? t : t′
	du[1] = (p′-p)/√2 + (p′+p)/√2 * u[1] - u[1]^3
	du[2] = u[1] - u[2]
end

@everywhere function σ(du,u,p,t)
	du .= fill(0.5,length(u))
end

#########################################################
#########################################################
#########################################################

set_theme!(Theme(
	fontsize=48,
    Axis = (
        xgridvisible  = false,
		ygridvisible  = false,
        xticksvisible  = false,
		yticksvisible  = false,
        xticklabelsvisible  = false,
        yticklabelsvisible  = false,
    )
))

figure = Figure(resolution = 72 .* (36,12) )
for (i,condition) ∈ enumerate([[2,2],[0,2],[2,0]])
	ax = Axis(figure[i,1],
		ylabel = i == 2 ? L"states $u(t)$" : "",
		xlabel = i == 3 ? L"$t$" : ""
	)
	for (prime,color) ∈ [([2,-1/2],RGBA(1.0,0.753,0.0,0.1)), ([-1/2,2],RGBA(0.0,0.69,0.941,0.1))]
		ensemble = solve( EnsembleProblem(
				SDEProblem(f,σ,randn(2),(0.0,5.0),(t=prime,t′=condition,τ=2.5)),
			prob_func = (x,_,_) -> remake(x,u0=randn(length(x.u0)))),
			SRIW1(), EnsembleDistributed(),trajectories=50
		)

		for trajectory ∈ ensemble
			lines!(ax,trajectory.t,map(u->u[1],trajectory.u),linewidth=5,color=color)
		end
	end
end

ax = Axis(figure[1:3,2],
    ylabel = L"p\prime",
    xlabel = L"p"
)

Δp = 25
Δ = Δp/9

U,V = Matrix{Float64}(undef,10,10),Matrix{Float64}(undef,10,10)
for condition ∈ Iterators.product(range(-1/2,2,length=10),range(-1/2,2,length=10))
	for (prime,color) ∈ [([2,-1/2],RGBA(1.0,0.753,0.0,0.1)), ([-1/2,2],RGBA(0.0,0.69,0.941,0.1))]
		ensemble = solve( EnsembleProblem(
				SDEProblem(f,σ,randn(2),(0.0,10.0),(t=prime,t′=condition,τ=5)),
			prob_func = (x,_,_) -> remake(x,u0=randn(length(x.u0)))),
			SRIW1(), EnsembleDistributed(),trajectories=200
		)

		density = StatsBase.normalize(fit(Histogram, (map(x->x.u[end][1]+Δp*condition[1],ensemble),map(x->x.u[end][2]+Δp*condition[2],ensemble)),nbins=10), mode=:density)
		contourf!(ax,density,colormap=[RGBA(color.r,color.g,color.b,0.5)],levels=[0.1,0.5,1.0],mode=:relative)

		x,y = Δp.*condition
		poly!(ax,[Point(x-Δ,y-Δ),Point(x-Δ,y+Δ),Point(x+Δ,y+Δ),Point(x+Δ,y-Δ)],
			color=RGBA(0,0,0,0),strokecolor=:black,strokewidth=1)
	end
end

ax = Axis(figure[1:3,3],
    ylabel = L"p\prime",
    xlabel = L"p"
)

figure

save("figures/hysteresis.pdf",figure)