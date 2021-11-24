using Pkg,Distributed,DifferentialEquations,Parameters,Polyhedra
using LinearAlgebra,StatsBase
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
		titlegap=25,
    ),
	Label = (
		tellheight=false,
		tellwidth=false
	)
))

figure = Figure(resolution = 72 .* (36,12) )
for (i,condition) ∈ enumerate([[2,2],[0,2],[2,0]])
	ax = Axis(figure[i,1],
		title = i == 1 ? L"Protocol $p(t),{p \prime}(t)$" : "",
		ylabel = i == 2 ? L"states $u(t)$" : "",
		xlabel = i == 3 ? L"time $t$" : "",
		ylabelpadding = 25,
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
	ylims!(ax,(-5,5)); xlims!(ax,(0,5))
	vlines!(ax, [2.5], linewidth=5,color=:gray,linestyle=:dash)
end

ax = Axis(figure[1:3,2],
	title = L"Steady State Distributions $P(u)$",
    ylabel = L"conditon $p\prime$",
    xlabel = L"conditon $p$"
)

Δp = 15
Δ = Δp/5

N = 6
for condition ∈ Iterators.product(range(-1/2,2,length=N),range(-1/2,2,length=N))
	for (prime,color) ∈ [([2,-1/2],RGBA(1.0,0.753,0.0,0.1)), ([-1/2,2],RGBA(0.0,0.69,0.941,0.1))]
		ensemble = solve( EnsembleProblem(
				SDEProblem(f,σ,randn(2),(0.0,10.0),(t=prime,t′=condition,τ=5)),
			prob_func = (x,_,_) -> remake(x,u0=randn(length(x.u0)))),
			SRIW1(), EnsembleDistributed(),trajectories=200
		)

		density = StatsBase.normalize(fit(Histogram, (map(x->-x.u[end][2]+Δp*condition[1],ensemble),map(x->x.u[end][1]+Δp*condition[2],ensemble)),nbins=10), mode=:density)
		contourf!(ax,density,colormap=[RGBA(color.r,color.g,color.b,0.5)],levels=[0.1,0.5,1.0],mode=:relative)
	end

	x,y = Δp.*condition
	poly!(ax,[Point(x-Δ,y-Δ),Point(x-Δ,y+Δ),Point(x+Δ,y+Δ),Point(x+Δ,y-Δ)], color=RGBA(0,0,0,0),
		strokecolor = all(condition .≈ [2,-1/2]) ? RGB(1.0,0.753,0.0) : all(condition .≈ [-1/2,2]) ? RGB(0.0,0.69,0.941) : :black,
		strokewidth = all(condition .≈ [-1/2,2]) ? 5 : all(condition .≈ [2,-1/2]) ? 5 : 1)
end

N = 24
U = Matrix{Float64}(undef,N,N)
for (j,condition) ∈ enumerate(Iterators.product(range(-1/2,2,length=N),range(-1/2,2,length=N)))

	u = Vector{Float64}(undef,2)
	for (i,(prime,color)) ∈ enumerate([([2,-1/2],RGBA(1.0,0.753,0.0,0.1)), ([-1/2,2],RGBA(0.0,0.69,0.941,0.1))])
		ensemble = solve( EnsembleProblem(
				SDEProblem(f,σ,randn(2),(0.0,10.0),(t=prime,t′=condition,τ=5)),
			prob_func = (x,_,_) -> remake(x,u0=randn(length(x.u0)))),
			SRIW1(), EnsembleDistributed(),trajectories=200
		)
		u[i] = norm(StatsBase.cov(map(x->-x.u[end],ensemble)))
	end
	U[j] = max(u[1],u[2])
end

ax = Axis(figure[1:3,3],
	title = L"Population Separation ${\Delta}u(p,{p \prime})$",
    ylabel = L"conditon $p\prime$",
    xlabel = L"conditon $p$"
)

contourf!(ax,range(-1/2,2,length=N),range(-1/2,2,length=N),U,colormap=:grayC)
contour!(ax,range(-1/2,2,length=N),range(-1/2,2,length=N),U, levels=[0.55], mode=:relative, linewidth=5, color=:black )
scatter!( ax, [0.13], [0.13], markersize=25, color=:white, strokecolor=:black, strokewidth=5 )

Label(figure[1,1], L"p \ll {p \prime}",
	color=RGB(0.0,0.69,0.941),
	height=-100, width=700
)

Label(figure[1,1], L"p \gg {p \prime}",
	color=RGB(1.0,0.753,0.0),
	height=200, width=700
)

Label(figure[1,1], L"p \approx {p \prime}",
	height=200, width=-50
)

Label(figure[2,1], L"p \ll {p \prime}",
	color=RGB(0.0,0.69,0.941),
	height=200, width=-50
)

Label(figure[3,1], L"p \gg {p \prime}",
	color=RGB(1.0,0.753,0.0),
	height=200, width=-50
)
figure