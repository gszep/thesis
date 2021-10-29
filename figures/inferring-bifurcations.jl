using Pkg,Distributed,DifferentialEquations,BifurcationInference
addprocs(4;exeflags="--project=$(Pkg.project().path)")

@everywhere using DifferentialEquations
using CairoMakie,Colors

######################################################## model
F(z::BorderedArray,θ::AbstractVector) = F(z.u,(θ=θ,p=z.p))
@everywhere function F(u::AbstractVector,parameters::NamedTuple)
	@unpack θ,p = parameters

	f = first(u)*first(p)*first(θ)
	F = similar(u,typeof(f))

	F[1] = θ[1] + p*u[1] + θ[2]*u[1]^3
    F[2] = u[1] - u[2]
	return F
end

@everywhere function f(du,u,p,t)
	du .= F(u,p)
end

@everywhere function σ(du,u,p,t)
  du .= fill(0.1,length(u))
end

######################################################### targets and initial guess
X = StateSpace( 2, -5:0.01:5, [-3,3] )
θ = SizedVector{2}(0.1,-π/2)

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
    )
))

figure = Figure(resolution = 72 .* (8,4) )
ax = Axis(figure[1:2,2:3],
    ylabel = L"steady states, $u(t\rightarrow\infty)$",
    xlabel = L"control condition, $p$"
)

for branch ∈ deflationContinuation(F,X.roots,(θ=θ,p=-5),getParameters(X))[1:2]

	parameter = map( s -> s.z.p, branch)
	states = map(s-> s.z.u[1], branch)

	lines!( ax, parameter, states,
		color=map(s->maximum(real(s.λ))<0 ? :darkblue : :lightblue, branch) )

	scatter!( ax, parameter[map(s->s.bif,branch)], states[map(s->s.bif,branch)],
		markersize=5, color=:black)
end

vlines!(ax, Vector(X.targets),linewidth=1,color=:gray,linestyle=:dash)

ensemble = solve( EnsembleProblem(
		SDEProblem(f,σ,randn(2),(0.0,1.0),(θ=θ,p=-3)),
	prob_func = (x,_,_) -> remake(x,u0=randn(length(x.u0)))),
	SRIW1(), EnsembleDistributed(),trajectories=500
)

ax = Axis(figure[1, 1],
	title=L"p<0",
    ylabel = L"state $u(t)$",
    xlabel = L"t"
)

for trajectory ∈ ensemble
	lines!(ax,trajectory.t,map(u->u[1],trajectory.u),linewidth=1,color=RGBA(0.0,0.0,0.545,0.1))
end

ax = Axis(figure[2, 1],
    ylabel = L"u_2",
    xlabel = L"u_1"
)

streamplot!( ax, (x,y)->Point2(F([x,y],(θ=θ,p=-3))), -3..3, -3..3, density=0.5, stepsize=0.1, colormap=[RGBA(0.0,0.0,0.545,0.2)])
scatter!( ax, [0], [0], markersize=5, color=:darkblue)

ensemble = solve( EnsembleProblem(
		SDEProblem(f,σ,randn(2),(0.0,1.0),(θ=θ,p=3)),
	prob_func = (x,_,_) -> remake(x,u0=randn(length(x.u0)))),
	SRIW1(), EnsembleDistributed(),trajectories=500
)

ax = Axis(figure[1, 4],
	title=L"p>0",
    ylabel = L"state $u(t)$",
    xlabel = L"t"
)

for trajectory ∈ ensemble
	lines!(ax,trajectory.t,map(u->u[1],trajectory.u),linewidth=1,color=RGBA(0.0,0.0,0.545,0.1))
end

ax = Axis(figure[2, 4],
    ylabel = L"u_2",
    xlabel = L"u_1"
)

streamplot!( ax, (x,y)->Point2(F([x,y],(θ=θ,p=3))), -3..3, -3..3, density=0.5, stepsize=0.1, colormap=[RGBA(0.0,0.0,0.545,0.2)])
scatter!( ax, [1.3,-1.3],[2.,-1.5], markersize=5, color=:darkblue)
scatter!( ax, [0], [0], markersize=5, color=:lightblue)

figure