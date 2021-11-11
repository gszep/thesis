using DifferentialEquations
using ForwardDiff,LinearAlgebra
using CairoMakie,Colors,InvertedIndices

CairoMakie.activate!()

function fixed_point(r::AbstractVector, r₀::AbstractVector; stable::Bool=true, α::Real=1)
    return -ForwardDiff.gradient( r -> exp(-α*norm(r-r₀)) * ( stable ? -1 : 1 ), r )
end

function limit_cycle(r::AbstractVector, f::Function; stable::Bool=true, ω::Real=1)

    normal = ForwardDiff.gradient(r->norm(f(r)),r)
    normal /= norm(normal)

    ∂f = ForwardDiff.jacobian( typeof(f(r))<:Number ? r->[f(r)] : f, r )
	tangent = [ (-1)^(i+1) * det(∂f[:,Not(i)]) for i ∈ 1:length(r) ]
    tangent /= norm(tangent)

    F = norm(f(r))
    return - ( stable ? 1 : -1 ) * ( tangent*ω*exp(-F) + normal*F*exp(-F) )
end

function limit_cycle(r::Point2{T}, f::Function; stable::Bool=true, ω::Real=1) where T<:Real
    ∇F, F = ForwardDiff.gradient(r->abs(f(r)),r), abs(f(r))

    normal = ∇F / norm(∇F)
    tangent = Point2{T}(-∇F[2],∇F[1]) / norm(∇F)

    return - ( stable ? 1 : -1 ) * ( tangent*ω*sign(f(r))*exp(-F) + normal*F*exp(-F) )
end

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

figure = Figure(resolution = 72 .* (7,7) )
ax = Axis(figure[1:2,1:3],
    ylabel = L"fixed points $u^*$",
    xlabel = L"control condition $p$"
)

x = range(-2,2,length=50)
lines!( ax, map(x->x^2,x), x, linewidth=3, color=:darkblue)
hlines!( ax, [0], linewidth=3, color=:lightblue, xmax=1, xmin=0.5)
hlines!( ax, [0], linewidth=3, color=:darkblue, xmin=0, xmax=0.5)
xlims!(-3,3)

ax = Axis(figure[3,1],
    ylabel = "",
    xlabel = ""
)

streamplot!( ax, x-> fixed_point(x,[0,0]),
    -2..2, -2..2, density=0.5, linewidth=3, colormap=[:lightblue], stepsize=0.1)
scatter!( ax, [0], [0], markersize=8, color=:darkblue)

ax = Axis(figure[3,2],
    ylabel = "",
    xlabel = ""
)

streamplot!( ax, x-> 0.2*fixed_point(x,[0,1/2])+ 0.2*fixed_point(x,[0,-1/2]),
    -2..2, -2..2, density=0.5, linewidth=3, colormap=[:lightblue], stepsize=0.05)
scatter!( ax, [0], [0], markersize=8, color=:gold)

ax = Axis(figure[3,3],
    ylabel = "",
    xlabel = ""
)

streamplot!( ax, x-> fixed_point(x,[0,-1])+ fixed_point(x,[0,1]),
    -2..2, -2..2, density=0.5, linewidth=3, colormap=[:lightblue], stepsize=0.1)
scatter!( ax, [0], [0], markersize=8, color=:lightblue)
scatter!( ax, [0], [-1], markersize=8, color=:darkblue)
scatter!( ax, [0], [1], markersize=8, color=:darkblue)

figure