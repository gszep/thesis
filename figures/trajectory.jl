using DifferentialEquations
using ForwardDiff,LinearAlgebra
using CairoMakie,Colors,InvertedIndices

function fixed_point(r::Vector; r₀::AbstractVector=zeros(length(r)), stable::Bool=true)
    return -ForwardDiff.gradient( r -> exp(-norm(r-r₀)) * ( stable ? -1 : 1 ), r )
end

function limit_cycle(r::Vector; f::Function, stable::Bool=true)
    U(r) = exp(-norm(f(r))) * ( stable ? -1 : 1 )
    normal = ForwardDiff.gradient(U,r)

    ∂f = ForwardDiff.jacobian( typeof(f(r))<:Number ? r->[f(r)] : f, r )
	tangent = similar(∂f[1,:])

	for i ∈ 1:length(r)
		tangent[i] = (-1)^(i+1) * det(∂f[:,Not(i)])
	end
    return - ( tangent/norm(tangent) + normal/norm(normal) ) * abs(U(r))
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

figure = Figure(resolution = 72 .* (8,8) )
ax = Axis(figure[1,1],
    ylabel = "",
    xlabel = L"vector field $F_\theta$"
)

streamplot!( ax, (x,y) -> Point2(
    limit_cycle([x,y]; f=r->norm(r+[2,-2])-1) + 5*fixed_point([x,y]; r₀=[1,-1])
), -5..2, -2..4, density=1, linewidth=3, colormap=[:lightblue],
arrow_size=12)

scatter!( ax, [1/10], [1.2], markersize=8, color=:lightblue)
scatter!( ax, [1], [-1], markersize=8, color=:darkblue)

figure