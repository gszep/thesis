using DifferentialEquations
using ForwardDiff,LinearAlgebra
using CairoMakie,Colors,InvertedIndices

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

figure = Figure(resolution = 72 .* (8,8) )
ax = Axis(figure[1,1],
    ylabel = "",
    xlabel = L"vector field $F_\theta$"
)

streamplot!( ax, x-> limit_cycle(x,r->norm(r+[2,-2])-1) + 5*fixed_point(x,[1,-1]),
    -5..2, -2..4, density=1, linewidth=3, colormap=[:lightblue])

scatter!( ax, [1/15], [1.2], markersize=8, color=:lightblue)
scatter!( ax, [1], [-1], markersize=8, color=:darkblue)
figure


using GLMakie
using DataStructures: CircularBuffer
GLMakie.activate!()

figure = Figure(resolution = 72 .* (8,8) )

x,y = -5..2,-2..4
minima, maxima = minimum.([x,y]), maximum.([x,y])
ax = Axis(figure[1,1],

    ylabel = "",
    xlabel = L"vector field $F_\theta$",

    xlims = (minima[1],maxima[1]),
    ylims = (minima[2],maxima[2])
)

play = Button(figure[2,1]; label="Play", tellwidth=false)
point,circle = select_point(ax.scene), select_line(ax.scene)
deactivate_interaction!(ax,:rectanglezoom)

point[] = [1,-1]
circle[] = [[-2,2],[-2,1]]

field = @lift( x-> limit_cycle(x,r->norm(r-first($circle))-norm(first($circle)-last($circle))) + fixed_point(x,$point) )
streamplot!( ax, field, x, y, density=1, linewidth=1, colormap=[:lightblue], arrow_size = 0)

scatter!( ax, point, markersize=8, color=:darkblue)
scatter!( ax, @lift(first($circle)), markersize=8, color=:lightblue)

tail, ensembleSize = 200, 200
ensemble = Vector{Observable{CircularBuffer{Point2{Float32}}}}(undef,ensembleSize)

for i ∈ 1:ensembleSize
    trajectory = CircularBuffer{Point2{Float32}}(tail)

    fill!(trajectory,(maxima-minima).*rand(2)+minima)
    ensemble[i] = Observable(trajectory)
end

color = map(i->RGBA{Float32}(to_color(:darkblue).r,to_color(:darkblue).g,to_color(:darkblue).b,(i/tail)^2),1:tail);
map( trajectory -> lines!(ax, trajectory; color=color, linewidth=3), ensemble)

ds = 0.01
playing = Observable(false)
on(play.clicks) do clicks; playing[] = !playing[]; end

on(play.clicks) do clicks
    @async while playing[]
        isopen(figure.scene) || break

        @async for trajectory ∈ ensemble

            xt = last(trajectory[])
            u = (maxima-minima).*rand(2)+minima

            rand() < 1e-3 ? push!.( Ref(trajectory[]), fill(Point2{Float32}(u),tail) ) : push!( trajectory[], xt + field[](xt)*ds )
            trajectory[] = trajectory[]
        end
        sleep(0.001)
    end
end

figure