using DifferentialEquations
using ForwardDiff,LinearAlgebra

using GLMakie,Colors,InvertedIndices
using DataStructures: CircularBuffer
GLMakie.activate!()

∧(u::AbstractVector,v::AbstractVector) = u*v' - v*u'
function attractor(r::AbstractVector, distance::Function; stable::Bool=true, ω::Union{Function,Nothing}=nothing)

    normal = ForwardDiff.gradient(r->distance(r),r)
    tangents = cross(normal,ω(r))
    return ( stable ? -1 : 1) * ( normal/norm(normal) * distance(r) + tangents/norm(tangents) ) * exp(-abs(distance(r)))
end

set_theme!(Theme(
    Axis3 = (
        xgridvisible  = false,
		ygridvisible  = false,
		zgridvisible  = false,
        xticksvisible  = false,
		yticksvisible  = false,
		zticksvisible  = false,
        xticklabelsvisible  = false,
        yticklabelsvisible  = false,
        zticklabelsvisible  = false,
    )
))

figure = Figure(resolution = 72 .* (8,8) )

x,y,z = -3..3,-3..3,-3..3
minima, maxima = minimum.([x,y,z]), maximum.([x,y,z])
ax = Axis3(figure[1,1],

    ylabel = "",
    xlabel = "",
    zlabel = "",

    xlims = (minima[1],maxima[1]),
    ylims = (minima[2],maxima[2]),
    zlims = (minima[3],maxima[3])
)

play = Button(figure[2,1]; label="Play", tellwidth=false)
field = x-> attractor(x,x->norm(x),ω=x->[1,0,0] )
field = x-> attractor(x,x->(sqrt(x[1]^2+x[2]^2)-2)^2+x[3]^2-1,ω=x->[0,0,sqrt(x[1]^2+x[2]^2)>2 ? 1 : -1] )
field = x-> attractor(x,x->(sqrt(x[1]^2+x[2]^2)-2)^2+x[3]^2-1,ω=x->[-x[2],x[1],0] )
field = x-> attractor(x,x->(sqrt(x[1]^2+x[2]^2)-2)^2+x[3]^2-1.9^2,ω=x->[-x[2],x[1],sqrt(x[1]^2+x[2]^2)>2 ? 2 : -2] )

tail, ensembleSize = 200, 200
ensemble = Vector{Observable{CircularBuffer{Point3{Float32}}}}(undef,ensembleSize)

for i ∈ 1:ensembleSize
    trajectory = CircularBuffer{Point3{Float32}}(tail)

    fill!(trajectory,(maxima-minima).*rand(3)+minima)
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
            u = (maxima-minima).*rand(3)+minima

            rand() < 1e-3 ? push!.( Ref(trajectory[]), fill(Point3{Float32}(u),tail) ) : push!( trajectory[], xt + field(xt)*ds )
            trajectory[] = trajectory[]
        end
        sleep(0.0001)
    end
end
