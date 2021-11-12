using LaTeXStrings,Plots
using Plots.Measures

using LinearAlgebra, NLsolve
using ForwardDiff,Zygote

α = 0.1
F(z::AbstractVector) = F(z[1:end-1],z[end])
function F(u::AbstractVector,p::Number)
	du = similar(u*p)
	du[1], du[2] = -u[1] + α*u[2] + u[2]*u[1]^2, p - α*u[2] - u[2]*u[1]^2
	return du
end

J(z::AbstractVector) = ForwardDiff.jacobian( z->F(z), z )
J(u::AbstractVector,p::Number) = ForwardDiff.jacobian( u->F(u,p), u )


λ(u::AbstractVector,p::Number) = eigvals(J(u,p))
function ∂λ(u::AbstractVector,p::Number)
	s = nullspace(J([u;p]))

	@assert size(s,2) == 1
	s = s[:,1] / norm(s[:,1])
	du,dp = s[1:end-1], s[end]

	∂, = Zygote.jacobian(a->real(λ(u+a*du,p+a*dp)),0.0)
	return ∂
end

function measure(u::AbstractVector,p::Number)
	λs, ∂λs = λ(u,p), ∂λ(u,p)
	x = @. real(λs) / ∂λs
	return sum( x -> 1/(1+abs(x)), x ) / length(x)
end

ps = range(0,1,step=0.0001)
u,p = randn(2), randn()
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]

layout = @layout [a;b{1.0w,0.5h}]
default(); default(grid=false,label="",margin=0mm)
plot(layout = layout, link = :x, size=(400,500), xlim=extrema(ps))

α = 0.13
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]
plot!(ps, measure.(us,ps), subplot=1, linewidth=2, color=:gold )

α = 0.1
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]
plot!(ps, measure.(us,ps), subplot=1, linewidth=2, color=:red,
	ylabel=L"\mathrm{Measure}\,\,\varphi_{\theta}(s)", ylim=(0,1), xlim=extrema(ps)
)

hline!([0],subplot=2,linewidth=1,ylim=(-1,1),color=:black,
	xlabel=L"\mathrm{arclength,}s", xmirror=true, topmargin=-5mm,
	ylabel=L"\mathrm{eigenvalues}\,\quad\lambda(s)"
)

α = 0.13
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]
for i ∈ 1:length(first(us))
	plot!(ps, map( (u,p) -> real(λ(u,p))[i], us,ps), subplot=2, linewidth=2, color=:gold )
	plot!(ps, map( (u,p) -> imag(λ(u,p))[i], us,ps), subplot=2, linestyle=:dot, linewidth=2, color=:gold )
end

α = 0.1
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]
for i ∈ 1:length(first(us))
	plot!(ps, map( (u,p) -> real(λ(u,p))[i], us,ps), subplot=2, linewidth=2, color=:red )
	plot!(ps, map( (u,p) -> imag(λ(u,p))[i], us,ps), subplot=2, linestyle=:dot, linewidth=2, color=:red )
end

plot!([NaN],[NaN], subplot=2, color=:gray, linewidth=2, label=L"\Real" )
plot!([NaN],[NaN], subplot=2, color=:gray, linestyle=:dot, linewidth=2, label=L"\Imag" )

yticks!([0,1],subplot=1)
yticks!([-2,-1,0],subplot=2)
xticks!([0,0.5,1],subplot=2)

xticks!([NaN],subplot=1)

plot!([0.16,0.16],[0,1],subplot=1,linewidth=2,color=:gray)
plot!([0.42,0.42],[0,1],subplot=1,linewidth=2,color=:gray)
plot!([0.79,0.79],[0,5],subplot=2,linewidth=2,color=:gray)

plot!([0.16,0.16],[0,5],subplot=2,linewidth=2,color=:gray)
plot!([0.79,0.79],[0,1],subplot=1,linewidth=2,color=:gray)
plot!([0.42,0.42],[0,5],subplot=2,linewidth=2,color=:gray)

scatter!([0.42],[0],subplot=2,marker=:star,color=:black,markersize=7)
scatter!([0.79],[0],subplot=2,marker=:star,color=:black,markersize=7)
scatter!([0.16],[0],subplot=2,marker=:circle,color=:black,markersize=4)

savefig("figures/hopf-measure.pdf")
