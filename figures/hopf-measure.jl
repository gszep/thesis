using LaTeXStrings,Plots
using Plots.Measures

using LinearAlgebra, NLsolve
using ForwardDiff,Zygote

α = 2
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
	re_part, = Zygote.jacobian(u->real(eigvals(J(u,p))), u)
	im_part, = Zygote.jacobian(u->imag(eigvals(J(u,p))), u)
	return re_part + im_part * 1im
end

function measure(u::AbstractVector,p::Number)
	x = λ(u,p)./∂λ(u,p)[1,:]
	return real( sum( x -> 1/(1+abs(x)), x ) / length(x) )
end

ps = range(0.75,1.25,step=0.005)
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]

layout = @layout [a;b{1.0w,0.5h}]
default(); default(grid=false,label="",margin=0mm)
plot(layout = layout, link = :x, size=(400,500), xlim=(0.75,1.25) )

α = 2.21
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]
plot!(ps, measure.(us,ps), subplot=1, linewidth=2, color=:gold )

α = 2.2
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]
plot!(ps, measure.(us,ps), subplot=1, linewidth=2, color=:red,
	ylabel=L"\mathrm{Measure}\,\,\varphi_{\theta}(s)", ylim=(0,1), xlim=(0.75,1.25)
)

hline!([0],subplot=2,linewidth=1,ylim=(-2,0.5),color=:black,
	xlabel=L"\mathrm{arclength,}s", xmirror=true, topmargin=-5mm,
	ylabel=L"\mathrm{eigenvalues}\,\quad\lambda(s)"
)

α = 2.21
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]
for i ∈ 1:length(first(us))
	plot!(ps, map( (u,p) -> real(λ(u,p))[i], us,ps), subplot=2, linestyle=:dot, linewidth=2, color=:gold )
	plot!(ps, map( (u,p) -> imag(λ(u,p))[i], us,ps), subplot=2, linewidth=2, color=:gold )
end

α = 2.2
us = [ nlsolve(u->F(u,p),u).zero for p ∈ ps ]
for i ∈ 1:length(first(us))
	plot!(ps, map( (u,p) -> real(λ(u,p))[i], us,ps), subplot=2, linestyle=:dot, linewidth=2, color=:red )
	plot!(ps, map( (u,p) -> imag(λ(u,p))[i], us,ps), subplot=2, linewidth=2, color=:red )
end

plot!([NaN],[NaN], subplot=2, color=:gray, linestyle=:dot, linewidth=2, label=L"\Re\mathrm{e}" )
plot!([NaN],[NaN], subplot=2, color=:gray, linewidth=2, label=L"\Im\mathrm{m}" )

yticks!([0,1],subplot=1)
yticks!([-2,-1,0],subplot=2)
xticks!([0.8,1,1.2],subplot=2)

xticks!([NaN],subplot=1)

scatter!([1.045],[0],subplot=2,marker=:star,color=:black,markersize=7)
scatter!([0.92],[0],subplot=2,marker=:star,color=:black,markersize=7)

plot!([1.045,1.045],[0,0.94],subplot=1,linewidth=2,color=:gray)
plot!([0.92,0.92],[0,5],subplot=2,linewidth=2,color=:gray)

plot!([0.92,0.92],[0,0.94],subplot=1,linewidth=2,color=:gray)
plot!([1.045,1.045],[0,5],subplot=2,linewidth=2,color=:gray)

savefig("figures/hopf-measure.pdf")
