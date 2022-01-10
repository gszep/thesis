using DifferentialEquations,DSP
using ForwardDiff,LinearAlgebra
using CairoMakie,Colors,InvertedIndices

CairoMakie.activate!()

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
ax = Axis(figure[1,1],
    ylabel = L"Eigenvalue Distribution $P(\lambda)$",
    xlabel = L"Scaled Eigenvalue $\frac{\lambda}{L}$"
)

λ = range(-5,1,length=200)
function P(λ::Number)
    return real(1 / ( 2π * √Complex(1-(λ/2+1)^2) ))
end

p = P.(λ)
df = range(-5,1,length=length(p))[2] - range(-5,1,length=length(p))[1]
band!(range(-5,1,length=length(p)),p/sum(p*df),0,color=RGBA(0.5,0,0,0.2))

p = conv(P.(λ),P.(λ))
df = range(-5,1,length=length(p))[2] - range(-5,1,length=length(p))[1]
band!(range(-5,1,length=length(p)),p/sum(p*df),0,color=RGBA(0.5,0,0,0.2))

p = conv(P.(λ),p)
df = range(-5,1,length=length(p))[2] - range(-5,1,length=length(p))[1]
band!(range(-5,1,length=length(p)),p/sum(p*df),0,color=RGBA(0.5,0,0,0.2))

for i ∈ 1:10
    p = conv(P.(λ),p)
end

df = range(-5,1,length=length(p))[2] - range(-5,1,length=length(p))[1]
band!(range(-5,1,length=length(p)),p/sum(p*df)/2,0,color=RGBA(0.5,0,0,0.2))

figure
save("figures/laplacian-spectrum.pdf",figure)