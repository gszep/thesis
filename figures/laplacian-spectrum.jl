using LinearAlgebra,StatsBase
using CairoMakie,Colors

⊗(A,B) = kron(A,B)
⊕(A,B) = A ⊗ I(size(B,1)) + I(size(A,1)) ⊗ B

function laplacian(N::Int)
    ∇² = zeros(Int8,N,N)

    # central differences
    ∇²[diagind(∇²)] .= -2 
    ∇²[diagind(∇²,1)] .= 1

    # periodic boundary conditions
    ∇²[1,end] = 1

    return Symmetric(∇²)
end

function eigenspectrum(A::AbstractMatrix; edges::AbstractVector=-5:0.05:2, mode=:pdf)
    return StatsBase.normalize( fit(Histogram,real(eigvals(A)),edges), mode=mode)
end

CairoMakie.activate!()
set_theme!(Theme(
    Axis = (
        xgridvisible  = false,
		ygridvisible  = false,
        xticksvisible  = false,
		yticksvisible  = false,
        # xticklabelsvisible  = false,
        yticklabelsvisible  = false,
    )
))

J = [1 1; -3 -2]
∇² = laplacian(2000)

stable = eigenspectrum( ∇² )
unstable = eigenspectrum( ∇² ⊗ Diagonal([1,100]) + I(size(∇²,1)) ⊗ J )
equaldiff = eigenspectrum( ∇² ⊗ Diagonal([1,1]) + I(size(∇²,1)) ⊗ J )

figure = Figure(resolution = 72 .* (7,7) )
ax = Axis(figure[1,1],
    ylabel = L"Eigenvalue Distribution $P(\lambda)$",
    xlabel = L"Eigenvalue $\Re$e$[\lambda]$"
)

scatter!(stable, color=:gray, markersize=3)
scatter!(equaldiff, color=:darkblue, markersize=3)
scatter!(unstable, color=:lightblue, markersize=3)

xlims!(-5,2)
figure

save("figures/laplacian-spectrum.pdf",figure)