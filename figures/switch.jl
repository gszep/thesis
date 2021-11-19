using DifferentialEquations, Catalyst
using Markdown,Latexify,Plots

LacReactions = @reaction_network begin

    # transcription
    α, LacI → LacI + LacR
    α, LacZ → LacZ + βGal

    # degredation
    μ, LacR → ∅
    μ, LacR⁻ → ∅
    μ, βGal → ∅

    # hydrolysis
    β*βGal, p → ∅

    # repressor binding
    (z,1), LacZ + LacR ↔ LacZ⁻
    (r,1), LacR + 2p   ↔ LacR⁻

end α β z r μ

rates,rates! = build_function(oderatelaw.(reactions(LacReactions)),
    states(LacReactions), parameters(LacReactions),
    independent_variable(LacReactions), expression=Val{false}
)

LacReactions

function F(u,p,θ)
    Z,R,G = u
    α,β,z,r,μ = θ

    return rates([1,R,Z,G,p*(β*G+p*r*R)/2,p,1-Z],[α,β,z,r,μ],0)
end

F(randn(3),randn(),randn(5))

LacSystem = convert(ODESystem,LacReactions)
calculate_tgrad(LacSystem)

p = (0.00166,0.0001,0.1)   # [c1,c2,c3]
tspan = (0., 100.)
u0 = [301., 100., 0., 0.]  # [S,E,SE,P]

# solve JumpProblem
dprob = DiscreteProblem(rs, species(rs) .=> u0, tspan, parameters(rs) .=> p)
jprob = JumpProblem(rs, dprob, Direct())
jsol = solve(jprob, SSAStepper())
plot(jsol,lw=2,title="Gillespie: Michaelis-Menten Enzyme Kinetics")