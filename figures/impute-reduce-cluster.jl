using CSV,DataFrames
using CairoMakie,Colors
using StatsBase, GigaSOM

set_theme!(Theme(
	fontsize=48,
    Axis = (
        xgridvisible  = false,
		ygridvisible  = false,
        xticksvisible  = false,
		yticksvisible  = false,
        xticklabelsvisible  = false,
        yticklabelsvisible  = false,
		titlegap=25,
    ),
	Label = (
		tellheight=false,
		tellwidth=false
	)
))

pbmc = CSV.read("figures/pbmc.csv",DataFrame;header=false)|>Matrix

pbmc = pbmc[:,findall(StatsBase.var(pbmc,dims=1)[1,:].≠0)]
sample = pbmc[:,findall(StatsBase.var(pbmc,dims=1)[1,:].>1/2)]

xdim,ydim = 20,20
som = initGigaSOM(sample,xdim,ydim)
som = trainGigaSOM(som,sample)

clusters = mapToGigaSOM(som,sample)
embedding = embedGigaSOM(som,sample)

sample = pbmc[1:100,1:100]
sample[75:end,25:40] .= NaN
sample[25:75,70:80] .= NaN

figure = Figure(resolution = 72 .* (36,12) )

ax = Axis(figure[1,1],
	title = L"Impute$\quad$",
    ylabel = L"feature $u$",
    xlabel = L"sample $n$"
)

heatmap!(ax,sample,colorrange=(0,1), colormap=:heat)

ax = Axis(figure[1,2],
	title = L"Reduce$\quad$",
	ylabel = L"embedding coordinate $y$",
	xlabel = L"embedding coordinate $x$"
)

sc = scatter!(ax,embedding[:,1],embedding[:,2],markersize=7,color=log10.(pbmc[:,findmax(StatsBase.var(pbmc,dims=1)[1,:])[2]]))
Colorbar(figure[1,3],sc,label = L"feature $u$")

ax = Axis(figure[1,4],
	title = L"Cluster$\quad$",
	ylabel = L"embedding coordinate $y$",
	xlabel = L"embedding coordinate $x$"
)

scatter!(ax,embedding[:,1],embedding[:,2],markersize=7,color=:darkred)

figure