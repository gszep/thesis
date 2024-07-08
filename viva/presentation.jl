### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ e5db1a5a-9fb9-11ec-2ae1-4b4a6d89e282
begin 
	using JSServe
	Page()
end

# ╔═╡ cae037ca-af83-44ef-90a3-602990f6c63c
begin
	using StatsBase, LinearAlgebra, BifurcationInference, BifurcationKit
	using Parameters, Setfield, StaticArrays
	using LaTeXStrings, WGLMakie, GeometryBasics, TriplotBase, Triangulate, Colors, ColorSchemes
	
	using ForwardDiff, Flux, Logging
	using BifurcationInference: mean, window_function, ∂Fu, realeigvals, tangent_field, derivative, ∇measure, ∇errors
	
	global_logger(NullLogger())

	const style = JSServe.Dependency( :style, ["assets/style.css"])
	const d3 = JSServe.Dependency( :d3, ["assets/d3.v6.min.js"])
	const venn = JSServe.Dependency( :venn, ["assets/venn.js"])
	html"""<center><button onclick="
		present();

	    elements = document.getElementsByTagName('pluto-shoulder');
	    for (var i = 0; i < elements.length; i++) {

	        console.log(elements[i].style.display == '')
	        elements[i].style.display = elements[i].style.display == '' ? 'none' : '';
		};

	    elements = document.getElementsByTagName('pluto-trafficlight');
	    for (var i = 0; i < elements.length; i++) {

	        console.log(elements[i].style.display == '')
	        elements[i].style.display = elements[i].style.display == '' ? 'none' : '';
		};

	    elements = document.getElementsByClassName('add_cell');
	    for (var i = 0; i < elements.length; i++) {

	        console.log(elements[i].style.display == '')
	        elements[i].style.display = elements[i].style.display == '' ? 'none' : '';
		};
	
	">Present Mode</button></center>"""
end

# ╔═╡ c91f9c9a-75be-4f59-97f2-cf7ccbabf726
HTML("""<center>
<h1>Inferring Bifurcations <br/> Between Phenotypes </h1>
Grisha Szep
<br/><br/>
<img style="width:25%"
		src="$(JSServe.Asset("assets/kcl-crest.png")
)">
</center>""")

# ╔═╡ e3319933-a7d3-4239-92a2-23ca3e29cf63
HTML("""<h2>Thesis Contents</h2>""")

# ╔═╡ 53dc1a67-1362-42b4-8fa4-b366e370ead9
App() do session::Session
    Venn = DOM.div(id = "venn", class = "svg-container",
		style="width:100%; height:400px" )
    JSServe.onload(session, Venn, js"""function (container){

var svg = $d3.select("div#venn").append("svg")
	.classed("svg-content", true)
	.style('width','100%').style('height','100%')

function relativeCoords ( event ) {
  var bounds = event.target.getBoundingClientRect();
  var x = event.clientX - bounds.left;
  var y = event.clientY - bounds.top;
  return {x: x, y: y};
}

const sets = [
  { sets: ['Cell Biology'], size: 12, fields: 'Lab Automation</br>Multiomics</br>Pathways' },
  { sets: ['Dynamical Systems Theory'], size: 12, fields: 'Differential Equations </br>Bifurcations</br>Pattern Formation'  },
  { sets: ['Machine Learning'], size: 12, fields: 'Automatic Differentiation</br>Deep Learning' },
  { sets: ['Cell Biology', 'Dynamical Systems Theory'], size: 2, fields: 'Synthetic Biology</br>Systems Biology'  },
  { sets: ['Dynamical Systems Theory', 'Machine Learning'], size: 2, fields: 'Differential Programming</br>Geometric Learning</br>SciML'  },
  { sets: ['Machine Learning', 'Cell Biology'], size: 2, fields: 'Bioinformatics</br>Design Biology'  },
  { sets: ['Dynamical Systems Theory', 'Machine Learning', 'Cell Biology'], size: 2, fields: 'Design-Build-Test-Learn'  },
];

const chart = venn.VennDiagram();
svg.datum(sets).call(chart);
svg.selectAll('text').style('fill', 'white');
svg.selectAll('.venn-circle path').style('fill','#ddd');
svg.selectAll('.venn-intersection text')
	.style('fill','#ffd700').style('font-size',24);

var tooltip = d3.select("div#venn").append("div")
    .attr("class", "venntooltip").style("position", "absolute");

svg.selectAll("path")
    .style("stroke-opacity", 0)
    .style("stroke", "#fff")
    .style("stroke-width", 3)

svg.selectAll("g")
    .on("mouseover", function(event,d) {
        // sort all the areas relative to the current item
        venn.sortAreas(svg, d);
        
        // Display a tooltip with the current size
        tooltip.transition().duration(400).style("opacity", .9);
        tooltip.html(d.fields);

        // highlight the current path
        var selection = d3.select(this).transition("tooltip").duration(400);
        selection.select("path")
            .style("fill-opacity", d.sets.length == 1 ? .4 : .1)
            .style("stroke-opacity", 1)
            .style("stroke", d.sets.length == 3 ? '#ffd700' : "#fff")
            .style("fill", d.sets.length == 3 ? '#ffd700' : "#fff")
    })

    .on("mousemove", function(event) {
	    coords = relativeCoords(event)
        tooltip.style("left", (coords.x - 28) + "px")
               .style("top", (coords.y - 40) + "px");
    })

    .on("mouseout", function(event,d) {
        tooltip.transition().duration(400).style("opacity", 0);
        var selection = d3.select(this).transition("tooltip").duration(400);
        selection.select("path")
            .style("fill-opacity", d.sets.length == 1 ? .25 : .0)
            .style("stroke-opacity", 0);
    })
    .on("click", function(event,d) {
        document.querySelectorAll(".venn-intersection text tspan").forEach( function(x) {
			if (x.innerHTML == '') {

				switch(true) {
				  case (x.getAttribute('x') == 299 && x.getAttribute('y') == 195):
					x.innerHTML = '3'
					break;
				  case (x.getAttribute('x') == 342 && x.getAttribute('y') == 170):
					x.innerHTML = '4'
					break;
				  case (x.getAttribute('x') == 257 && x.getAttribute('y') == 170):
					x.innerHTML = '5'
					break;
				  default:
					x.innerHTML = ''
				}
			} else {
				x.innerHTML = ''
			}
    	})
    })

}""")
    return DOM.div(venn,style,Venn)
end

# ╔═╡ eedf0b49-cca7-4683-b205-05eb27104d69
HTML("""<h2>Interpretation of Morphogen Gradients</br>by a Synthetic Bistable Circuit</h2>""")

# ╔═╡ 0126a530-7322-4e75-be00-6e73420c5ec9
HTML("""<ul>
<li><b>How do dynamic morphogen gradients lead to robust gene expression domains?</b></li>
</br>
<center><img src=$(JSServe.Asset("assets/gap.jpg")) width=400px>
</br>
Figure 1: Expression patterns of pair-rule <i>Gap</i> genes in the <i>Drosophila</i> embryo.
</center>
</br>
</ul>
""")

# ╔═╡ 17220d8c-f64a-4707-a879-a2eed1875e9b
HTML("""<h2>Interpretation of Morphogen Gradients</br>by a Synthetic Bistable Circuit</h2>""")

# ╔═╡ 75a6e25c-2301-4fb7-a791-665d3d81363e
HTML("""<ul>
<li><b>How do dynamic morphogen gradients lead to robust gene expression domains?</b></li>
</br>
<center><img src=$(JSServe.Asset("assets/agar.svg")) width=400px>
</br>
Figure 2: Morphogen gradients <i>C<sub>6</sub></i> and <i>C<sub>12</sub></i> diffuse through Agar</br>acting on genetic switches in synthetic <i>E. coli</i> colony
</center>
</br>
</ul>
""")

# ╔═╡ 3f0e4940-9887-4827-9105-ef1c0082700c
HTML("""<h2>Interpretation of Morphogen Gradients</br>by a Synthetic Bistable Circuit</h2>""")

# ╔═╡ a0728fba-f63c-48cf-a2a3-43dd5a33c6d3
begin
	function switch(u,p)
		@unpack α, β, θ, origin = p
		
		F = similar(u)
		x,y = origin

		α′ = (α-x)*cos(θ) - (β-y)*sin(θ)
		β′ = (α-x)*sin(θ) + (β-y)*cos(θ)

		F[1] = α′ + β′*u[1] - u[1]^3
		return F
	end

	∂F(u,p) = ForwardDiff.jacobian(u->switch(u,p),u)
	N = 50; ∇² = SymTridiagonal([-1;fill(-2,N-2);-1],fill(1,N-1))
	
	function switch!(u,p)
		@unpack α, β, θ, origin, Dα, Dβ, kα, kβ, dt = p

		# reaction-diffusion
		u[] += @. ( (α[]-origin[][1])*cos(θ[]) - (β[]-origin[][2])*sin(θ[]) + ((α[]-origin[][1])*sin(θ[]) + (β[]-origin[][2])*cos(θ[]))*u[] - u[]^3 ) * dt[]
		α[] += ( Dα[]*∇²*α[] .- min.(u[],0)*kα[] ) * dt[]
		β[] += ( Dβ[]*∇²*β[] .+ max.(u[],0)*kβ[] ) * dt[]
	end

	options = ContinuationPar(dsmin = 0.01, dsmax = 0.1, ds = -0.02,
		pMin = -1.0, pMax = 1.0, maxSteps = 100,
		newtonOptions = NewtonPar(tol = 1e-6, maxIter=100), detectBifurcation=3
	)

	set_theme!(Theme(
		backgroundcolor = colorant"#1f1f1f",
		textcolor = :white,
	    Axis = (
	        xgridvisible  = false,
			ygridvisible  = false,
	        xticksvisible  = false,
			yticksvisible  = false,
	        xticklabelsvisible  = false,
	        yticklabelsvisible  = false,
	        ylabelpadding = 10,
			xlabelpadding = 10,
			titlesize = 24,
			xlabelsize = 24,
			ylabelsize = 24,
		    bottomspinecolor = :white,
			leftspinecolor = :white,
			rightspinecolor = :white,
			topspinecolor = :white,
		)
	))
	HTML("""<ul><li><b>
	How do dynamic morphogen gradients lead to robust gene expression domains?
	</b></li></ul>""")
end

# ╔═╡ cd1e5b29-cfb0-4de3-ac24-772ebb6c4b76
App() do session::Session
	u,x = Observable(zeros(N)), range(0,1,length=N)
	p = ( α = Observable(zeros(N)), β = Observable(zeros(N)),
		  Dα = Observable(0.1), Dβ = Observable(0.1),
	
	      kα = JSServe.Slider(range(0.0,0.01,length=50)),
		  kβ = JSServe.Slider(range(0.0,0.01,length=50)),

		  θ = Observable(π/4), origin = Observable([0.25,0.25]),
	      dt = JSServe.Slider(range(0.2,1.0,length=50)), u₀ = [1.0],
	)

	bistable_region = @lift begin
		branches = nothing
	
		branches, z = continuation( switch, ∂F, p.u₀,
			(α = 1.0, β = 1.0, θ = $(p.θ), origin = $(p.origin)), 
			(@lens _.α), options )

		if length(branches.specialpoint) ≠ 0
	
			branches, z = continuation( switch, ∂F, branches, 1,
				(@lens _.β), ContinuationPar(options, saveSolEveryStep=1))

			return map( s -> Point(s.x.p,s.p), branches.sol)

		else
			branches, z = continuation( switch, ∂F, p.u₀,
				(α = 1.0, β = 1.0, θ = $(p.θ), origin = $(p.origin)), 
				(@lens _.β), options )

			branches, z = continuation( switch, ∂F, branches, 1,
				(@lens _.α), ContinuationPar(options, saveSolEveryStep=1))

			return map( s -> Point(s.p,s.x.p), branches.sol)
		end
	end

	figure = Figure(resolution=(4.5*256,2*256))

	ax_space = figure[1:3,1] = Axis(figure, backgroundcolor=colorant"#2a2928",
		xlabel = L"Space $x$",
		ylabel = L"States $u(x,t)$")

	pressed = Observable("")
    on(events(ax_space.scene).mousebutton) do event
        if Makie.is_mouseinside(ax_space.scene) && event.action == Mouse.press
			pressed[] = event.button == Mouse.right ? "right" : "left"
        end
        if event.action ≠ Mouse.press && pressed[] ≠ ""
            pressed[] = ""
        end
    end

    on(events(ax_space.scene).mouseposition) do event
        if Makie.is_mouseinside(ax_space.scene) && pressed[] ≠ ""
	
            xp,yp = mouseposition(ax_space.scene)
			i = findfirst(x.>xp)

			if !isnothing(i)

				pressed[] == "left" ? p.α[][i] = yp : p.β[][i] = yp
				p.α[] = p.α[]; p.β[] = p.β[]
			end
        end
    end
	
	deactivate_interaction!(ax_space,:rectanglezoom)
	deactivate_interaction!(ax_space,:scrollzoom)
	deactivate_interaction!(ax_space,:dragpan)
	
	xlims!(ax_space,0,1)
	ylims!(ax_space,0,3)

	cfp = @lift(-min.($u,0))
	yfp = @lift(max.($u,0))
	
	c = band!(ax_space, x, zeros(N), cfp , linewidth=5, color=colorant"#00b0f055")
	y = band!(ax_space, x, yfp, zeros(N), linewidth=5, color=colorant"#ffd70055")

	b = lines!(ax_space, x, p.β, linewidth=5, color=colorant"#000099")
	a = lines!(ax_space, x, p.α, linewidth=5, color=colorant"#ff6600")

	Legend(figure[1,2], [c,y], [L"CFP",L"YFP"], titlesize=24,
		L"Fluorescence $u_1,u_2$", framevisible = false, labelsize=18)
	Legend(figure[2,2], [b,a], [L"C_6",L"C_{12}"], titlesize=24,
		L"Morphogens $u_3,u_4$", framevisible = false, labelsize=20)

	ax_states = figure[1:3,3] = Axis(figure, backgroundcolor=colorant"#2a2928",
		xlabel = L"Signal $C_{12}(x,t)$",
		ylabel = L"Signal $C_6(x,t)$")

	edge_length = 10
    on(events(ax_states.scene).mousebutton) do event
        if event.button == Mouse.right && Makie.is_mouseinside(ax_states.scene)
			
            xp,yp = mouseposition(ax_states.scene)
            if event.action == Mouse.press
				
					p.α[] .= 0
					p.α[][1:edge_length] .= xp*N/edge_length
					
					p.β[] .= 0
					p.β[][end-edge_length+1:end] .= yp*N/edge_length
			
					p.α[] = p.α[]
					p.β[] = p.β[]
				
                return Consume(false)
            end
        end
        return Consume(false)
    end

	circle = select_line(ax_states.scene)
	on(circle) do (r,dr)

		p.θ[] = atan((dr-r)...)
		p.origin[] = r
	end
	
	deactivate_interaction!(ax_states,:rectanglezoom)
	deactivate_interaction!(ax_states,:scrollzoom)
	
	xlims!(ax_states,0.0,1.0)
	ylims!(ax_states,0.0,1.0)
	
	region = lines!(ax_states, bistable_region, linewidth=5, color=:white)
	colors = @lift(map( ui -> ui > 0 ? colorant"#ffc000" : colorant"#00b0f0", $u))
	
	arrows!(ax_states, p.α, p.β, @lift( 10* $cfp * $(p.kα) ), @lift( 10* $yfp * $(p.kβ)),
		color=colors, linewidth=3 )
	scatter!(ax_states, p.α, p.β, color=colors)
	
	avg = scatter!(ax_states, @lift([StatsBase.mean($(p.α))]), @lift([StatsBase.mean($(p.β))]),
		marker=:x, color=:white, markersize=15)

	Legend(figure[3,2], [avg,region], ["Spatial Average","Limit Points"],
		framevisible = false)

	play = JSServe.Button("Play")
	playing = Observable(false)
	
    on(play) do click; playing[] = !playing[]; end
	on(play) do click
	    @async while playing[]
	
	        switch!(u,p)
	        sleep(0.01)
	    end
	end
	return DOM.div( style="text-align: center", figure, play, html"&nbsp;&nbsp;&nbsp;&nbsp;", p.dt, p.kα, p.kβ )
end

# ╔═╡ 8cddf5fb-7b93-43c6-8c61-b9931d69fd83
HTML("""<h2>Interpretation of Morphogen Gradients</br>by a Synthetic Bistable Circuit</h2>""")

# ╔═╡ 0fd1d551-49b0-4194-b668-fec7aa3ed9a7
HTML("""<ul>
<li><b>Which genetic designs satisfy a target cusp bifurcation in flow cytometry data?</b></li>
</br>
<center><img src=$(JSServe.Asset("assets/hysteresis.svg")) width=600px>
</br>
Figure 3: Extracting limit points with respect to conditions <i>p,p\'</i>. Steady state distributions and population separation
</center>
</br>
</ul>
""")

# ╔═╡ e2c520ec-2717-42c4-b4f0-0da56cf59927
HTML("""<h2>Parameter Inference</br>with Bifurcation Diagrams</h2>""")

# ╔═╡ 03619862-196e-419d-9c37-2e9ada04c84d
HTML("""<ul><li><b>Which differential equations satisfy a target bifurcation diagram?</b></li><li><b>How do we organise models in terms of geometric and topological equivalence?</b></li></ul>""")

# ╔═╡ 55e55af6-f0c2-4fbd-96a5-dcc1bd2f1dd4
L"""
$\frac{\partial u}{\partial t} = F_\theta(u,p) \quad \mathrm{where} \quad F_\theta(u,p) := p + \theta_1 u + \theta_2 u^3$
"""

# ╔═╡ 9d69479a-e21b-4316-b46d-c36bbf634f0b
begin
	F(z::BorderedArray,θ::AbstractVector) = F(z.u,(θ=θ,p=z.p))
	function F(u::AbstractVector,parameters::NamedTuple)

		@unpack θ,p = parameters
	
		f = first(u)*first(p)*first(θ)
		F = similar(u,typeof(f))
	
		F[1] = p + θ[1]*u[1] + θ[2]*u[1]^3
		return F
	end
end

# ╔═╡ 2c67684f-392d-4a48-abba-be0fa8079344
HTML("""<h2>Parameter Inference</br>with Bifurcation Diagrams</h2>""")

# ╔═╡ 10ee39be-371b-487c-ac7c-9cbe60d2b367
HTML("""<ul><li><b>Which differential equations satisfy a target bifurcation diagram?</b></li><li><b>How do we organise models in terms of geometric and topological equivalence?</b></li></ul>""")

# ╔═╡ a127bc38-a03c-494a-a8c0-142d1674dc6c
L"""
$L(\theta) := \Big\langle |p(\theta)-p'| \Big\rangle_{\! p,\,p'} - \mu\log\Psi(\theta)
"""

# ╔═╡ 65bc7068-0cb1-4d54-9b81-8ec307893db8
L"""
\Psi(\theta):=\frac{\int_{F_\theta=0}\varphi_\theta(s)\,\mathrm{d}s}{\int_{F_\theta=0}\mathrm{d}s}
\quad\mathrm{where}\quad
\varphi_\theta(s) := \sum_{\lambda\in\frac{\partial F_\theta}{\partial u}}\left(\left|\frac{d}{ds}\log\Re\mathrm{e}[\lambda(s)]\right|^{-1}+1\right)^{-1}
"""

# ╔═╡ d88dddba-d9e1-4aaf-86ca-b107bb937fd8
begin 
	import BifurcationInference: measure
	function measure(F::Function,z::BorderedArray,θ::AbstractVector,t::StateSpace)
		
		λ = realeigvals(∂Fu(F,z,θ))
		dλ = derivative( z -> realeigvals(∂Fu(F,z,θ)), z, tangent_field(F,z,θ) )
		
		return window_function(t.parameter,z) * mapreduce(
			(λ,dλ) -> 1 / ( 1 + abs(λ/dλ)^2 ), +, λ, dλ
		)
	end
end

# ╔═╡ cc0c9d38-eb52-4f9e-87bf-70aeef1847db
App() do session::Session

	X = StateSpace( 1, -2:0.01:2, [-1.0,1.0] )
	parameters = ( θ=SizedVector{2}(5.0,-0.93), p=minimum(X.parameter) )
	hyperparameters = getParameters(X)

	r,α = range(0.01,5,length=10), range(0.01-π,π-0.01,length=40)
	θ₁ = vec(@. r*cos(α')); θ₂ = vec(@. r*sin(α'))

	triangles = first(triangulate("Q",
		Triangulate.TriangulateIO( pointlist = [θ₁'; θ₂'] ))
	).trianglelist

	state_landscape = map( (x,y) -> deflationContinuation( F, X.roots,
		( θ=SizedVector{2}(x,y), p=parameters.p ), hyperparameters ), θ₁, θ₂)

	predictions = map( branches -> unique([ s.z for branch ∈ branches for s ∈ branch if s.bif ]; atol=2*step(X.parameter)), state_landscape)

	function norm( F::Function, z::BorderedArray, θ::AbstractVector, targets::StateSpace )
		return window_function(targets.parameter,z)
	end

	function norm( F::Function, branch::Vector{<:NamedTuple}, θ::AbstractVector, targets::StateSpace )
		return sum( s -> norm(F,s.z,θ,targets)*s.ds, branch )
	end
	
	function norm( F::Function, branches::AbstractVector{<:AbstractVector}, θ::AbstractVector, targets::StateSpace )
		return sum( branch -> norm(F,branch,θ,targets), branches )
	end

	measures = map( (x,y,s,p) -> (length(p)==0)*(log(measure(F,s,SizedVector{2}(x,y),X))-log(norm(F,s,SizedVector{2}(x,y),X))), θ₁, θ₂, state_landscape, predictions)
	
	optimiser = Flux.Optimise.Momentum(0.01)
	targets, θ = Observable([X.targets...]), Observable(parameters.θ)

	levels = [0.1,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
	contours = @lift begin
		errors = map( x -> mean( p′-> mean( z->(z.p-p′)^2, x; type=:geometric ), $targets; type=:arithmetic ), predictions )
		return TriplotBase.tricontour(θ₁,θ₂,asinh.(errors-measures),triangles,levels)
	end

	steady_states = @lift begin
		return deflationContinuation( F, X.roots, ( θ=$θ, p=parameters.p ), hyperparameters )
	end
	
	∂L = @lift begin

		p = unique([ s.z for branch ∈ $steady_states for s ∈ branch if s.bif ]; atol=2*step(X.parameter))
		t = StateSpace{1,Float64}(X.roots,X.parameter,$targets)
	
		if length(p) == 0 
			return -∇measure(F,$steady_states,$θ,t)/measure(F,$steady_states,$θ,t)
		else
			return ∇errors(F,p,$θ,t)
		end
	end

	bifurcation_diagram = @lift begin
		return vcat(map( branch -> [map( s -> Point(s.z.p,s.z.u[1]), branch )..., Point(NaN,NaN)], $steady_states )...)
	end

	bifurcation_colors = @lift begin
		return vcat(map( branch -> [map( s -> real(s.λ[1]) < 0 ? :blue : :lightblue, branch )..., :transparent], $steady_states )...)
	end

	smin,smax = extrema(levels)
	colormap = cgrad([:white,colorant"#2a2928",:black], [0,0.5,1])
	landscape = @lift begin
		
		contour_lines = Point{2,Float32}[]
		contour_colors = RGB{Float32}[]
		
		for contour ∈ $contours
			s = (contour.level - smin) / (smax - smin)

			for x ∈ contour.polylines
				
				points = [map(Point,x)..., Point(NaN,NaN)]
				append!(contour_lines,points)
				
				colors = get(colormap,fill(s,size(points)))
				append!(contour_colors,colors)
			end
		end

		return contour_lines,contour_colors
	end

	figure = Figure(resolution=(4.5*256,2*256))
	ax_cost = figure[1:2,1] = Axis(figure, backgroundcolor=colorant"#2a2928",

		title = L"Cost Landscape $L(\theta)$",
		xlabel = L"parameter $\theta_1$",
		ylabel = L"parameter $\theta_2$")

	lines!(ax_cost, lift(first,landscape), color=lift(last,landscape), linewidth=5 )
	scatter!(ax_cost, @lift([ $θ ]), color=get(colormap,0))
	Colorbar(figure[1:2,2], limits = (smin,smax), colormap = colormap)

    on(events(ax_cost.scene).mousebutton) do event
        if event.button == Mouse.left && Makie.is_mouseinside(ax_cost.scene)
            if event.action == Mouse.press
				θ[] = mouseposition(ax_cost.scene)
            end
        end
    end

	deactivate_interaction!(ax_cost,:rectanglezoom)
	deactivate_interaction!(ax_cost,:scrollzoom)
	deactivate_interaction!(ax_cost,:dragpan)
	
	xlims!(ax_cost,-5,5)
	ylims!(ax_cost,-5,5)

	ax_diagram = figure[1:2,3] = Axis(figure, backgroundcolor=colorant"#2a2928",
		title = L"Bifurcation Diagram $F_\theta(u,p)=0$",
		xlabel = L"condition $p$",
		ylabel = L"fixed points $u$")

	bistable_select = select_line(ax_diagram.scene, linewidth=1, color=:gold)
	on(bistable_select) do (r,dr)
		targets[] = [r[1],dr[1]]
		@show targets[]
	end

	deactivate_interaction!(ax_diagram,:rectanglezoom)
	deactivate_interaction!(ax_diagram,:scrollzoom)
	deactivate_interaction!(ax_diagram,:dragpan)

	lines!(ax_diagram, bifurcation_diagram, color=bifurcation_colors, linewidth=5)
	t = vlines!(ax_diagram, targets, linewidth=1, color=:gold )

	xlims!(ax_diagram,-2,2)
	ylims!(ax_diagram,-4,4)

	s = lines!(ax_diagram, [Point(NaN,NaN)], color=:blue, linewidth=5)
	u = lines!(ax_diagram, [Point(NaN,NaN)], color=:lightblue, linewidth=5)

	Legend(figure[1:2,4], [s,u,t],
		[L"stable $\Re$e$[\lambda]<0$",L"unstable $\Re$e$[\lambda]>0$",L"targets $p$"],
		framevisible = false, labelsize=20)
	
	play = JSServe.Button("Play")
	playing = Observable(false)
	
    on(play) do click; playing[] = !playing[]; end
	on(play) do click
	    @async while playing[]

	        Flux.update!(optimiser, θ[], ∂L[] )
			θ[] = θ[]
	        sleep(0.01)
	    end
	end
	return DOM.div( style="text-align: center", figure, play )
end

# ╔═╡ 492b8dbf-3aac-4966-ba4c-03c78fea787d
HTML("""<h2>Exploring Bifurcations</br>Between Phenotypes</h2>""")

# ╔═╡ c3501463-f43a-465c-a81d-4bd5629793a3
HTML("""<ul>
<li><b>How do we identify qualitatively distinct cell populations in flow cytometry data?</b></li>
</br>
<center><img src=$(JSServe.Asset("assets/impute-reduce-cluster.svg")) width=600px>
</br>
Figure 4 : Impute, Reduce, Cluster pipeline
</center>
</br>
</ul>
""")

# ╔═╡ 47ecff8a-c30d-40b5-87a0-324e39493f1e
HTML("""<h2>Conclusions</h2>""")

# ╔═╡ f7a1eb21-e615-4209-baa2-d682aefa1f1b
HTML("""
<details><summary>
<b style='color:#ffd700'>6.1.1</b> A design–learn workflow for synthetic biology
</summary>
</br>
<center><img src=$(JSServe.Asset("assets/design-learn.svg")) width=600px>
</br>
Figure 5 : Overview for a <i>design-learn</i> workflow, developed in hindsight, for genetic design of the <i>double exclusive reporter</i>
</center>
</br>
</details>
<details><summary>
<b style='color:#ffd700'>6.1.2</b> Bifurcations and model order reduction
</summary>
</br>
<center><img src=$(JSServe.Asset("assets/model-reduction.svg")) width=600px>
</br>
Figure 6 : Overview of how <i>FlowAtlas.jl</i> can be adapted to explore the space of models with respect to a given </br> control condition. The parameter embedding would align clusters of equivalent bifurcations in different models
</center>
</br></details>""")

# ╔═╡ fc60be23-b455-4803-930e-1738f70a3c76
HTML("""<h2>Future Work</h2>""")

# ╔═╡ d1eb323f-f72f-4538-b4de-1c9541acf453
App() do session::Session
	function fixed_point(r::AbstractVector, r₀::AbstractVector; stable::Bool=true, α::Real=1)
	    return -ForwardDiff.gradient( r -> exp(-α*LinearAlgebra.norm(r-r₀)) * ( stable ? -1 : 1 ), r )
	end
	
	function limit_cycle(r::Point2{T}, f::Function; stable::Bool=true, ω::Real=1) where T<:Real
	    ∇F, F = ForwardDiff.gradient(r->abs(f(r)),r), abs(f(r))
	
	    normal = ∇F / LinearAlgebra.norm(∇F)
	    tangent = Point2{T}(-∇F[2],∇F[1]) / LinearAlgebra.norm(∇F)
	
	    return - ( stable ? 1 : -1 ) * ( tangent*ω*sign(f(r))*exp(-F) + normal*F*exp(-F) )
	end
	
	figure = Figure(resolution = 72 .* (5,5) )
	x,y = -5..2,-2..4
	minima, maxima = minimum.([x,y]), maximum.([x,y])
	ax_field = Axis(figure[1,1], backgroundcolor=colorant"#2a2928",
	
	    ylabel = "",
	    xlabel = L"vector field $F_\theta$",
	
	    xlims = (minima[1],maxima[1]),
	    ylims = (minima[2],maxima[2])
	)
	
	play = JSServe.Button("Play")
	point, circle = Observable(Point(1.0,-1.0)), select_line(ax_field.scene, color=:gold, linewidth=5)
	circle[] = [[-2,2],[-2,1]]

    on(events(ax_field.scene).mousebutton) do event
        if event.button == Mouse.right && Makie.is_mouseinside(ax_field.scene)
			
            xp,yp = mouseposition(ax_field.scene)
            if event.action == Mouse.press
				point[] = [xp,yp]
            end
        end
    end
	deactivate_interaction!(ax_field,:rectanglezoom)
	deactivate_interaction!(ax_field,:scrollzoom)
	deactivate_interaction!(ax_field,:dragpan)
	
	field = @lift( x-> limit_cycle(x,r->LinearAlgebra.norm(r-first($circle))-LinearAlgebra.norm(first($circle)-last($circle))) + fixed_point(x,$point) )
	streamplot!( ax_field, field, x, y, density=1, linewidth=5, colormap=[:white], arrow_size = 0)
	
	scatter!( ax_field, point, markersize=16, color=:blue)
	scatter!( ax_field, @lift(first($circle)), markersize=16, color=:lightblue)
	
	return DOM.div(HTML("""<details><summary> <b style='color:#ffd700'>6.3.1</b> Designing Limit Cycles</summary>
		</br> <center><img src=$(JSServe.Asset("assets/design-limit-cycles.svg")) width=300px>"""), figure,
		HTML("""Figure 7 : Using global constraints with basis functions in state space</br> and local constraints can help design global bifurcations of limit cycles</center></br>""") )
end

# ╔═╡ 440651dd-2ee5-46c0-a9ea-a5bce1093643
HTML("""<h2>Future Work</h2>""")

# ╔═╡ 60d407f7-0b36-4355-b4e0-ef0c6d83f800
HTML("""
<details><summary>
<b style='color:#ffd700'>6.3.2</b> Spatially Extended Systems
</summary>
</details>""")

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BifurcationInference = "7fe238d6-d31e-4646-aa16-9d8429fd6da8"
BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
Flux = "587475ba-b771-5e3f-ad9e-33799f191a9c"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
JSServe = "824d6782-a2ef-11e9-3a09-e5662e0c26f9"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Logging = "56ddb016-857b-54e1-b83d-db4d58db5568"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
Setfield = "efcf1570-3423-57d1-acb7-fd33fddbac46"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Triangulate = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
TriplotBase = "981d1d27-644d-49a2-9326-4793e63143c3"
WGLMakie = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"

[compat]
BifurcationInference = "~0.1.3"
BifurcationKit = "~0.1.11"
ColorSchemes = "~3.17.1"
Colors = "~0.12.8"
Flux = "~0.12.9"
ForwardDiff = "~0.10.25"
GeometryBasics = "~0.3.10"
JSServe = "~1.2.5"
LaTeXStrings = "~1.3.0"
Parameters = "~0.12.3"
Setfield = "~0.8.2"
StaticArrays = "~1.4.2"
StatsBase = "~0.33.16"
Triangulate = "~2.1.2"
TriplotBase = "~0.1.0"
WGLMakie = "~0.4.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cde29ddf7e5726c9fb511f340244ea3481267608"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.2"

[[AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[Arpack_jll]]
deps = ["Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "e214a9b9bd1b4e1b4f15b22c0994862b66af7ff7"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.0+3"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "81f0cb60dc994ca17f68d9fb7c942a5ae70d9ee4"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "5.0.8"

[[ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[ArrayInterfaceStaticArraysCore]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "01a9f8e6cfc2bfdd01d333f70b8014a04893103c"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.4"

[[ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "734d86dfee527f71f61370f410fb3f6c26740d33"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.19"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c06a868224ecba914baa6942988e2f2aade419be"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "0.1.0"

[[Automa]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "ef9997b3d5547c48b41c7bd8899e812a917b409d"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.4"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[BFloat16s]]
deps = ["LinearAlgebra", "Printf", "Random", "Test"]
git-tree-sha1 = "a598ecb0d717092b5539dbbe890c98bac842b072"
uuid = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"
version = "0.2.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BifurcationInference]]
deps = ["BifurcationKit", "Flux", "ForwardDiff", "InvertedIndices", "LaTeXStrings", "LinearAlgebra", "Parameters", "Plots", "Setfield", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "624bf942aab1c8eb63efa235b02b97d9b8e0b356"
uuid = "7fe238d6-d31e-4646-aa16-9d8429fd6da8"
version = "0.1.3"

[[BifurcationKit]]
deps = ["ArnoldiMethod", "Arpack", "BlockArrays", "DataStructures", "Dates", "DocStringExtensions", "FastGaussQuadrature", "ForwardDiff", "IterativeSolvers", "KrylovKit", "LinearAlgebra", "LinearMaps", "Parameters", "Printf", "RecipesBase", "RecursiveArrayTools", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StructArrays"]
git-tree-sha1 = "d9e4f3d358d966013403317363dfbfc2f13d18b6"
uuid = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
version = "0.1.12"

[[BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "3b15c61bcece7c426ea641d143c808ace3661973"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "0.16.25"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[CUDA]]
deps = ["AbstractFFTs", "Adapt", "BFloat16s", "CEnum", "CompilerSupportLibraries_jll", "ExprTools", "GPUArrays", "GPUCompiler", "LLVM", "LazyArtifacts", "Libdl", "LinearAlgebra", "Logging", "Printf", "Random", "Random123", "RandomNumbers", "Reexport", "Requires", "SparseArrays", "SpecialFunctions", "TimerOutputs"]
git-tree-sha1 = "6717cb9a3425ebb7b31ca4f832823615d175f64a"
uuid = "052768ef-5323-5732-b1bb-66c8b64840ba"
version = "3.13.1"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[ChainRules]]
deps = ["ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "Statistics"]
git-tree-sha1 = "f9d6dd293ed05233d37a5644f880f5def9fdfae3"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.42.0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"

[[ChangesOfVariables]]
deps = ["InverseFunctions", "LinearAlgebra", "Test"]
git-tree-sha1 = "2fba81a302a7be671aefe194f0525ef231104e7f"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.8"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "b8fe8546d52ca154ac556809e10c75e6e7430ac8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.5"

[[ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "d476eaeddfcdf0de15a67a948331c69a585495fa"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.47.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["AliasTables", "ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9c405847cc7ecda2dc921ccf18b47ca150d7317e"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.109"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "03b753748fd193a7f2730c02d880da27c5a24508"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.6.0"

[[EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[EnzymeCore]]
deps = ["Adapt"]
git-tree-sha1 = "3a3177ba05b4763234819060fb6c2e1613379ca6"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.7.6"

[[EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "58d83dd5a78a36205bdfddb82b1bb67682e64487"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "0.4.9"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "7072f1e3e5a8be51d525d64f63d3ec1287ff2790"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.11"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[Flux]]
deps = ["AbstractTrees", "Adapt", "ArrayInterface", "CUDA", "CodecZlib", "Colors", "DelimitedFiles", "Functors", "Juno", "LinearAlgebra", "MacroTools", "NNlib", "NNlibCUDA", "Pkg", "Printf", "Random", "Reexport", "SHA", "SparseArrays", "Statistics", "StatsBase", "Test", "ZipFile", "Zygote"]
git-tree-sha1 = "511b7c48eebb602a8f63e7d6c63e25633468dc16"
uuid = "587475ba-b771-5e3f-ad9e-33799f191a9c"
version = "0.12.10"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[Formatting]]
deps = ["Logging", "Printf"]
git-tree-sha1 = "fb409abab2caf118986fc597ba84b50cbaf00b87"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.3"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"

[[FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "StaticArrays"]
git-tree-sha1 = "d51e69f0a2f8a3842bca4183b700cf3d9acce626"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.1"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[Functors]]
git-tree-sha1 = "223fffa49ca0ff9ce4f875be001ffe173b2b7de4"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.2.8"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "3f74912a156096bd8fdbef211eff66ab446e7297"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+0"

[[GPUArrays]]
deps = ["Adapt", "GPUArraysCore", "LLVM", "LinearAlgebra", "Printf", "Random", "Reexport", "Serialization", "Statistics"]
git-tree-sha1 = "2e57b4a4f9cc15e85a24d603256fe08e527f48d1"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "8.8.1"

[[GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[GPUCompiler]]
deps = ["ExprTools", "InteractiveUtils", "LLVM", "Libdl", "Logging", "TimerOutputs", "UUIDs"]
git-tree-sha1 = "19d693666a304e8c371798f4900f7435558c7cde"
uuid = "61eb1bfa-7361-4325-ad38-22787b887f55"
version = "0.17.3"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "cf0a9940f250dc3cb6cc6c6821b4bf8a4286cf9c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.66.2"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "bc9f7725571ddb4ab2c4bc74fa397c1c5ad08943"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.69.1+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "82853ebc70db4f5a3084853738c68fd497b22c7c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.3.10"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43ba3d3c82c18d88471cfd2924931658838c9d8f"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+4"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Match", "Observables"]
git-tree-sha1 = "d44945bdc7a462fa68bb847759294669352bd0a4"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.5.7"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[IRTools]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "950c3717af761bc3ff906c2e8e52bd83390b6ec2"
uuid = "7869d1d1-7146-5819-86e3-90919afe41df"
version = "0.4.14"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[ImageIO]]
deps = ["FileIO", "Netpbm", "PNGFiles"]
git-tree-sha1 = "0d6d09c28d67611c68e25af0c2df7269c82b73c7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.4.1"

[[ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils"]
git-tree-sha1 = "8e2eae13d144d545ef829324f1f0a5a4fe4340f3"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.3.1"

[[ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "8d2e786fd090199a91ecbf4a66d03aedd0fb24d4"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.11+4"

[[ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14eb2b542e748570b56446f4c50fbfb2306ebc45"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.0+0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "bcf640979ee55b652f3b01650444eb7bbe3ea837"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.4"

[[InverseFunctions]]
deps = ["Dates", "Test"]
git-tree-sha1 = "e7cbed5032c4c397a6ac23d1493f3289e01231c4"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.14"

[[InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "59545b0a2b27208b0650df0a46b8e3019f85055b"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.4"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "eb3edce0ed4fa32f75a0a11217433c31d56bd48b"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.0"

[[JSServe]]
deps = ["Base64", "CodecZlib", "Colors", "HTTP", "Hyperscript", "JSON3", "LinearAlgebra", "Markdown", "MsgPack", "Observables", "RelocatableFolders", "SHA", "Sockets", "Tables", "Test", "UUIDs", "WebSockets", "WidgetsBase"]
git-tree-sha1 = "8f90ecadd55a7020da1929d4b0ea7a53aad2f404"
uuid = "824d6782-a2ef-11e9-3a09-e5662e0c26f9"
version = "1.2.6"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[KernelAbstractions]]
deps = ["Adapt", "Atomix", "EnzymeCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "PrecompileTools", "Requires", "SparseArrays", "StaticArrays", "UUIDs", "UnsafeAtomics", "UnsafeAtomicsLLVM"]
git-tree-sha1 = "d0448cebd5919e06ca5edc7a264631790de810ec"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.22"

[[KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "49b0c1dd5c292870577b8f58c51072bd558febb9"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.4"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "f044a2796a9e18e0531b9b3072b0019a61f264bc"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.17.1"

[[LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "070e4b5b65827f82c16ae0916376cb47377aa1b5"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.18+0"

[[LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "8c57307b5d9bb3be1ff2da469063628631d4d51e"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.21"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LinearMaps]]
deps = ["ChainRulesCore", "LinearAlgebra", "SparseArrays", "Statistics"]
git-tree-sha1 = "ee79c3208e55786de58f8dcccca098ced79f743f"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.11.3"

[[LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Showoff", "SignedDistanceFields", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "UnicodeFun"]
git-tree-sha1 = "d03c5a4056707bb8d43e349bc2cb49fc1cfa8b9f"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.15.1"

[[MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "7bcc8323fb37523a6a51ade2234eee27a11114c8"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.1.3"

[[MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "Test"]
git-tree-sha1 = "69b565c0ca7bf9dae18498b52431f854147ecbf3"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.1.2"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "f5db02ae992c260e4826fe78c942954b48e1d9c2"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.1"

[[NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Pkg", "Random", "Requires", "Statistics"]
git-tree-sha1 = "72240e3f5ca031937bd536182cb2c031da5f46dd"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.8.21"

[[NNlibCUDA]]
deps = ["Adapt", "CUDA", "LinearAlgebra", "NNlib", "Random", "Statistics"]
git-tree-sha1 = "b05a082b08a3af0e5c576883bc6dfb6513e7e478"
uuid = "a00861dc-f156-4864-bf3c-e6376f28a68d"
version = "0.2.6"

[[NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

[[Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "f4049d379326c2c7aa875c702ad19346ecb2b004"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.4.1"

[[PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "3f9b0706d6051d8edf9959e2422666703080722a"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.32.0"

[[PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[PtrArrays]]
git-tree-sha1 = "f011fbb92c4d401059b2212c05c0601b70f8b759"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.0"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "c860e84651f58ce240dd79e5d9e055d55234c35a"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.2"

[[RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"

[[RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "Tables", "ZygoteRules"]
git-tree-sha1 = "a5ce741acddc02f0d4fc6505463ca89697d7fb23"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.32.3"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d483cd324ce5cf5d61b77930f0bbd6cb61927d21"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.2+0"

[[RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "StaticArraysCore", "Statistics", "Tables"]
git-tree-sha1 = "fe89a8113ea445bcff9ee570077830674babb534"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.81.0"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "63c6b8796d28a1f942c29659e5519e2ef9ef4a59"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.2.7"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

[[StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "5d2c08cef80c7a3a8ba9ca023031a85c263012c5"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.6"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2bbd9f2e40afd197a1379aef05e0d85dba649951"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.7"

[[StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5950925ff997ed6fb3e985dcce8eb1ba42a0bbe7"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.18"

[[StructArrays]]
deps = ["DataAPI", "Tables"]
git-tree-sha1 = "ad1f5fd155426dcc879ec6ede9f74eb3a2d582df"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.4.2"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "ca4bccb03acf9faaf4137a9abc1881ed1841aa70"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.10.0"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "5a13ae8a41237cff5ecf34f73eb1b8f42fff6531"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.24"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[Triangle_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "fe28e9a4684f6f54e868b9136afb8fd11f1734a7"
uuid = "5639c1d2-226c-5e70-8d55-b3095415a16a"
version = "1.6.2+0"

[[Triangulate]]
deps = ["DocStringExtensions", "Libdl", "Printf", "Test", "Triangle_jll"]
git-tree-sha1 = "796a9c0b02a3414af6065098bb7cf0e88dfa450e"
uuid = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
version = "2.1.3"

[[TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

[[UnsafeAtomicsLLVM]]
deps = ["LLVM", "UnsafeAtomics"]
git-tree-sha1 = "ead6292c02aab389cb29fe64cc9375765ab1e219"
uuid = "d80eeb9a-aca5-4d75-85e5-170c8b632249"
version = "0.1.1"

[[Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[WGLMakie]]
deps = ["Colors", "FileIO", "FreeTypeAbstraction", "GeometryBasics", "Hyperscript", "ImageMagick", "JSServe", "LinearAlgebra", "Makie", "Observables", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "56d9ca8f1e1c2c09c1ccdf19c4ce5d45efd66097"
uuid = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"
version = "0.4.5"

[[Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "f91a602e25fe6b89afc93cf02a4ae18ee9384ce3"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.5.9"

[[WidgetsBase]]
deps = ["Observables"]
git-tree-sha1 = "30a1d631eb06e8c868c559599f915a62d55c2601"
uuid = "eead4739-05f7-45a1-878c-cee36b57321c"
version = "0.1.4"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "d9717ce3518dc68a99e6b96300813760d887a01d"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.1+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "3593e69e469d2111389a9bd06bac1f3d730ac6de"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.4"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[Zygote]]
deps = ["AbstractFFTs", "ChainRules", "ChainRulesCore", "DiffRules", "Distributed", "FillArrays", "ForwardDiff", "GPUArrays", "GPUArraysCore", "IRTools", "InteractiveUtils", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NaNMath", "Random", "Requires", "SparseArrays", "SpecialFunctions", "Statistics", "ZygoteRules"]
git-tree-sha1 = "8ac61a92a33b3fd2a4cbf92951817831e313a004"
uuid = "e88e6eb3-aa80-5325-afca-941959d7151f"
version = "0.6.44"

[[ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "27798139afc0a2afa7b1824c206d5e87ea587a00"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.5"

[[isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─e5db1a5a-9fb9-11ec-2ae1-4b4a6d89e282
# ╟─cae037ca-af83-44ef-90a3-602990f6c63c
# ╟─c91f9c9a-75be-4f59-97f2-cf7ccbabf726
# ╟─e3319933-a7d3-4239-92a2-23ca3e29cf63
# ╟─53dc1a67-1362-42b4-8fa4-b366e370ead9
# ╟─eedf0b49-cca7-4683-b205-05eb27104d69
# ╟─0126a530-7322-4e75-be00-6e73420c5ec9
# ╟─17220d8c-f64a-4707-a879-a2eed1875e9b
# ╟─75a6e25c-2301-4fb7-a791-665d3d81363e
# ╟─3f0e4940-9887-4827-9105-ef1c0082700c
# ╟─a0728fba-f63c-48cf-a2a3-43dd5a33c6d3
# ╟─cd1e5b29-cfb0-4de3-ac24-772ebb6c4b76
# ╟─8cddf5fb-7b93-43c6-8c61-b9931d69fd83
# ╟─0fd1d551-49b0-4194-b668-fec7aa3ed9a7
# ╟─e2c520ec-2717-42c4-b4f0-0da56cf59927
# ╟─03619862-196e-419d-9c37-2e9ada04c84d
# ╟─55e55af6-f0c2-4fbd-96a5-dcc1bd2f1dd4
# ╠═9d69479a-e21b-4316-b46d-c36bbf634f0b
# ╟─2c67684f-392d-4a48-abba-be0fa8079344
# ╟─10ee39be-371b-487c-ac7c-9cbe60d2b367
# ╟─a127bc38-a03c-494a-a8c0-142d1674dc6c
# ╟─cc0c9d38-eb52-4f9e-87bf-70aeef1847db
# ╟─65bc7068-0cb1-4d54-9b81-8ec307893db8
# ╠═d88dddba-d9e1-4aaf-86ca-b107bb937fd8
# ╟─492b8dbf-3aac-4966-ba4c-03c78fea787d
# ╟─c3501463-f43a-465c-a81d-4bd5629793a3
# ╟─47ecff8a-c30d-40b5-87a0-324e39493f1e
# ╟─f7a1eb21-e615-4209-baa2-d682aefa1f1b
# ╟─fc60be23-b455-4803-930e-1738f70a3c76
# ╟─d1eb323f-f72f-4538-b4de-1c9541acf453
# ╟─440651dd-2ee5-46c0-a9ea-a5bce1093643
# ╟─60d407f7-0b36-4355-b4e0-ef0c6d83f800
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
