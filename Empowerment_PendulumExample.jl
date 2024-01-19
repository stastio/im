# run by 'julia> include("Empowerment_PendulumExample.jl")'

using Plots
using LinearAlgebra
using StatsBase
using Distributions 
using LaTeXStrings

ScalarType = Float64
ArrayType  = Array{ScalarType, 1}
Num = ScalarType

abstract type Model end

mutable struct SinglePendulum <: Model
	g::Float64
	l::Float64
	m::Float64
	Id
	empType::String
end


function SinglePendulum(; g=9.8, l=1.0, m=1.0, Id = diagm(0=>ones(2)), empType="EMP")
	SinglePendulum(g, l, m, Id, empType)
end

function RungeKutta(dyn::Model, f::Function, dt::ScalarType, x₀::ArrayType)

	k₁ = f(x₀        ,  dyn)*dt
	k₂ = f(x₀ + 0.5k₁,  dyn)*dt
	k₃ = f(x₀ + 0.5k₂,  dyn)*dt
	k₄ = f(x₀ +    k₃,  dyn)*dt

	x  = x₀ + (k₁+2k₂+2k₃+k₄)/6

end

function RandomShoot(dyn, 𝒟, 𝒱, 𝒫, dt, T, x₀, a₀)

	T̄  = Int(T/1)
	dt̄ = dt

	X = Array{ScalarType, 2}(undef, size(x₀, 1), T̄)
	X[:, 1] = x₀

	for i in 2:T̄
		X[:, i] = RungeKutta(dyn, 𝒟[:f], dt̄, X[:, i-1]) + 𝒟[:g](X[:, i-1], dyn)*a₀*dt̄
	end

	[X]

end


function Sensitivity(dyn::Model, f::Function, ∇f::Function, g::Function, x₀::ArrayType, dt::ScalarType, T::Int)

	if T == 1
		return g
	end

	M̃ = Vector{Array{ScalarType, 2}}(undef, T)
	g̃ = Vector{Array{ScalarType, 1}}(undef, T)

	xₜ  = x₀

	for n = 1:T-1
		M̃[n] = dyn.Id + ∇f(xₜ, dyn)*dt
		g̃[n] = g(xₜ, dyn)*dt
		xₜ = RungeKutta(dyn, f, dt, xₜ)
	end
	M̃[T] = dyn.Id
	g̃[T] = g(xₜ, dyn)*dt

	if dyn.empType == "EMP"
		M̂  = cumprod( M̃ |> reverse ) |> reverse
		M = hcat(map( t ->  M̂[t] * g̃[t], 1:T)...)
	elseif dyn.empType == "CEF"
		M̂  = cumprod( M̃ )
		M = vcat(map( t ->  M̂[t] * g̃[1], 1:T)...)
	elseif dyn.empType == "ALE"
		M̂  = cumprod(M̃)[end]
		M = M̂ * g̃[1]
	end
	M
end

function Gramian(dyn::Model, 𝒟::Dict, x₀::ArrayType, dt::ScalarType, T::Int)

	f, ∇f, g = 𝒟[:f], 𝒟[:∇f], 𝒟[:g]

	ℱ =  Sensitivity(dyn, f, ∇f, g, x₀, dt, T)

	ℱ*ℱ'

end

function Zᵢ(dyn::Model, 𝒟::Dict, x₀::ArrayType, dt::ScalarType, T::Int)

	𝒢 = Gramian(dyn, 𝒟, x₀, dt, T)

	E, V = eigen(𝒢)

	Ē = sqrt.(abs.(E))

	Z = Ē' .* V

end

function Σᵤ_ₓ(i::Union{Int, Nothing}, σ::Array{ScalarType, 1}, Zᵢ::Array{ScalarType, 2}, Ση::Array{ScalarType, 2})

	N = size(Zᵢ, 2)

	Ĩ = filter( x -> !isequal(x, i), 1:N)

	Z̃ = Zᵢ[:, Ĩ]

	if dyn.empType =="CEF"
		return σ[1]*Z̃  * Z̃' + Ση
	end

	σ̃ = σ[Ĩ]
	Z̃ * diagm(σ̃) * Z̃' + Ση

end


function Ση(𝒱::Dict, dyn::Model, 𝒟::Dict, x₀::ArrayType, dt::ScalarType, T::Int)

	f, ∇f, h = 𝒟[:f], 𝒟[:∇f], 𝒟[:h]

	𝒢 =  Sensitivity(dyn, f, ∇f, h, x₀, dt, T)

	# system noise: 𝒢*𝒱*𝒢'
	# system noise+ sensor noise: 𝒢*𝒱*𝒢' + Σν*𝒱
	#sensor noise
	𝒢*𝒱[:dynamics]*𝒢' + 𝒱[:sensor]*Σν
end

function waterfillingIterative(𝒱::Dict, P::ScalarType, dyn::Model, 𝒟::Dict, x₀::ArrayType, dt::ScalarType, T::Int)

	Σ̄η = Ση(𝒱, dyn, 𝒟, x₀, dt, T)
	Z̄ᵢ = Zᵢ(dyn, 𝒟, x₀, dt, T)

	N = length(x₀)
	σ² = P*ones(ScalarType, N)/N

	if dyn.empType == "CEF"
		Σ′ᵤ_ₓ = P*Z̄ᵢ * Z̄ᵢ' + Σ̄η
		return 0.5log(abs(det(Σ′ᵤ_ₓ))) - 0.5log(abs(det(Σ̄η)))
	end



	μ  = -1  
	h² = -1
	for k in 1:10

		Σ̄ᵤ_ₓ = [Σᵤ_ₓ(î, σ², Z̄ᵢ, Σ̄η)^(-1) for î in 1:N]

		h²  = vcat([ (Z̄ᵢ[:, i]' * Σ̄ᵤ_ₓ[i] * Z̄ᵢ[:, i])  for i in 1:N ]...)

		σ², μ = waterfilling(P, h²)

	end

	0.5log(abs(det(Σᵤ_ₓ(nothing, σ², Z̄ᵢ, Σ̄η) .+ 0))) - 0.5log(abs(det(Σ̄η)))

end

function waterfilling(P̃, h̃²::Vector{ScalarType})


	h² = 1.0 ./h̃²  
	sort!(h²)     

	P = Vector{Float64}(undef,length(h²))

	ν = h²                      # alias for ν-axis
	P[1] = 0.0

	for i in 2:length(h²)
		P[i] = P[i-1] + (i-1) * (ν[i]-ν[i-1])
	end

	# find the right ν

	bot = 1
	top = length(ν)+1           # ∞, not in array

	while top-bot > 1           # still to search
		mid = (bot+top) ÷ 2
		# mid is never outside of array, and is always at least 1
		# larger than bot and 1 smaller than top
		if P̃ >= P[mid]
			bot = mid
		else
			top = mid
		end
	end

	# now P̃ is between bot and top index, bot possibly included.
	# Compare with the bot and add the residual Δν

	ν̃ = ν[bot] + (P̃-P[bot])/bot

	# we now have our ν, compute the powers of each channel
	p = zeros(length(P))
	for i in 1:bot
		p[i] = ν̃-ν[i]
	end

	return p |> z -> z[end:-1:1], ν̃
end

function Landscape(dyn, 𝒟, 𝒱, 𝒫, dt, T, K::Int, stateSpaceBox)

	numDims = size(stateSpaceBox, 1)
	Δ 	= diff(stateSpaceBox, dims=2) / K |> vec
	E   	= Array{ScalarType, numDims}(undef, ntuple(i->K, numDims)...)

	Threads.@threads for i in CartesianIndices( ntuple(i->K, numDims) )
		xᵢ = stateSpaceBox[:, 1] .+  (i.I .- 1) .* Δ
		E[i] = waterfillingIterative(𝒱, 𝒫, dyn, 𝒟, xᵢ, dt, T)
	end

	E
end


include("est_max.jl")

function intrinsic_controller(dyn, 𝒟, 𝒱, 𝒫, dt, T, dt̄, T̄, xᵢ, M, repeatedAction, Amax)
	Amax = Amax
	X = Array{ScalarType, 2}(undef, xᵢ |> length, T̄*repeatedAction)
	X̃  = Vector{Array{ScalarType, 2}}(undef, 0) 
	E′ = []
	ã = 0.0
	Ã = Array{ScalarType, 2}(undef, 1, T̄*repeatedAction)
	for t̄ in 1:T̄
		X̄ = vcat([RandomShoot(dyn, 𝒟, 𝒱, 𝒫, dt, T̄, xᵢ, aᵢ) for aᵢ in [-Amax, 0, +Amax]]...)

		Ē = hcat([(x̄ -> waterfillingIterative(𝒱, 𝒫, dyn, 𝒟, x̄[:, i], dt, T)).(X̄) for i ∈ T]...)

		E′ = Ē[:, end]

		pairs = Pairs(-Amax, E′[1], 0, E′[2], +Amax, E′[3])

		ã   = parabola_max_interval(pairs)

		foreach(x̃ -> push!(X̃, x̃), X̄)

		for k in 1:repeatedAction
			X[:, (t̄-1)*repeatedAction + k] = xᵢ
			xᵢ = RungeKutta(dyn, 𝒟[:f], dt̄, xᵢ)+ 𝒟[:g](xᵢ, dyn)*ã*dt̄
		end
		Ã[(t̄-1)*repeatedAction+1:t̄*repeatedAction] .= ã*ones(repeatedAction)
	end
	X, Ã, X̃
end


### load dynamics  ###
include("SinglePendulum.jl")
dyn = SinglePendulum(; empType="EMP")
stateSpaceBox = [-2pi/1 +2pi/1; -8pi 8pi]
Σν =  diagm(0=>ones(2))

### init for control ###
x₀ = [0.0, 0.0]

### grid size for plots ###
K = 55*7 



# empowerment params
Ttotal = 0.50
dt = 1e-1 
T = Int(floor(Ttotal/dt))
Energy = 5.0

### pole length ###
L= 1.0
### power ###
𝒫 = 2*Energy/(Ttotal^2)###
### noise variance ###
𝓥 = Dict(:dynamics=>L^2, :sensor=>1/L^2)
𝒟 = Dict(:f=>f, :∇f=>∇f, :h=>h, :g=>g, "type"=>:pend)

# control params 
T̄ = 1500 # control trajectory length
M = 2  # number of random shooting 
dt̄ = 0.01 
repeatedAction = 1
### maximal action ###
Amax = 1.0

X, Ã, X̃ = intrinsic_controller(dyn, 𝒟, 𝓥, 𝒫, dt, T, dt̄, T̄, x₀, M, repeatedAction, Amax)

### plots ###

θaxis = range(stateSpaceBox[1, :]..., length=K)
θ̇axis = range(stateSpaceBox[2, :]..., length=K)
E = Landscape(dyn, 𝒟, 𝓥, 𝒫, dt, T, K, stateSpaceBox); 

p0 = plot(heatmap(θaxis, θ̇axis, E'), legend=:none, colorbar=:none,
xticks = ([θaxis[1], θaxis[(K+1)/4 |> ceil |> Int ],  θaxis[(K+1)/2 |> Int], θaxis[3(K+1)/4 |> ceil |> Int ],   θaxis[end]], ["-2π", "π", 0, "π", "2π"]), xtickfontsize=16,
yticks = ([θ̇axis[1],  θ̇axis[(K+1)/2 |> Int], θ̇axis[end]], ["-5π", "0", "5π"]), ytickfontsize=16, 
ylabel = L"\dot{\!\!\!\!\theta}\si{\rad}"*" in rad/sec", xlabel=L"\theta"*" in rad", labelfontsize=16)
plot!(size=(1400,1000))

h0 = plot(heatmap(θaxis, θ̇axis, E'), legend=:none, 
xticks = ([θaxis[1], θaxis[(K+1)/4 |> ceil |> Int ],  θaxis[(K+1)/2 |> Int], θaxis[3(K+1)/4 |> ceil |> Int ],   θaxis[end]], ["-2π", "-π (top)", "0 (bottom)", "π (top)", "2π"]), xtickfontsize=18,
yticks = ([θ̇axis[1],  θ̇axis[(K+1)/2 |> Int], θ̇axis[end]], ["-8π", "0", "8π"]), ytickfontsize=18, right_margin=5Plots.mm, left_margin=5Plots.mm)
plot!(size=(1600,1400))


# X̃_ = reshape(X̃, (3, :)) 
# for i in 1:size(X̃_, 2)
# 	plot!( X̃_[1, i][1, 1:2:end], X̃_[1, i][2, 1:2:end], lc=:green, lw=2, ms=1, msa=0, msw=0)
# 	plot!( X̃_[3, i][1, 1:2:end], X̃_[3, i][2, 1:2:end], lc=:red, lw=2, ms=1, msa=0, msw=0)
# 	plot!( X̃_[2, i][1, 1:2:end], X̃_[2, i][2, 1:2:end], lc=:white, lw=2, ms=1, msa=0, msw=0)
# 
# end

plot!(p0, X[1, :], X[2, :], lc=:white, lw=2, mc=:white, ms=2, label="Intrinsic Trajectory")
plot!(h0, X[1, :], X[2, :], lc=:white, lw=2, mc=:white, ms=2, label="Intrinsic Trajectory")
plot!(h0, [X[1, 1]], [X[2, 1]], lt=:scatter, mc=:red, ms=7, markershape=:circle)
plot!(h0, [X[1, end]], [X[2, end]], lt=:scatter, mc=:green, ms=7, markershape=:circle)
p01 = plot( X[1, :], X[2, :], lc=:black, lw=1, mc=:black, ms=2, label="Traj")

p1 = plot( 1:(T̄*repeatedAction), [X[1, :] X[2, :]], label=["θ" "θ̇"])
p2 = plot( Ã |> vec, lt=:scatter, msw=0, mca=0, ms=3, color=:red, labelfontsize=18, tickfontsize=18, label="Control, Ē=$(round(sum(Ã .^ 2)/length(Ã), digits=3))")
h1 = plot((1:length(Ã))*dt̄, Ã |> vec, yticks=([-1, 0, +1], ["-1", "0", "+1"]), markershape=:circle, lt=:scatter, msc=:match, msw=0, mca=0, ms=3, color=:black, label=:none, xlabel=L"t"*" in sec", labelfontsize=18, tickfontsize=18)
plot( h0, h1, size=(1100, 1100), layout = @layout [a{0.8h}; b])

savefig("Empowerment_PendulumExample.png")
