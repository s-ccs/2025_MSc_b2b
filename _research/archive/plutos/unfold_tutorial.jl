### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 3bb6f77e-70da-49d0-8d24-e28f8971be2f
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 3d527906-09a9-49eb-b6d1-fbda847d7f75
begin
	using DataFrames
	using Unfold
	using UnfoldMakie, CairoMakie # for plotting
	using UnfoldSim
end

# ╔═╡ 4d87d74d-a718-49f3-ae9a-c8f9db186c67
data, evts = UnfoldSim.predef_eeg()

# ╔═╡ d481dc8b-a247-4c79-abd1-6b999e9d2871
evts

# ╔═╡ 63daa402-4ff6-4e8a-941a-1e14df2b8c54
typeof(data)

# ╔═╡ 10175e71-eada-45b2-b11f-3454a909a0cf
size(data)

# ╔═╡ c8613e73-22aa-4605-b09e-6f8ecbd16d53
data[1:20]

# ╔═╡ 07c4bf5d-ad5d-4f10-aa69-a166f9e02331
begin
	times_cont = range(0,length=200,step=1/100) # we simulated with 100hz for 0.5 seconds
	
	fig,ax,h = plot(times_cont,data[1:200])
	vlines!(evts[evts.latency .<= 200, :latency] ./ 100;color=:black) # show events, latency in samples!
	ax.xlabel = "time [s]"
	ax.ylabel = "voltage [µV]"
	fig
end

# ╔═╡ 0dd2ae56-1b47-4f05-a8b2-9b783cc34a05
begin
	# Unfold supports multi-channel, so we could provide matrix ch x time, which we can create like this from a vector:
	data_r = reshape(data, (1,:))
	# cut the data into epochs
	data_epochs, times = Unfold.epoch(data = data, tbl = evts, τ = (-0.4, 0.8), sfreq = 100); # channel x timesteps x trials
	size(data_epochs)
end

# ╔═╡ 448f9999-0b88-4553-86a0-dbbc32e9cfd0
f = @formula 0 ~ 1 + condition + continuous # note the formulas left side is `0 ~ ` for technical reasons`

# ╔═╡ 78f641cb-b800-4998-9aa9-49fc499d8743
m

# ╔═╡ cb4e3746-0f45-4c44-ac67-34a57b47f24c


# ╔═╡ 47b53e29-246e-4d4a-aa4b-88d1c2fc46d2
first(coeftable(m), 6)

# ╔═╡ 5598be42-54ac-4966-ab0e-7bf7d7dce73a
basisfunction = firbasis(τ=(-0.4,.8),sfreq=100)

# ╔═╡ 8518b019-03b4-4fcf-88a4-6c1e5cc07eb3
bf_vec = [Any=>(f,basisfunction)]

# ╔═╡ a10b8585-83d9-436e-921c-76a275efac70
m_fir = fit(UnfoldModel,bf_vec,evts,data;)

# ╔═╡ 24fbc058-2cbc-49fa-b38b-3599b07687cd


# ╔═╡ 912b1799-2be2-49bd-b5d6-e88eb63fa7d5
# ╠═╡ disabled = true
#=╠═╡
m = fit(UnfoldModel, f, evts, data_epochs, times);
  ╠═╡ =#

# ╔═╡ 6588229f-3cde-40e7-a2b9-e914899dba9f
begin
	results = coeftable(m_fir)
	plot_erp(results)
end

# ╔═╡ af2e0d7d-1f22-4a5f-a252-b194b97e940d
# ╠═╡ disabled = true
#=╠═╡
begin
	results = coeftable(m)
	plot_erp(results)
end
  ╠═╡ =#

# ╔═╡ cac3aadb-9424-4d81-9614-01c9358b10d3
m = fit(UnfoldModel, [Any=>(f, times)], evts, data_epochs);

# ╔═╡ Cell order:
# ╠═3bb6f77e-70da-49d0-8d24-e28f8971be2f
# ╠═3d527906-09a9-49eb-b6d1-fbda847d7f75
# ╠═4d87d74d-a718-49f3-ae9a-c8f9db186c67
# ╠═d481dc8b-a247-4c79-abd1-6b999e9d2871
# ╠═63daa402-4ff6-4e8a-941a-1e14df2b8c54
# ╠═10175e71-eada-45b2-b11f-3454a909a0cf
# ╠═c8613e73-22aa-4605-b09e-6f8ecbd16d53
# ╠═07c4bf5d-ad5d-4f10-aa69-a166f9e02331
# ╠═0dd2ae56-1b47-4f05-a8b2-9b783cc34a05
# ╠═448f9999-0b88-4553-86a0-dbbc32e9cfd0
# ╠═912b1799-2be2-49bd-b5d6-e88eb63fa7d5
# ╠═cac3aadb-9424-4d81-9614-01c9358b10d3
# ╠═78f641cb-b800-4998-9aa9-49fc499d8743
# ╠═cb4e3746-0f45-4c44-ac67-34a57b47f24c
# ╠═47b53e29-246e-4d4a-aa4b-88d1c2fc46d2
# ╠═af2e0d7d-1f22-4a5f-a252-b194b97e940d
# ╠═5598be42-54ac-4966-ab0e-7bf7d7dce73a
# ╠═8518b019-03b4-4fcf-88a4-6c1e5cc07eb3
# ╠═a10b8585-83d9-436e-921c-76a275efac70
# ╠═6588229f-3cde-40e7-a2b9-e914899dba9f
# ╠═24fbc058-2cbc-49fa-b38b-3599b07687cd
