using Agents
using BenchmarkTools
using Distributions
using AgentsPlots
using Plots
using NaNMath; nm=NaNMath

include("covid-model.jl")

#initialize potential model parameters
model25ig = initialize(extents=(1, 1), N=1000, dt=1.0/24, speed=.02, percIgnore=0.25, interaction_radius=.012, detection_threshold=20)
model50ig = initialize(extents=(1, 1), N=5000, dt=1.0/24, speed=.02, percIgnore=0.5, interaction_radius=.012, detection_threshold=20)
model10ig = initialize(extents=(1, 1), N=1000, dt=1.0/24, speed=.02, percIgnore=0.1, interaction_radius=.012, detection_threshold=20)

#run potential model parameters
data1, mdata1 = run!(model25ig, agent_step!, model_step!, 2400, adata=adata, mdata=mdata);
mdata1[!, :deceased] = map(n -> 1000 - n, mdata1[:, :nagents]);
@time data2, mdata2 = run!(model50ig, agent_step!, model_step!, 2400, adata=adata, mdata=mdata);
mdata2[!, :deceased] = map(n -> 5000 - n, mdata2[:, :nagents]);
data3, mdata3 = run!(model10ig, agent_step!, model_step!, 2400, adata=adata, mdata=mdata);
mdata3[!, :deceased] = map(n -> 1000 - n, mdata3[:, :nagents]);

last(mdata2[:, :deceased]) / sum(data2[:, :detected_quarantined])
nm.mean(data2[:, :detected_quarantined] ./ data2[:, :infected_status])


#produce graphs of detected cases
gr()

p = plot(
    data1[:, aggname(:quarantined, detected)],
    label="25% ignore stay-at-home",
    ylabel = "Infected",
)
plot!(p, data2[:, aggname(:quarantined, detected)], label="50% ignore stay-at-home")
plot!(p, data3[:, aggname(:quarantined, detected)], label="10% ignore stay-at-home")
vline!([21*24, 42*24, 64*24, 86*24], label="")

#product graphs of covid deaths
p = plot(
    mdata1[:, :deceased],
    label="25% ignore stay-at-home",
    ylabel = "Deaths",
)
plot!(p, mdata2[:, :deceased], label="50% ignore stay-at-home")
plot!(p, mdata3[:, :deceased], label="10% ignore stay-at-home")

