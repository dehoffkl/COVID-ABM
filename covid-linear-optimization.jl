using Distributed

procs = addprocs(5)

@everywhere using CSV, DataFrames, Query
@everywhere using Statistics, Dates, TimeSeries
@everywhere using Distributions
@everywhere using ProgressMeter
@everywhere using GLMNet



"""
    get_data(file, county)

    Reads the CSV file of all COVID cases defined in `file`, and 
    returns a daily aggregation from the county defined in `county`.
"""
@everywhere function get_data(file, county)
    #Read all cases from FL
    FL_cases = CSV.read(file);

    #filter to Orange County
    cases = FL_cases |> 
        @filter(_.County == county) |>
        @map({_.EDvisit, _.Hospitalized, _.Died, _.Contact, _.Case1, _.EventDate, _.ChartDate}) |>
        DataFrame;
    cases[!, :HospBool] = cases[:, :Hospitalized].=="YES";
    cases[!, :DiedBool] = cases[:, :Died].=="Yes";
    cases[!, :Date] = Date.(cases[:, "EventDate"], Dates.DateFormat("y/m/d H:M:S"));
    cases=combine(groupby(cases, :Date), nrow=>:detected, :HospBool=>sum=>:hosp_cases, :DiedBool=>sum=>:deaths)
    cases[!, :detected] = cases[:, :detected].*1.0;
    cases[!, :hosp_cases] = cases[:, :hosp_cases].*1.0;
    cases[!, :deaths] = cases[:, :deaths].*1.0;
    sort!(cases);

    tot_deaths = sum(cases[:, :deaths])
    tot_hosp = sum(cases[:, :hosp_cases])
    tot_cases = sum(cases[:, :detected])
    mort = tot_deaths / tot_cases

    ts = TimeArray(cases, timestamp=:Date)
    ts = moving(mean, ts, 10)

    return Dict(:ts=>DataFrame(ts), :deaths=>tot_deaths, :cases=>tot_cases,
        :hosp=>tot_hosp, :mortality=>mort) 
end

#get the targets (Orange county values)
@everywhere case_data = get_data("Florida COVID Case Line 20200602.csv", "Orange")

@everywhere include("covid-model.jl")

@everywhere function run_model(params; rseed=-1)
    steps = 2400;
    dt = Dates.Hour(1);
    N = 1000;

    if rseed < 0
        rseed = rand(1:1000000)
    end

    adata = [:status, :hospitalized, :quarantined]
    mdata = [:tick, nagents]

    Kβ_min = params[11] < params[12] ? params[11] : params[12]
    Kβ_max = params[11] < params[12] ? params[12] : params[11]

    model50ig = initialize(extents=(params[1], params[1]), N=N, dt=dt, speed=params[2], percIgnore=params[3], 
        interaction_radius=params[4], detection_threshold=params[5], percInfected=params[6], critical_threshold=params[7],
        mask_effect=params[8], mask_perc=params[9], reinfect_prob=params[10], Kβ_min=Kβ_min,
        Kβ_max=Kβ_max, βmax=params[13], initial_qlevel=2, rseed=rseed)
    adata_out, mdata_out = run!(model50ig, agent_step!, model_step!, steps, adata=adata, mdata=mdata);

    #determine the total number of deaths
    mdata_out[!, :deceased] = map(n -> 1000 - n, mdata_out[:, :nagents]);

    #create a date column for aggregation
    mdata_out[!, :date] = convert.(Date, mdata_out[:,:tick]);

    #determine where new infections, recoveries, hospitalizations, and detections fall
    sort!(adata_out, [:id, :step])
    adata_out[!, :new_infected] .= false
    adata_out[!, :new_recovered] .= false
    adata_out[!, :new_hosp] .= false
    adata_out[!, :new_detect] .= false

    prev_stat = adata_out[1, :status];
    prev_hosp = adata_out[1, :hospitalized];
    prev_detected = adata_out[1, :quarantined];
    for i in 2:nrow(adata_out)
        #check for id change
        if adata_out[i, :id] != adata_out[i-1, :id]
            prev_stat = adata_out[i, :status];
            prev_hosp = adata_out[i, :hospitalized];
            prev_detected = adata_out[i, :quarantined];
            continue
        end

        if adata_out[i, :status] == :I && prev_stat != :I
            adata_out[i, :new_infected] = true
        end
        if adata_out[i, :status] == :R && prev_stat != :R
            adata_out[i, :new_recovered] = true
        end
        if adata_out[i, :hospitalized] && !prev_hosp
            adata_out[i, :new_hosp] = true
        end
        if adata_out[i, :quarantined] && !prev_detected
            adata_out[i, :new_detect] = true
        end
        prev_stat = adata_out[i, :status]
        prev_hosp = adata_out[i, :hospitalized]
        prev_detected = adata_out[i, :quarantined]
    end;

    #join with step-level data and aggregate to the daily level
    adata_out = innerjoin(adata_out, mdata_out, on=:step);
    data_agg = combine(groupby(adata_out, :date),
        :new_infected=>sum=>:inf,
        :new_recovered=>sum=>:rec,
        :new_hosp=>sum=>:hosp,
        :new_detect=>sum=>:detect,
        :deceased=>maximum=>:tot_deceased
    );

    #determine cumulative total of infections, recoveries, hospitalizations, and daily count of deaths
    data_agg[!, :tot_inf] .= 0
    data_agg[!, :tot_rec] .= 0
    data_agg[!, :tot_hosp] .= 0
    data_agg[!, :tot_detect] .= 0
    data_agg[!, :deaths] .= 0
    data_agg[1, :tot_inf] = data_agg[1, :inf]
    data_agg[1, :tot_rec] = data_agg[1, :rec]
    data_agg[1, :tot_hosp] = data_agg[1, :hosp]
    data_agg[1, :tot_detect] = data_agg[1, :detect]
    data_agg[1, :deaths] = data_agg[1, :tot_deceased]
    for i in 2:nrow(data_agg)
        data_agg[i, :tot_inf] = data_agg[i-1, :tot_inf] + data_agg[i, :inf]
        data_agg[i, :tot_rec] = data_agg[i-1, :tot_rec] + data_agg[i, :rec]
        data_agg[i, :tot_hosp] = data_agg[i-1, :tot_hosp] + data_agg[i, :hosp]
        data_agg[i, :tot_detect] = data_agg[i-1, :tot_detect] + data_agg[i, :detect]
        data_agg[i, :deaths] = data_agg[i, :tot_deceased] - data_agg[i-1, :tot_deceased]
    end;

    max_case_dt = maximum(case_data[:ts].timestamp)
    data_agg = filter(row->row[:date] <= max_case_dt, data_agg)

    #mod_scale = case_data[:cases] / sum(data_agg.detect)

    modeled = data_agg;
    #modeled[!,:detect] .*= mod_scale;
    #modeled[!,:hosp] .*= mod_scale;
    #modeled[!,:deaths] .*= mod_scale;

    return modeled
end

@everywhere function evaluate(params)
    id = popfirst!(params)
    test = run_model(params)

    return(format_output(test, id=id))
end

@everywhere function format_output(test; id=-1.)
    test2 = test[!, [:date, :detect, :hosp, :deaths]];
    test2[!, :detect] = Float64.(coalesce(test2[!, :detect],0))
    test2[!, :hosp] = Float64.(coalesce(test2[!, :hosp],0))
    test2[!, :deaths] = Float64.(coalesce(test2[!, :deaths],0))
    test2[!, :detectxhosp] = test2[!,:detect] .* test[!,:hosp];
    test2[!, :detectxdeaths] = test2[!,:detect] .* test[!,:deaths];
    test2[!, :hospxdeaths] = test2[!,:hosp] .* test[!,:deaths];
    test2[!, :detect_hosp] = test2[!,:detect] .+ test[!,:hosp];
    test2[!, :detect_deaths] = test2[!,:detect] .+ test[!,:deaths];
    test2[!, :hosp_deaths] = test2[!,:hosp] .+ test[!,:deaths];
    test2 = stack(test2, Not(:date))
    test2[!, :variable] = Dates.format.(test2[!, :date], "yyyy-mm-dd") .* "_" .* test2[!, :variable]
    test2[!, :id] .= id;
    for i in 1:nrow(test2)
        test2[i, :value] += rand(Normal(0, 0.01))
    end
    test2 = select!(test2, Not(:date));

    test2 = unstack(test2, :id, :variable, :value)
    return disallowmissing(test2)
end

#Set limits on ABM parameters
RandExtent = Uniform(0.5, 5);
RandSpeed = Uniform(0.0001, 5.0);
RandPerc = Uniform(0.0, 1.0);
RandRad = Uniform(1.e-10, 1.0);
RandDetect = Uniform(5.0, 50.0);
RandCrit = Uniform(25.0, 95.0);
RandInf = Uniform(1.0e-10, 0.25);
RandReinf = Uniform(1.0e-10, 0.50);
RandKβ_min = Uniform(1.0, 7.0);
RandKβ_max = Uniform(3.0, 10.0);
Randβ_max = Uniform(20.0, 100.0);

"""
params: extent, 
speed, percIgnore, 
interaction_radius, detection_threshold,
percInfected, critical_threshold,
mask_effect, mask_perc,
reinfect_prob, Kβ_min,
Kβ_max, βmax
"""
input = [rand(RandExtent), rand(RandExtent), rand(RandSpeed),
    rand(RandPerc), rand(RandRad), rand(RandDetect),
    rand(RandInf), rand(RandCrit), rand(RandPerc),
    rand(RandPerc), rand(RandReinf), rand(RandKβ_min),
    rand(RandKβ_max), rand(Randβ_max)];


n = 1000
plist = []
for i in 1:n
    push!(plist, [i, rand(RandExtent), rand(RandSpeed),
    rand(RandPerc), rand(RandRad), rand(RandDetect),
    rand(RandInf), rand(RandCrit), rand(RandPerc),
    rand(RandPerc), rand(RandReinf), rand(RandKβ_min),
    rand(RandKβ_max), rand(Randβ_max)])
end

param_frame = DataFrame(Matrix(DataFrame(plist))')
names!(param_frame, [:id, :ext, :speed, :percIgnore, :int_rad, :det_thresh, :percInit, :crit_thresh,
                      :mask_eff, :mask_perc, :reinf_prob, :Kβ_min, :Kβ_max, :β_max])

CSV.write("inputs.csv", param_frame)

bunched = @showprogress pmap(x->evaluate(x), plist);
output = DataFrame()
for i in 1:length(bunched)
    append!(output, bunched[i]);
end
CSV.write("simulations.csv", output)

y = param_frame[!, :crit_thresh]
X = disallowmissing(output[!, Not(:id)])

models = Dict()
for i in 2:ncol(param_frame)
    println("Fitting ", names(param_frame)[i])
    y = param_frame[!, i]
    cv = glmnetcv(Matrix(X), y, alpha=0.8, parallel=true)
    models[names(param_frame)[i]] = Dict("beta"=>cv.path.betas[:,argmin(cv.meanloss)],
                                        "int"=>cv.path.a0[argmin(cv.meanloss)])
end

cd = DataFrame(case_data[:ts]);
names!(cd, [:date, :detect, :hosp, :deaths])
exist = Matrix(format_output(cd)[!, Not(:id)])

"""
params: extent, 
speed, percIgnore, 
interaction_radius, detection_threshold,
percInfected, critical_threshold,
mask_effect, mask_perc,
reinfect_prob, Kβ_min,
Kβ_max, βmax
"""
fitted = [
    models["ext"]["int"] + (exist*models["ext"]["beta"])[1],
    models["speed"]["int"] + (exist*models["speed"]["beta"])[1],
    models["percIgnore"]["int"] + (exist*models["percIgnore"]["beta"])[1],
    models["int_rad"]["int"] + (exist*models["int_rad"]["beta"])[1],
    models["det_thresh"]["int"] + (exist*models["det_thresh"]["beta"])[1],
    models["percInit"]["int"] + (exist*models["percInit"]["beta"])[1],
    models["crit_thresh"]["int"] + (exist*models["crit_thresh"]["beta"])[1],
    models["mask_eff"]["int"] + (exist*models["mask_eff"]["beta"])[1],
    models["mask_perc"]["int"] + (exist*models["mask_perc"]["beta"])[1],
    models["reinf_prob"]["int"] + (exist*models["reinf_prob"]["beta"])[1],
    models["Kβ_min"]["int"] + (exist*models["Kβ_min"]["beta"])[1],
    models["Kβ_max"]["int"] + (exist*models["Kβ_max"]["beta"])[1],
    models["β_max"]["int"] + (exist*models["β_max"]["beta"])[1]
]

for i in 1:length(fitted)
    if fitted[i] < 0.0
        i == 6 ? fitted[i] = 0.001 : fitted[i] = 0.0
    end
end

f = open("linear-fitted.txt", "w")
for i in fitted
    println(f, i)
end
close(f)

fitted = readdlm("linear-fitted.txt", '\n')


n = 50
test_params = []
for i in 1:n
    push!(test_params, fitted)
end

bunched = @showprogress pmap(x->run_model(x), test_params)

test = DataFrame()
for i in 1:length(bunched)
    append!(test, bunched[i])
end
result = combine(groupby(test, :date), nrow, :detect=>mean, :detect=>std, :hosp=>mean, :hosp=>std, :deaths=>mean, :deaths=>std)
result[!, :detect_z] = 1.96.*result[:,:detect_std]./sqrt.(result[:,:nrow])
result[!, :hosp_z] = 1.96.*result[:,:hosp_std]./sqrt.(result[:,:nrow])
result[!, :deaths_z] = 1.96.*result[:,:deaths_std]./sqrt.(result[:,:nrow])

CSV.write("fitted_runs.csv", result)