using Distributed
@everywhere using ProgressMeter

procs = addprocs(5)

@everywhere using CSV, DataFrames, Query
@everywhere using Statistics, Dates, TimeSeries
@everywhere using BlackBoxOptim

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

@everywhere function calcRMSE(truth, vals)
    return sqrt(mean( (truth .- vals).^2 ))
end

#get the targets (Orange county values)
@everywhere case_data = get_data("Florida COVID Case Line 20200602.csv", "Orange")

@everywhere include("covid-model.jl")

"""
ext = (1,1);
N = 1000;
dt = Dates.Hour(1);
speed = 0.02;
percIgnore = 0.5;
int_rad = 0.012;
detect_thresh = 20;
steps = 2400;
"""

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

    mod_scale = case_data[:cases] / sum(data_agg.detect)

    modeled = data_agg;
    #modeled[!,:detect] .*= mod_scale;
    #modeled[!,:hosp] .*= mod_scale;
    #modeled[!,:deaths] .*= mod_scale;

    return modeled
end

@everywhere function evaluate(params)

    modeled = run_model(params, rseed=42)

    detect_err = calcRMSE(case_data[:ts].detected, modeled.detect)
    hosp_err = calcRMSE(case_data[:ts].hosp_cases, modeled.hosp)
    death_err = calcRMSE(case_data[:ts].deaths, modeled.deaths)

    return detect_err + 10.0*hosp_err + 20.0*death_err
end;

"""
params: extent, 
speed, percIgnore, 
interaction_radius, detection_threshold,
percInfected, critical_threshold,
mask_effect, mask_perc,
reinfect_prob, Kβ_min,
Kβ_max, βmax
"""
setup = bbsetup(evaluate;
    SearchRange = [ (0.5, 5.), 
        (0.0001, 5.), (0.0, 1.0), 
        (1.0E-10, 1.0), (5., 50.),
        (1.0E-10, 0.25), (25., 95.),
        (0.0, 1.0), (0.0, 1.0),
        (1E-10, 0.5), (1., 7.,), 
        (3., 10.), (20., 100.)],
    Method=:dxnes,
    MaxTime=28800,
    Workers=workers()
);

res = bboptimize(setup)

n = 50
bp = best_candidate(res)
best_params = []
for i in 1:n
    push!(best_params, bp)
end

f = open("bboptim-fitted.txt", "w")
for i in bp
    println(f, i)
end
close(f)

bunched = @showprogress pmap(x->run_model(x), best_params)

test = DataFrame()
for i in 1:length(bunched)
    append!(test, bunched[i])
end
result = combine(groupby(test, :date), nrow, :detect=>mean, :detect=>std, :hosp=>mean, :hosp=>std, :deaths=>mean, :deaths=>std)
result[!, :detect_z] = 1.96.*result[:,:detect_std]./sqrt.(result[:,:nrow])
result[!, :hosp_z] = 1.96.*result[:,:hosp_std]./sqrt.(result[:,:nrow])
result[!, :deaths_z] = 1.96.*result[:,:deaths_std]./sqrt.(result[:,:nrow])

CSV.write("bboptim-fitted_runs.csv", result)



modeled = run_model(best_candidate(res))
best_candidate(res)
#evaluate([1,1,0.02,0.5,0.012,20])
run_model([1,1,0.02,0.5,0.012,20,0.01,40,0.5,0.3,0.001,5,10,75], rseed=42)