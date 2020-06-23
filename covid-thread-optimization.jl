using CSV, DataFrames, Query
using Statistics, Dates, TimeSeries
using BlackBoxOptim

"""
    get_data(file, county)

    Reads the CSV file of all COVID cases defined in `file`, and 
    returns a daily aggregation from the county defined in `county`.
"""
function get_data(file, county)
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

function calcRMSE(truth, vals)
    return sqrt(sum( (truth .- vals).^2 ))
end

#get the targets (Orange county values)
case_data = get_data("Florida COVID Case Line 20200602.csv", "Orange")

include("covid-model.jl")

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

function run_model(params)
    steps = 2400;
    dt = Dates.Hour(1);
    N = 1000;

    adata = [:status, :hospitalized, :quarantined]
    mdata = [:tick, nagents]

    Kβ_min = params[11] < params[12] ? params[11] : params[12]
    Kβ_max = params[11] < params[12] ? params[12] : params[11]

    model50ig = initialize(extents=(params[1], params[2]), N=N, dt=dt, speed=params[3], percIgnore=params[4], 
        interaction_radius=params[5], detection_threshold=params[6], percInfected=params[7], critical_threshold=params[8],
        mask_effect=params[9], mask_perc=params[10], Kβ_min=params[11],
        Kβ_max=params[12], βmax=params[13])
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
    modeled[!,:detect] .*= mod_scale;
    modeled[!,:hosp] .*= mod_scale;
    modeled[!,:deaths] .*= mod_scale;

    return modeled
end

function evaluate(params)

    modeled = run_model(params)

    detect_err = calcRMSE(case_data[:ts].detected, modeled.detect)
    hosp_err = calcRMSE(case_data[:ts].hosp_cases, modeled.hosp)
    death_err = calcRMSE(case_data[:ts].deaths, modeled.deaths)

    return (detect_err, hosp_err, death_err)
end;

"""
params: x_extent, y_extent, 
speed, percIgnore, 
interaction_radius, detection_threshold,
percInfected, critical_threshold,
mask_effect, mask_perc,
reinfect_prob, Kβ_min,
Kβ_max, βmax
"""
range = [ (1., 5.), (1., 5.), 
        (0.001, 1.0), (0.01, 0.99), 
        (0.0001, 0.2), (5., 50.),
        (0.01, 0.2), (40., 80.),
        (0.1, 0.8), (0.1, 0.95),
        (0.001, 0.1), (1., 7.,), 
        (3., 10.), (40., 99.)]

    Method=:borg_moea,
    FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
    MaxTime=28800,
    TraceMode=:verbose,
    TraceInterval=1.,
    Workers=workers()
);

res = bboptimize(evaluate, 
SearchRange=range, 
Method=:borg_moea, 
FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
TraceMode=:verbose,
TraceInterval=1.,
MaxTime=600,
NThreads=Threads.nthreads()-1
)

modeled = run_model(best_candidate(res))
best_candidate(res)
#evaluate([1,1,0.02,0.5,0.012,20])
evaluate([1,1,0.02,0.5,0.012,20,0.01,40,0.5,0.3,0.001,5,10,75])