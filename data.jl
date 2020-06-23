using CSV
using DataFrames
using Query
using Dates
using TimeSeries
using Statistics

using Plots

FL_cases = CSV.read("Florida COVID Case Line 20200602.csv");
orange_cases = FL_cases |> 
    @filter(_.County == "Orange") |>
    @map({_.EDvisit, _.Hospitalized, _.Died, _.Contact, _.Case1, _.EventDate, _.ChartDate}) |>
    DataFrame;
orange_cases[!, :HospBool] = orange_cases[:, :Hospitalized].=="YES";
orange_cases[!, :DiedBool] = orange_cases[:, :Died].=="Yes";
orange_cases[!, :Date] = Date.(orange_cases[:, "EventDate"], Dates.DateFormat("y/m/d H:M:S"));
orange_cases=combine(groupby(orange_cases, :Date), nrow=>:detected, :HospBool=>sum=>:hosp_cases, :DiedBool=>sum=>:deaths)
orange_cases[!, :detected] = orange_cases[:, :detected].*1.0;
orange_cases[!, :hosp_cases] = orange_cases[:, :hosp_cases].*1.0;
orange_cases[!, :deaths] = orange_cases[:, :deaths].*1.0;
sort!(orange_cases);
cases = TimeArray(orange_cases, timestamp=:Date)
cases = moving(mean, cases, 10)

gr()
plot(cases.detected, label="Known cases")
plot!(cases.hosp_cases, label="Hospitalized")
plot!(cases.deaths, label="Deaths")
