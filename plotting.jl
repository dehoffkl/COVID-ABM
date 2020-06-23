using Plots

gr()

plot(result[:date], result[:detect_mean], ribbon=[result[:detect_z], result[:detect_z]], alpha=0.4, label="Modeled Detected Cases");
plot!(case_data[:ts][!,:timestamp], case_data[:ts][!,:detected], label="Smoothed Actual Detected Cases")

plot(result[:date], result[:hosp_mean], ribbon=[result[:hosp_z], result[:hosp_z]], alpha=0.4, label="Modeled Hospitalizations");
plot!(case_data[:ts][!,:timestamp], case_data[:ts][!,:hosp_cases], label="Smoothed Actual Hospitalizations")

plot(result[:date], result[:deaths_mean], ribbon=[result[:deaths_z], result[:deaths_z]], alpha=0.4, label="Modeled Deaths");
plot!(case_data[:ts][!,:timestamp], case_data[:ts][!,:deaths], label="Smoothed Actual Deaths")

plot(result[:date], result[:detect], label="Modeled Detected Cases")
