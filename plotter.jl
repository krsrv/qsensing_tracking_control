using CSV;
using DataFrames;

# Only 60 seconds to create output for vy
df = CSV.read("vmax_output.csv", DataFrame);
rename!(df,[:t1,:t2,:detuning,:vx,:detuned_vy_max,:ramsey_vy_max,:detuned_vy_end,:ramsey_vy_end,:sim_time]);

df.detuning_log = log10.(df.detuning);
df.t1 = ceil.(df.t1*10^6);
df.t2 = ceil.(df.t2*10^6);
df.t1t2 = [(df.t1[i], df.t2[i]) for i in range(1,size(df)[1])];
plot(df.detuning_log, df.detuned_vy_end, group=df.t1t2, seriestype=:scatter,
	xlabel="Log10(Detuning Frequency (Hz))", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);

# Either group in for loop, or group in plot
for detuning_log_df in groupby(df, [:detuning_log])
	plot(detuning_log_df.t2, detuning_log_df.detuned_vy_end, title="$(10^detuning_log_df.detuning_log[1])", group=detuning_log_df.t1,
		xlabel="T2", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);
	# plot(detuning_log_df.t2, detuning_log_df.ramsey_vy_end, group=detuning_log_df.t1,
	# 	xlabel="T2", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);
end

for t1_df in groupby(df, [:t1])
	plot(t1_df.t2, t1_df.detuned_vy_end, title="$(t1_df.t1[1])", group=t1_df.detuning_log,
		xlabel="T2", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);
	# plot(detuning_log_df.t2, detuning_log_df.ramsey_vy_end, group=detuning_log_df.t1,
	# 	xlabel="T2", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);
end

for t2_df in groupby(df, [:t2])
	plot(t2_df.t1, t2_df.detuned_vy_end, title="$(t2_df.t2[1])", group=t2_df.detuning_log,
		xlabel="T1", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);
	# plot(detuning_log_df.t2, detuning_log_df.ramsey_vy_end, group=detuning_log_df.t1,
	# 	xlabel="T2", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);
end
# Count unique values in column
# combine(groupby(df, :t1t2), nrow=>:count)

# List each column
# show(names(df))

# Check element types
# eltype.(eachcol(df))

# Number of rows
# size(df)

# Delete at row
# delete!(df, []);

# Recast column
# df.t1 = parse.(Float64, df.t1);

