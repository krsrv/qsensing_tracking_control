using CSV;
using DataFrames;
using Plots;
plotly();

# Only 60 seconds to create output for vy
df = CSV.read("detuning_to_t2_fixed_v.csv", DataFrame; header=1);
# rename!(df,[:t1,:t2,:detuning,:vx,:sim_z,:sim_h,:sim_dis,:cm_vy_max,:ramsey_vy_max,:cm_vy_max_t,:ramsey_vy_max_t,:cm_vy_end,:ramsey_vy_end,:sim_ratio]);

# df.detuning_log = log10.(df.detuning);
df.t1 = ceil.(df.t1*10^6);
df.t2 = ceil.(df.t2*10^6);

# filter([:ramsey_vy_max,:cm_vy_end]=>(x,y)->x<y,df)
function calculate_predicted_vy_end(t1,t2,detuning_t2,vx)
	dephasing_gamma = 0.5*10^6*((1/t2)-(1/(2*t1)));
	thermal_gamma = 10^6 / t1;
	detuning = detuning_t2/t2*10^6;
	# alpha = dephasing_gamma + thermal_gamma/4;
	# vz_inf = 0.5 + 0.5*sqrt(1-8*vx^2*alpha/thermal_gamma); # Predicted stable z for ideal solution
	# H = alpha * vx / vz_inf; # Predicted stable H for ideal solution
	# angular_detuning = detuning * pi;
	# vy_end = H * angular_detuning / (alpha^2 + angular_detuning^2 + 2*H*alpha/thermal_gamma);
	# vx_end = alpha*vy_end/angular_detuning;
	# vz_end = 1-2*H*vx_end/thermal_gamma;
	# return [vz_inf, H, vy_end, vx_end^2+vy_end^2+vz_end^2<=1];
	vy_end = detuning*sqrt()
	return [vy_end];
end
transform!(df, [:t1,:t2,:detuning_t2,:vx] => ByRow(calculate_predicted_vy_end) => [:predicted_vz,:predicted_H,:predicted_vy_end,:valid]);
df[:,[:t1,:t2,:detuning_t2,:vx,:predicted_vy_end,:cm_vy_end]]


# Plot vy as a function of detuning * t2
for buff_df in groupby(df, [:t1])
	if buff_df.t1[1] != 60
		continue
	end
	sort!(buff_df,[:detuning_t2]);
	labels=["T2 - $(x.t2) us" for x in eachrow(unique(buff_df,:t2))];
	
	graph = plot(buff_df.detuning_t2, log10.(buff_df.cm_vy_max),
		title="T1 = $(buff_df.t1[1]) us", group=buff_df.t2,
		xlabel="Detuning × T2", ylabel="Log10(vy)",
		labels=reshape([x*" (CM max)" for x in labels], 1, length(labels)),
		size=(1500,800));
	# graph = plot(1 ./ buff_df.detuning_t2, log10.(buff_df.cm_vy_max), title="T1 = $(buff_df.t1[1]) us", group=buff_df.t2,
	# 	xlabel="Gamma/Detuning", ylabel="Log(vy)", labels=reshape([x*" (CM max)" for x in labels], 1, length(labels)), size=(1500,800));
	
	plot_only_one_ramsey = true
	if plot_only_one_ramsey
		another_buff_df = filter([:t2]=> (x) -> buff_df.t2[1]==x, buff_df)
		graph = plot!(another_buff_df.detuning_t2, log10.(another_buff_df.ramsey_vy_max),
			title="T1 = $(another_buff_df.t1[1]) us",
			xlabel="Detuning × T2", ylabel="Log10(vy)",
			labels="Ramsey",
			size=(1500,800), linestyle=:dash, linecolor=:black, markershape=:circle);
	else
		graph = plot!(buff_df.detuning_t2, log10.(buff_df.ramsey_vy_max),
			title="T1 = $(buff_df.t1[1]) us", group=buff_df.t2,
			xlabel="Detuning × T2", ylabel="Log10(vy)",
			labels=reshape([x*" (Ramsey)" for x in labels],1, length(labels)),
			size=(1500,800), linestyle=:dash);
	end
	display(graph);
end

# Plot # experiment runs as a function of detuning*t2
for buff_df in groupby(df, [:t1])
	if buff_df.t1[1] != 60
		continue
	end
	std_ratio = 0.1;
	new_df = filter([:ramsey_vy_max,:cm_vy_max]=> (x,y) -> x<y, buff_df)
	sort!(new_df,[:detuning_t2]);
	transform!(new_df, [:ramsey_vy_max,:cm_vy_max] => ByRow((x,y)->[(1-x^2)/((std_ratio*x/2)^2),(1-y^2)/((std_ratio*y/2)^2)]) => [:nruns_ramsey,:nruns_cm_max]);
	labels=["T2 - $(x.t2) us" for x in eachrow(unique(new_df,:t2))];
	
	graph = plot(new_df.detuning_t2, log10.(new_df.nruns_cm_max),
		title="T1 = $(new_df.t1[1]) us, Relative stddev = $(std_ratio)", group=new_df.t2,
		xlabel="Detuning × T2", ylabel="Log10(# Experiment runs)",
		labels=reshape([x*" (CM max)" for x in labels], 1, length(labels)),
		size=(1500,800));

	plot_only_one_ramsey = true
	if plot_only_one_ramsey
		another_buff_df = filter([:t2]=> (x) -> maximum(new_df.t2)==x, new_df)
		graph = plot!(another_buff_df.detuning_t2, log10.(another_buff_df.nruns_ramsey),
			title="T1 = $(another_buff_df.t1[1]) us, Relative stddev = $(std_ratio)",
			xlabel="Detuning × T2", ylabel="Log10(# Experiment runs)",
			labels="Ramsey",
			size=(1500,800), linestyle=:dash, linecolor=:black, markershape=:circle);
	else
		graph = plot!(new_df.detuning_t2, new_df.nruns_ramsey,
			title="T1 = $(new_df.t1[1]) us, Target stddev = $(std)", group=new_df.t2,
			xlabel="Detuning × T2", ylabel="# Experiment runs",
			labels=reshape([x*" (Ramsey)" for x in labels], 1, length(labels)),
			size=(1500,800), linestyle=:dash);
	end
	
	graph = plot!(new_df.detuning_t2, [log10(15000) for x in eachrow(new_df)], labels=reshape(["15000"],1,1), linestyle=:dashdot);
	display(graph);
end


# Plot vy as a function of x, fixing T1 and T2 - use "breakdown.csv" files which has data for varying vx
for buff_df in groupby(df, [:t1,:t2])
	if buff_df.t1[1] != 60 && buff_df.t2[1]!=90
		continue
	end
	sort!(buff_df,[:vx]);
	labels=["Detuning*T2 - $(x.detuning_t2)" for x in eachrow(unique(buff_df,:detuning_t2))];
	
	graph = plot(buff_df.vx, log10.(buff_df.cm_vy_max),
		title="T1 = $(buff_df.t1[1]) us, T2 = $(buff_df.t2[1])", group=buff_df.detuning_t2,
		xlabel="vx", ylabel="Log10(vy)",
		labels=reshape([x*" (Ramsey)" for x in labels], 1, length(labels)),legend=:bottomright,
		size=(1500,800));
	graph = vline!([0.5*sqrt(buff_df.t2[1]/buff_df.t1[1])]);
	for new_df in groupby(buff_df, [:detuning_t2])
		sort!(new_df, [:cm_vy_max]);
		print("$(new_df.t1[1]), $(new_df.t2[1]), $(new_df.sim_ratio[end]/new_df.t2[end]))\n");
	end
	display(graph);
end


# Plot vy as a function of x, fixing T1 and Detuning*T2 - use "breakdown.csv" files which has data for varying vx
for buff_df in groupby(df, [:t1,:detuning_t2])
	if buff_df.t1[1] != 60 && buff_df.t2[1]!=90
		continue
	end
	sort!(buff_df,[:vx]);
	labels=["T2 - $(x.t2)" for x in eachrow(unique(buff_df,:t2))];
	
	graph = plot(buff_df.vx, log10.(buff_df.cm_vy_max),
		title="T1 = $(buff_df.t1[1]) us, Detuning*T2 = $(buff_df.detuning_t2[1])", group=buff_df.t2,
		xlabel="vx", ylabel="Log10(vy)",
		labels=reshape([x*" (Ramsey)" for x in labels], 1, length(labels)),legend=:bottomright,
		size=(1500,800));
	for x in eachrow(unique(buff_df,:t2))
		graph = plot!([0.5*sqrt(x.t2[1]/x.t1[1])], seriestype=:vline, labels="vmax for ∞ - T2=$(x.t2[1])")
	end
	display(graph);
end


# Plot simulation time for max vy as a function of T2, fixing T1 and Detuning*T2 - use "breakdown.csv" files which has data for varying vx
for buff_df in groupby(df, [:t1])
	# https://stackoverflow.com/questions/65024962/select-rows-of-a-dataframe-containing-minimum-of-grouping-variable-in-julia
	new_df = combine(sdf -> sdf[argmax(sdf.cm_vy_max), :], groupby(buff_df, [:t2, :detuning_t2]))
	title="T1 = $(buff_df.t1[1]) us";
	labels=["Detuning*T2 - $(x.detuning_t2)" for x in eachrow(unique(new_df,:detuning_t2))];
	sort!(new_df, [:t2]);
	graph = plot(new_df.t2, new_df.sim_ratio, group=new_df.detuning_t2,
			title=title, xlabel="T2", ylabel="Simulation time/T2",
			size=(1500,800), labels=reshape(labels, 1, length(labels)))
	display(graph);
end


# Plot max vy as a function of simulation time for max vy, for a fixed T1, T2, Detuning*T2 - use "breakdown.csv" files which has data for varying vx
for buff_df in groupby(df, [:t1,:t2,:detuning_t2])
	if buff_df.t1[1] != 50 || buff_df.t2[1] != 60 || abs(buff_df.detuning_t2[1]-0.127) > 0.001
		continue
	end
	# https://stackoverflow.com/questions/65024962/select-rows-of-a-dataframe-containing-minimum-of-grouping-variable-in-julia
	title="T1 = $(buff_df.t1[1]) us, T2 = $(buff_df.t2[1]) us, Detuning*T2 = $(buff_df.detuning_t2[1])";
	sort!(buff_df, [:sim_ratio]);
	graph = plot(buff_df.sim_ratio, log10.(buff_df.cm_vy_max), group=buff_df.t2,
			title=title, xlabel="Simulation time/T2", ylabel="Log10(vy)",
			size=(1500,800), show=false)
	display(graph);
end

# Plot breakdown time for which max vy happens as a function of detuning*t2, for a fixed T1, T2, Detuning*T2 - use "breakdown.csv" files which has data for varying vx
for buff_df in groupby(df, [:t1])
	# https://stackoverflow.com/questions/65024962/select-rows-of-a-dataframe-containing-minimum-of-grouping-variable-in-julia
	title="T1 = $(buff_df.t1[1]) us";
	new_df = combine(sdf -> sdf[argmax(sdf.cm_vy_max), :], groupby(buff_df, [:t2, :detuning_t2]))
	labels=["T2 - $(x.t2)" for x in eachrow(unique(new_df,:t2))];
	sort!(new_df, [:detuning_t2]);
	new_df.detuning = 1e6*new_df.detuning_t2 ./ new_df.t2;
	graph = plot(new_df.detuning_t2, new_df.sim_ratio, group=new_df.t2,
			title=title, xlabel="Detuning*T2", ylabel="Simulation time/T2 corr to max vy",
			size=(1500,800), labels=reshape(labels, 1, length(labels)), marker=:circle)
	display(graph);
end

# Plot vy_max_time as a function of vx, for a fixed T1, T2, Detuning*T2 - use "breakdown.csv" files which has data for varying vx
for buff_df in groupby(df, [:t1,:t2])
	# https://stackoverflow.com/questions/65024962/select-rows-of-a-dataframe-containing-minimum-of-grouping-variable-in-julia
	title="T1 = $(buff_df.t1[1]) us, T2 = $(buff_df.t2[1]) us";
	new_df = filter([:detuning_t2]=>x->(abs(x-1/3)<0.03 || (abs(x-0.01)<0.001)), buff_df);
	labels=["Detuning*T2 - $(x.detuning_t2)" for x in eachrow(unique(new_df,:detuning_t2))];
	sort!(new_df, [:vx]);
	# new_df.detuning = 1e6*new_df.detuning_t2 ./ new_df.t2;
	graph = plot(new_df.vx, 1e6*new_df.cm_vy_time ./ new_df.t2, group=new_df.detuning_t2,
			title=title, xlabel="vx", ylabel="Peak vy time/T2",
			size=(1500,800), labels=reshape(labels, 1, length(labels)), marker=:circle, show=false)
	graph = plot!([(0.5*sqrt(new_df.t2[1] / new_df.t1[1]))], seriestype=:vline, labels="vx_max for ∞ breakdown");
	display(graph);
end

# Plot vx vs breakdown time - for 0 detuning ideally
for buff_df in groupby(df, [:t1,:t2])
	if buff_df.t1[1] != 60 && buff_df.t2[1]!=90
		continue
	end
	labels=["Detuning*T2 - $(x.detuning_t2)" for x in eachrow(unique(buff_df,:detuning_t2))];
	sort!(buff_df,[:vx]);
	graph = plot(buff_df.vx, log10.(buff_df.sim_ratio),
		title="T1 = $(buff_df.t1[1]) us, T2 = $(buff_df.t2[1])",
		xlabel="vx", ylabel="Simulation time",
		labels=reshape([x*" (Ramsey)" for x in labels], 1, length(labels)), legend=:bottomright,
		size=(1500,800));
	graph = plot!(buff_df.vx, [log10(30*buff_df.t2[1]/10^6) for i in eachrow(buff_df)], labels="30 * T2", size=(1500,800));
	display(graph);
end


##              ##
## Earlier Code ##
##              ##

function calculate_predicted_vy_end(t1,t2,detuning,vx)
	dephasing_gamma = 0.5*10^6*((1/t2)-(1/(2*t1)));
	thermal_gamma = 10^6 / t1;
	alpha = dephasing_gamma + thermal_gamma/4;
	vz_inf = 0.5 + 0.5*sqrt(1-8*vx^2*alpha/thermal_gamma); # Predicted stable z for ideal solution
	H = alpha * vx / vz_inf; # Predicted stable H for ideal solution
	angular_detuning = detuning * pi;
	vy_end = H * angular_detuning / (alpha^2 + angular_detuning^2 + 2*H*alpha/thermal_gamma);
	vx_end = alpha*vy_end/angular_detuning;
	vz_end = 1-2*H*vx_end/thermal_gamma;
	return [vz_inf, H, vy_end, vx_end^2+vy_end^2+vz_end^2<=1];
end
transform!(df, [:t1,:t2,:detuning,:vx] => ByRow(calculate_predicted_vy_end) => [:predicted_vz,:predicted_H,:predicted_vy_end,:valid]);
df[:,[:t1,:t2,:detuning,:vx,:predicted_vy_end,:cm_vy_end]]

# Either group in for loop, or group in plot
for buff_df in groupby(df, [:detuning_log])
	labels=["T1 - $(x.t1)us" for x in eachrow(unique(buff_df,:t1))];
	plot(buff_df.t2, buff_df.cm_vy_max, title="Detuning - $(10^buff_df.detuning_log[1]) Hz", group=buff_df.t1,
		xlabel="T2", ylabel="V_y at end", labels=reshape(labels, 1, length(labels)), ylimits=(0,1), size=(1500,800), show=true);
	# plot(buff_df.t2, buff_df.ramsey_vy_end, group=buff_df.t1,
	# 	xlabel="T2", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);
end

for buff_df in groupby(df, [:t1])
	labels=["T2 - $(x.t2) us" for x in eachrow(unique(buff_df,:t2))];
	# plot(buff_df.detuning_log, buff_df.cm_vy_end, title="Max v_y, (CM) T1 = $(buff_df.t1[1]) us", group=buff_df.t2,
	# 	xlabel="Log(Detuning)", ylabel="V_y at end", labels=reshape(labels, 1, length(labels)), ylimits=(0,1), size=(1500,800), show=true);
	plot(buff_df.detuning_log, buff_df.cm_vy_max_t, title="Time to max v_y, T1 = $(buff_df.t1[1]) us (CM)", group=buff_df.t2,
		xlabel="Log(Detuning)", ylabel="t", labels=reshape(labels, 1, length(labels)), size=(1500,800), show=true);
	# plot(buff_df.t2, buff_df.ramsey_vy_end, group=buff_df.t1,
	# 	xlabel="T2", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);
end

for buff_df in groupby(df, [:t2])
	if buff_df.t2[1] > 200
		continue
	end
	labels=["T1 - $(x.t1) us" for x in eachrow(unique(buff_df,:t1))];
	# plot(buff_df.detuning_log, buff_df.ramsey_vy_max, title="Max v_y, T2 = $(buff_df.t2[1]) us (CM)", group=buff_df.t1,
	# 	xlabel="Log(Detuning)", ylabel="V_y at end", labels=reshape(labels, 1, length(labels)), ylimits=(0,1), size=(1500,800), show=true);
	plot(buff_df.detuning_log, buff_df.ramsey_vy_max_t, title="Time to max v_y, T2 = $(buff_df.t2[1]) us (Ramsey)", group=buff_df.t1,
		xlabel="Log(Detuning)", ylabel="t", labels=reshape(labels, 1, length(labels)), size=(1500,800), show=true);
	# plot(buff_df.t2, buff_df.ramsey_vy_end, group=buff_df.t1,
	# 	xlabel="T2", ylabel="V_y at end", ylimits=(0,1), size=(1500,800), show=true);
end

function quadratic_fit_sample(df)
	fit_df = filter([:t1,:t2,:detuning_log] => (x,y,z)-> (x==70 && y==70), df);
	plot(fit_df.detuning, fit_df.cm_vy_end, seriestype=:scatter, title="V_y end vs detuning (T1=T2=70us)", xlabel="Detuning (Hz)", ylabel="End v_y", label="Simulation", size=(1500,800));
	detuning_list = range(minimum(fit_df.detuning),maximum(fit_df.detuning));
	predictions = [calculate_predicted_vy_max(70,70,x,fit_df.vx[1]) for x in detuning_list];
	plot!(detuning_list, predictions, label="Prediction", size=(1500,800), show=true);
end

function filter_t(df,t1=nothing,t2=nothing)
	if t1 == nothing
		return filter([:t1,:t2]=>(x,y)->y==t2, df)
	elseif t2 == nothing
		return filter([:t1,:t2]=>(x,y)->x==t1, df)
	else
		return filter([:t1,:t2]=>(x,y)->x==t1 && y==t2, df)
	end
end

## Simulation choice
# vx = vxmax/2 to increase accuracy of last point in theoretic framework
#
#
## Plotting
# labels=["T1 - $(x.t1) us" for x in eachrow(unique(buff_df,:t1))];
# plot(x, y, title="text $(variable) us", ylims=(0,20000),
#	xlabel="X axis", ylabel="Y axis", labels=reshape([x*" (CM max)" for x in labels], 1, length(labels)),
#	linestyle=:dash, seriestype=:scatter, size=(1500,800), show=true);

## DF Manipulation

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
#
# Filter df
# filter([:ramsey_vy_max,:cm_vy_end]=>(x,y)->x<y,df)
#
# Transform/create new column
# transform!(df, [:t1,:t2,:detuning,:vx] => ByRow(calculate_predicted_vy_end) => [:predicted_vz,:predicted_H,:predicted_vy_end,:valid]);
