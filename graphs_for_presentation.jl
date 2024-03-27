using DifferentialEquations;
using LinearAlgebra;

# The setup for the experiment is:
# 1. Z-axis dephasing
# 2. Pure relaxation, i.e. relaxation to ground state
# The initial state is the bloch vector (1/sqrt(2),0,1/sqrt(2))
# The output is a Hamiltonian, with units in angular frequency. Rabi frequency
# is linear frequency, so divide by 2pi to get in terms of Rabi frequency.

# Experiment setup
h = 6.62607015 * 1e-34; # Planck's constant - Joule per Hertz
kb = 1.380649 * 1e-23; # Boltzmann's constant - Joule per Kelvin

# 1/t2 = 1/(2t1) + 1/t_phi.
# t_phi - corresponds to dephasing. Equal to 1/gamma
# t1 - corresponds to thermal relaxation.
t1 = 50 * 1e-6;
t2 = 60 * 1e-6;
dephasing_gamma = 0.5*((1/t2)-(1/(2*t1)));
thermal_gamma = 1 / t1;

@enum SimulationType begin
   ramsey_detuned = 1
   ideal = 2
   detuned = 3
end;

# Arguments
if length(ARGS) < 0
	print("Needs detuning frequency as argument\n");
	exit();
end
detuning_freq = 1/(2*t2)#parse(Float64, ARGS[1]);;

# Sampling rate is 2.4 giga samples per sec
sampling_rate = 2.4 * 1e9;
qbit_freq = 3877496000;
bath_temp = 50 * 10e-3;

# Simulation code
F4 = Matrix{ComplexF64}[];
push!(F4, [1 0;0 1]);
push!(F4, [0 1;1 0]);
push!(F4, [0 -1im;1im 0]);
push!(F4, [1 0;0 -1]);

matrix_to_coeff(matrix) = [real(tr(x' * matrix)) for x in F4];
coeff_to_matrix(coeff) = 0.5 * reduce(+, [x*y for (x,y) in zip(cat([1],coeff,dims=1), F4)])

commutator(x, y) = x*y - y*x;

function get_relaxation_coeff(temp, qbit_freq)
	# Returns the rate associated with relaxation.
	#
	# Assume qubit energy is -(1/2)ωσ. Then ground state population in Gibbs state is
	# 1/(1+exp(-βH)). The probability p of the excitation channel is simply the population
	# of excited state at equilibrium.
	βH = (h * qbit_freq) / (kb * temp);
	return 1/(1+exp(-βH))
end

function dissipator(v)
	dephasing_dissipator = -2 * dephasing_gamma * [v[1], v[2], 0];
	relaxation_coefficient = 1 # get_relaxation_coeff(bath_temp, qbit_freq);
	thermal_dissipator = thermal_gamma * (-relaxation_coefficient * [v[1]/2, v[2]/2, v[3]-1] - (1-relaxation_coefficient) * [v[1]/2, v[2]/2, v[3]+1]);
	return dephasing_dissipator + thermal_dissipator;
end

function target(v)
	return v[1]^2 + v[2]^2;
end

function grad(v)
	return 2 * [v[1], v[2], 0];
end

function get_hamiltonian(v,p,t)
	if p[1] == ramsey_detuned::SimulationType
		return [0,0,(detuning_freq/2) * 2 * pi];
	elseif p[1] == ideal::SimulationType
		dephasing_hamiltonian = -dephasing_gamma / v[3] * [v[2], -v[1], 0];
		thermal_hamiltonian = -thermal_gamma / (4 * v[3]) * [v[2], -v[1], 0];
		return dephasing_hamiltonian + thermal_hamiltonian;
	elseif p[1] == detuned::SimulationType
		buff_v = p[2](t);
		dephasing_hamiltonian = -dephasing_gamma / buff_v[3] * [buff_v[2], -buff_v[1], 0];
		thermal_hamiltonian = -thermal_gamma / (4 * buff_v[3]) * [buff_v[2], -buff_v[1], 0];
		return dephasing_hamiltonian + thermal_hamiltonian + [0,0,detuning_freq/2 * 2 * pi];
	end
	throw(Error("Undefined control sequence in get_hamiltonian: unexpected value of SimulationType"))
end

# Parameters (tuple):
# 1. First element - is_ramsey_setup
# 2. Second element - ideal solution
function lindblad(v, p, t)
	hamiltonian = get_hamiltonian(v,p,t);
	if any(isnan, hamiltonian)
		return [Inf, Inf, Inf]
	end
	return 2 * cross(hamiltonian, v) + dissipator(v)
end

function get_hamiltonian_matrix(u, t, int)
	h_vec = get_hamiltonian(u, (ideal::SimulationType, nothing), 0)
	return h_vec[1] * F4[2] + h_vec[2] * F4[3] + h_vec[3] * F4[4]
end

function save_target(u, t, int)
	return target(u)
end

function integral(solution)
	prev_time = 0
	rectangle_estimates = [solution.u[i] * (solution.t[i] - solution.t[i-1]) for i in range(2,length(solution))]
	return sum(rectangle_estimates)
end

is_ramsey_setup = false;
simulation_type = ideal::SimulationType;
v = [1/sqrt(2),0,1/sqrt(2)]
# v = [0.41,0,0.7]
if !is_ramsey_setup
	x = 0.5*sqrt(2 * thermal_gamma / (4*dephasing_gamma + thermal_gamma))-0.01
	# x = 0.66;
	# x = 0.3;
	v = [x,0,sqrt(1-x^2)];
	# v=v/2;
else
	simulation_type = ramsey_detuned::SimulationType;
	v = [1,0,0]
end
abstol, reltol = 1e-8,1e-6;
tend = 2*t2;

problem = ODEProblem(lindblad, v, (0.0, tend), (simulation_type, nothing));
ideal_solution = solve(problem, alg_hints=[:stiff], abstol=abstol, reltol=reltol);

# using JLD2;
# @save "ideal_solution.jld2" ideal_solution
# print("Initial state: ", v, "\n");
# print("T1, T2: ", t1, ", ", t2, "\n");
# print("Simulation time for ideal solution: ", ideal_solution.t[end], "\n");
# print("Is ramsey? ", is_ramsey_setup, "\n");
# z = ideal_solution.u[end][3]
using Plots;
plotly();
# plot(ideal_solution, size=(1500,800), show=true, title="T1=$(t1*10^6)us, T2=$(t2*10^6)us, Ideal", label=["vx ideal" "vy ideal" "vz ideal"]);
# plot(ideal_solution.t, [target(u) for u in ideal_solution.u], show=true, ylim=(0,1), label="target ideal")
# plot(ideal_solution.t, [[get_hamiltonian(u)[1] for u in ideal_solution.u], [get_hamiltonian(u)[2] for u in ideal_solution.u], [get_hamiltonian(u)[3] for u in ideal_solution.u]], show=true, label=["hx ideal" "hy ideal" "hz ideal"])

# problem = ODEProblem(lindblad, v, (0.0, tend), (detuned::SimulationType, ideal_solution));
# detuned_solution = solve(problem, alg_hints=[:stiff], abstol=abstol, reltol=reltol);
# @save "detuned_solution.jld2" detuned_solution
# detuned_hamiltonian = [[get_hamiltonian(u)[1] for u in ideal_solution.u], [get_hamiltonian(u)[2] for u in ideal_solution.u], [get_hamiltonian(u)[3] for u in ideal_solution.u]]

# plot(detuned_solution, size=(1500,800), show=true, title="T1=$(t1*10^6)us, T2=$(t2*10^6)us, Detuning=$(detuning_freq)Hz (CM)", label=["vx" "vy" "vz"]);
# plot(detuned_solution.t, [target(u) for u in detuned_solution.u], show=true, ylim=(0,1), label="target")
# print("Ideal solution - area under v_y: ", integral(ideal_solution), "\n")
# print("Detuned solution - area under v_y: ", integral(detuned_solution), "\n")
# detuned_vy=[abs(detuned_solution.u[i][2]) for i in range(1,length(detuned_solution))];
# print("v_y max in detuned solution: ", maximum(vy), "\n");
# print(detuned_vy[end], "\t", ideal_solution.u[end][3], "\n");


# is_ramsey_setup = true
problem = ODEProblem(lindblad, [1,0,0], (0.0, 2*tend), (ramsey_detuned::SimulationType, ideal_solution));
ramsey_solution = solve(problem, alg_hints=[:stiff], abstol=abstol, reltol=reltol);
# ramsey_vy=[abs(ramsey_solution.u[i][2]) for i in range(1,length(ramsey_solution))];
plot(ramsey_solution, title="High detuning",  linewidth=5, size=(1500,800), show=true, xlabel="Time", ylim=(-0.5,1), label=["vx" "vy" "vz"], xticks=([4*t2*x/16 for x in 1:1:16], ["$(x/8) *T2" for x in 1:1:16]));
# print("$(2/detuning_freq),$(detuning_freq),$(maximum(ramsey_vy)),$(ramsey_vy[end])\n");[latexstring("\\frac{$(x)}{8\\delta}") for x in 1:1:16]




