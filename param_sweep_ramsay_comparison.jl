using DifferentialEquations;
using LinearAlgebra;

# This program runs a parameter sweep for coherence preservation enhanced quantum sensing.
# The setup for the experiment is:
# 1. Z-axis dephasing
# 2. Pure relaxation, i.e. relaxation to ground state
# The initial state is the bloch vector (1/sqrt(2),0,1/sqrt(2))
# The output is a Hamiltonian, with units in angular frequency. Rabi frequency
# is linear frequency, so divide by 2pi to get in terms of Rabi frequency.

# Experiment setup
h = 6.62607015 * 1e-34; # Planck's constant - Joule per Hertz
kb = 1.380649 * 1e-23; # Boltzmann's constant - Joule per Kelvin

# Sampling rate is 2.4 giga samples per sec
sampling_rate = 2.4 * 1e9;
qbit_freq = 3877496000;
bath_temp = 50 * 10e-3;

SIMULATION_TYPE = 1;
DETUNING_FREQ = 2;
DEPHASING_GAMMA = 3;
THERMAL_GAMMA = 4;
IDEAL_TRAJECTORY = 5;

@enum SimulationType begin
	# For detuning experiments, the static field is assumed to be (detuning_freq)/2 * σz
	# in linear frequency units. With no dissipation and external control, this is
	# expected to give oscillations with time period (1/detuning_freq).
	ramsey_interferometry = 1
	ideal_tracking_control = 2
	detuned_tracking_control = 3
end;

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

function dissipator(v,p,t)
	dephasing_dissipator = -2 * p[DEPHASING_GAMMA] * [v[1], v[2], 0];
	assume_pure_relaxation = true;
	relaxation_coefficient = assume_pure_relaxation ? 1 : get_relaxation_coeff(bath_temp, qbit_freq);
	thermal_dissipator = p[THERMAL_GAMMA] * (-relaxation_coefficient * [v[1]/2, v[2]/2, v[3]-1] - (1-relaxation_coefficient) * [v[1]/2, v[2]/2, v[3]+1]);
	return dephasing_dissipator + thermal_dissipator;
end

function target(v)
	# Coherence magnitude = vx^2+vy^2;
	return v[1]^2 + v[2]^2;
end

function grad(v)
	return 2 * [v[1], v[2], 0];
end

function get_hamiltonian(v,p,t)
	# For conversion from Hamiltonian in Pauli basis to vectors, we assume
	# H = hx σx + hy σy + hz σz. This gives vec(h) = (hx, hy, hz). The resulting
	# differential equation is dv/dt = 2 h x v. Note that the detuning Hamiltonian
	# is detuning_freq/2 in linear frequency units.
	detuning_hamiltonian = [0,0,p[DETUNING_FREQ]/2 * 2 * pi];
	# dephasing_hamiltonian(x) = -dephasing_gamma / x[3] * [x[2], -x[1], 0];
	# thermal_hamiltonian(x) = -thermal_gamma / (4 * x[3]) * [x[2], -x[1], 0];
	if p[SIMULATION_TYPE] == ramsey_interferometry::SimulationType
		return detuning_hamiltonian;
	elseif p[SIMULATION_TYPE] == ideal_tracking_control::SimulationType
		dephasing_hamiltonian = -p[DEPHASING_GAMMA] / v[3] * [v[2], -v[1], 0];
		thermal_hamiltonian = -p[THERMAL_GAMMA] / (4 * v[3]) * [v[2], -v[1], 0];
		return dephasing_hamiltonian + thermal_hamiltonian;
	elseif p[SIMULATION_TYPE] == detuned_tracking_control::SimulationType
		# Use the ideal case hamiltonian at time t to apply here. The
		# ideal case hamiltonian is a function of the state in ideal case at time t.
		if t > p[IDEAL_TRAJECTORY].t[end]*0.9
			return detuning_hamiltonian;
		end
		buff_v = p[IDEAL_TRAJECTORY](t);
		dephasing_hamiltonian = -p[DEPHASING_GAMMA] / buff_v[3] * [buff_v[2], -buff_v[1], 0];
		thermal_hamiltonian = -p[THERMAL_GAMMA] / (4 * buff_v[3]) * [buff_v[2], -buff_v[1], 0];
		return dephasing_hamiltonian + thermal_hamiltonian + detuning_hamiltonian;
	end
	throw(Error("Undefined control sequence in get_hamiltonian: unexpected value of SimulationType"))
end

# Parameters p (tuple):
# 1. First element - is_ramsey_setup
# 2. Second element - ideal solution
function lindblad(v, p, t)
	hamiltonian = get_hamiltonian(v,p,t);
	if any(isnan, hamiltonian)
		return [Inf, Inf, Inf]
	end
	return 2 * cross(hamiltonian, v) + dissipator(v,p,t)
end

# Parameters
trial_run = true;
fine_sampling = false;
extra_precision = true;

param_space = Dict(
		"t1" => range(10,200,step=10), # remember to convert to microseconds
		"t2" => range(10,400,step=10), # remember to convert to microseconds
		"detune_ratio" => append!(collect(range(1/10,1/2,length=30)), 1 ./ collect(range(10,100,length=30)))
	);
if trial_run
	param_space = Dict(
		"t1" => range(50,50,step=10), # remember to convert to microseconds
		"t2" => range(50,60,step=10), # remember to convert to microseconds
		"detune_ratio" => [0,0.2],#append!(collect(range(1/10,1/2,length=30)), 1 ./ collect(range(10,100,length=30)))#range(1/20,1/2,length=5)
	);
end

sampling = [];
if fine_sampling
	sampling = 1/sampling_rate;
end

abstol, reltol = 1e-6, 1e-3;
if extra_precision
	abstol, reltol = 1e-8,1e-10;
end

using NLopt, ForwardDiff;

# Find the max vy using the solution to the differential equation.
function max_vy_objective(x::Vector, grad::Vector, diffEqSolution)
	if length(grad) > 0
    # use ForwardDiff for the gradient for vy
    grad[1] = ForwardDiff.derivative((x)->diffEqSolution(first(x), idxs=2), x)
  end
  diffEqSolution(first(x), idxs=2)	
end

# Find the time to be less than the total simulation time for the differential equation solution.
function t_constraint(x::Vector, grad::Vector, diffEqSolution)
    if length(grad) > 0
        grad[1] = 1
    end
    x[1]-diffEqSolution.t[end]
end

function get_max_vy(detuned_solution, ideal_solution, ramsey_solution)
	opt = NLopt.Opt(:GN_ORIG_DIRECT_L, 1);
	NLopt.lower_bounds!(opt, [0.0]);
	NLopt.upper_bounds!(opt, [1.0]);
	NLopt.max_objective!(opt, (x,g)->max_vy_objective(x,g,detuned_solution));
	NLopt.inequality_constraint!(opt, (x,g) -> t_constraint(x,g,detuned_solution), 1e-8);
	(max_vy,argmax_t,ret) = NLopt.optimize(opt,[0.5]);

	# Now find maximum vy for t < Breakdown
	opt = NLopt.Opt(:GN_ORIG_DIRECT_L, 1);
	NLopt.lower_bounds!(opt, [0.0]);
	NLopt.upper_bounds!(opt, [1.0]);
	# NLopt.xtol_rel!(opt); # ,1e-8
	NLopt.max_objective!(opt, (x,g)->max_vy_objective(x,g,detuned_solution));
	NLopt.inequality_constraint!(opt, (x,g) -> t_constraint(x,g,ideal_solution), 1e-8);
	(max_vy_b,argmax_t_b,ret) = NLopt.optimize(opt,[0.5]);

	# Now find maximum vy for ramsey interferometry
	opt = NLopt.Opt(:GN_ORIG_DIRECT_L, 1);
	NLopt.lower_bounds!(opt, [0.0]);
	NLopt.upper_bounds!(opt, [1.0]);
	# NLopt.xtol_rel!(opt); # ,1e-8
	NLopt.max_objective!(opt, (x,g)->max_vy_objective(x,g,ramsey_solution));
	NLopt.inequality_constraint!(opt, (x,g) -> t_constraint(x,g,ramsey_solution), 1e-8);
	(max_vy_r,argmax_t_r,ret) = NLopt.optimize(opt,[0.5]);

	return (max_vy, argmax_t, max_vy_b, argmax_t_b, max_vy_r, argmax_t_r);
end

header = "";
header = header * "t1"; # T1
header = header * ",t2"; # T2
header = header * ",detune_ratio"; # Detuning freq * T2
header = header * ",vx"; # Starting state
header = header * ",ideal_t"; # Breakdown time or max simulation time
header = header * ",max_vy"; # Max vy - coherence magnitude
header = header * ",argmax_t"; # Time for max vy
header = header * ",max_vy_b"; # Max vy - before Breakdown
header = header * ",argmax_t_b"; # Time for max vy - before Breakdown
header = header * ",max_vy_r"; # Max vy - ramsey
header = header * ",vy_end"; # vy end - coherence magnitude
header = header * ",vy_end_r"; # vy end - ramsey
header = header * "\n";
print(header);
@time for t1us in param_space["t1"]
	t1 = t1us * 10^-6;
	for t2us in param_space["t2"]
		t2 = t2us * 10^-6;
		dephasing_gamma = 0.5*((1/t2)-(1/(2*t1)));
		if dephasing_gamma < 0
			continue
		end
		thermal_gamma = 1 / t1;
		# vx_maximum for stable regions is supposed to be
		# 0.5*sqrt(2 * thermal_gamma / (4*dephasing_gamma + thermal_gamma)),
		# which is the same as 0.5*sqrt(t2:t1)
		vx_maximum = 0.98; # minimum((1,0.5*sqrt(t2us/t1us)))-0.0001;
		vx_minimum = 0.97#0.5*sqrt(t2/t1)+0.01; #maximum((0.1,vx_maximum-floor(10*vx_maximum)/10));
		for detune_ratio in param_space["detune_ratio"]
			detuning_freq = detune_ratio/t2;
			tend = 10*t2;
			# Ramsey interference setup
			simulation_params = (
				ramsey_interferometry::SimulationType, # simulation type
				detuning_freq,
				dephasing_gamma,
				thermal_gamma,
				nothing, # ideal solution/trajectory data
			);
			problem = ODEProblem(lindblad, [1,0,0], (0.0, tend), simulation_params);
			# For some reason, there's a problem with sampling rate
			ramsey_solution = solve(problem, alg_hints=[:stiff], saveat=sampling, abstol=abstol, reltol=reltol);
			ramsey_vy=[abs(ramsey_solution.u[i][2]) for i in range(1,length(ramsey_solution))];

			for buff_vx in range(vx_maximum,vx_minimum,step=-0.01)
				vx = buff_vx;
				vz = sqrt(1-vx^2);
				# Ideal coherence magnitude preserving control
				v = [vx,0,vz];
				simulation_params = (
					ideal_tracking_control::SimulationType, # simulation type
					detuning_freq,
					dephasing_gamma,
					thermal_gamma,
					nothing, # ideal solution/trajectory data
				);
				problem = ODEProblem(lindblad, v, (0.0, tend), simulation_params);
				ideal_solution = solve(problem, alg_hints=[:stiff], saveat=sampling, abstol=abstol, reltol=reltol);

				# Detuned coherence magnitude preserving control
				simulation_params = (
					detuned_tracking_control::SimulationType, # simulation type
					detuning_freq,
					dephasing_gamma,
					thermal_gamma,
					ideal_solution, # ideal solution/trajectory data
				);
				problem = ODEProblem(lindblad, v, (0.0, ideal_solution.t[end]), simulation_params);
				detuned_solution = solve(problem, alg_hints=[:stiff], saveat=sampling, abstol=abstol, reltol=reltol);
				
				(max_vy, argmax_t, max_vy_b, argmax_t_b, max_vy_r, argmax_t_r) = get_max_vy(detuned_solution, ideal_solution, ramsey_solution);
				
				# alpha = vx^2*t1/t2-1/4;
				# predicted = -t1*( 0.5*log((0.25+alpha)/((vz-0.5)^2+alpha)) + 1/(2*sqrt(alpha)) * (atan(-0.5/sqrt(alpha))-atan((vz-0.5)/sqrt(alpha))));
				output = "";
				output = output * "$(t1)"; # T1
				output = output * ",$(t2)"; # T2
				output = output * ",$(detune_ratio)"; # Detuning freq * T2
				output = output * ",$(vx)"; # Starting state
				output = output * ",$(ideal_solution.t[end])"; # Breakdown time or max simulation time
				output = output * ",$(max_vy)"; # Max vy - coherence magnitude
				output = output * ",$(argmax_t)"; # Time for max vy
				output = output * ",$(max_vy_b)"; # Max vy - before Breakdown
				output = output * ",$(argmax_t_b)"; # Time for max vy - before Breakdown
				output = output * ",$(max_vy_r)"; # Max vy - ramsey
				output = output * ",$(detuned_vy[end])"; # vy end - coherence magnitude
				output = output * ",$(ramsey_vy[end])"; # vy end - ramsey
				print("$(output)\n");
			end
		end
	end
end













