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
# t1 = nothing;
# t2 = nothing;
# dephasing_gamma = nothing;
# thermal_gamma = nothing;

# Arguments
# detuning_freq = nothing;

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
   ramsey_detuned = 1
   ideal = 2
   detuned = 3
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

function dissipator(v,p)
	dephasing_dissipator = -2 * p[DEPHASING_GAMMA] * [v[1], v[2], 0];
	relaxation_coefficient = 1 # get_relaxation_coeff(bath_temp, qbit_freq);
	thermal_dissipator = p[THERMAL_GAMMA] * (-relaxation_coefficient * [v[1]/2, v[2]/2, v[3]-1] - (1-relaxation_coefficient) * [v[1]/2, v[2]/2, v[3]+1]);
	return dephasing_dissipator + thermal_dissipator;
end

function target(v)
	return v[1]^2 + v[2]^2;
end

function grad(v)
	return 2 * [v[1], v[2], 0];
end

function get_hamiltonian(v,p,t)
	if p[SIMULATION_TYPE] == ramsey_detuned::SimulationType
		return [0,0,p[DETUNING_FREQ]/2 * 2 * pi];
	elseif p[SIMULATION_TYPE] == ideal::SimulationType
		dephasing_hamiltonian = -p[DEPHASING_GAMMA] / v[3] * [v[2], -v[1], 0];
		thermal_hamiltonian = -p[THERMAL_GAMMA] / (4 * v[3]) * [v[2], -v[1], 0];
		return dephasing_hamiltonian + thermal_hamiltonian
	elseif p[SIMULATION_TYPE] == detuned::SimulationType
		if p[IDEAL_TRAJECTORY] == nothing
			throw(Error("Ideal solution missing for detuned simulation"));
		end
		buff_p = (
			ideal::SimulationType, # simulation type
			0,
			p[DEPHASING_GAMMA],
			p[THERMAL_GAMMA],
			nothing, # ideal solution/trajectory data
		);
		return get_hamiltonian(p[IDEAL_TRAJECTORY](t), buff_p, t) + [0,0,p[DETUNING_FREQ]/2 * 2 * pi];
	else
		throw(Error("Undefined control flow - Simulation type does not match one of ramsey, ideal or detuned"));
	end
end

function lindblad(v, p, t)
	hamiltonian = get_hamiltonian(v,p,t);
	if any(isnan, hamiltonian)
		return [Inf, Inf, Inf]
	end
	return 2 * cross(hamiltonian, v) + dissipator(v, p)
end

# Parameters
trial_run = false;
fine_sampling = false;
extra_precision = true;

param_space = Dict(
		"t1" => range(10,200,step=10), # remember to convert to microseconds
		"t2" => range(10,400,step=10), # remember to convert to microseconds
		"detune_index" => range(1,4,step=0.3)
	);
if trial_run
	param_space = Dict(
		"t1" => range(50,60,step=10), # remember to convert to microseconds
		"t2" => range(50,50,step=10), # remember to convert to microseconds
		"detune_index" => range(1,4,step=0.3)
	);
end

sampling = [];
if fine_sampling
	sampling = 1/sampling_rate;
end

abstol, reltol = 1e-6, 1e-3;
if extra_precision
	abstol, reltol = 1e-8,1e-6;
end

@time for t1us in param_space["t1"]
	t1 = t1us * 10^-6;
	for t2us in param_space["t2"]
		t2 = t2us * 10^-6;
		dephasing_gamma = ((1/t2)-(1/(2*t1)));
		if dephasing_gamma < 0
			continue
		end
		thermal_gamma = 1 / t1;
		vx_maximum = minimum((1,0.5*sqrt(2 * thermal_gamma / (4*dephasing_gamma + thermal_gamma))-0.01));
		vx_minimum = maximum((0.1,vx_maximum-floor(10*vx_maximum)/10));
		for detune_index in param_space["detune_index"]
			detuning_freq = 10^detune_index;
			tend = 30*t2;
			# Ramsey interference setup
			simulation_params = (
				ramsey_detuned::SimulationType, # simulation type
				detuning_freq,
				dephasing_gamma,
				thermal_gamma,
				nothing, # ideal solution/trajectory data
			);
			problem = ODEProblem(lindblad, [1,0,0], (0.0, tend), simulation_params);
			# For some reason, problem with sampling rate
			ramsey_solution = solve(problem, alg_hints=[:stiff], saveat=sampling, abstol=abstol, reltol=reltol);
			ramsey_vy=[abs(ramsey_solution.u[i][2]) for i in range(1,length(ramsey_solution))];

			for buff_vx in range(vx_maximum,vx_maximum,step=-0.1)
				vx = buff_vx/2;
				# Ideal coherence magnitude preserving control
				v = [vx,0,sqrt(1-vx^2)]; # Checking vx/2
				simulation_params = (
					ideal::SimulationType, # simulation type
					detuning_freq,
					dephasing_gamma,
					thermal_gamma,
					nothing, # ideal solution/trajectory data
				);
				problem = ODEProblem(lindblad, v, (0.0, tend), simulation_params);
				ideal_solution = solve(problem, alg_hints=[:stiff], saveat=sampling, abstol=abstol, reltol=reltol);

				# Detuned coherence magnitude preserving control
				simulation_params = (
					detuned::SimulationType, # simulation type
					detuning_freq,
					dephasing_gamma,
					thermal_gamma,
					ideal_solution, # ideal solution/trajectory data
				);
				problem = ODEProblem(lindblad, v, (0.0, tend), simulation_params);
				detuned_solution = solve(problem, alg_hints=[:stiff], saveat=sampling, abstol=abstol, reltol=reltol);
				detuned_vy=[abs(detuned_solution.u[i][2]) for i in range(1,length(detuned_solution))];
				stable_h_ideal = get_hamiltonian(ideal_solution.u[end],(ideal::SimulationType,detuning_freq,dephasing_gamma,thermal_gamma,ideal_solution),ideal_solution.t[end]);
				stable_z_ideal = ideal_solution.u[end][3];
				lind_end = lindblad(ideal_solution.u[end],(ideal::SimulationType,detuning_freq,dephasing_gamma,thermal_gamma,ideal_solution),ideal_solution.t[end])[3];
				print("$(t1),$(t2),$(detuning_freq),$(vx),$(stable_z_ideal),$(stable_h_ideal[2]),$(lind_end),$(maximum(detuned_vy)),$(maximum(ramsey_vy)),$(detuned_solution.t[argmax(detuned_vy)]),$(ramsey_solution.t[argmax(ramsey_vy)]),$(detuned_vy[end]),$(ramsey_vy[end]),$(detuned_solution.t[end]/tend)\n");
			end
		end
	end
end














