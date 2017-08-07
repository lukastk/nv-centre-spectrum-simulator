% Parameters for the NV-centre in contact with a spin-bath

spin_e_triplet = 1;
spin_N = 1;

NV_spins = [spin_e_triplet, spin_N]; % Spin of the electron triplet and the N-atom
NV_zfs = [1, 2]; % Zero-field splitting of triplet and N respectively.
NV_zeeman = [0, 0]; % Zeeman interactions
A_perp = 0.5; % hyperfine coupling constant in xy-plane
A_parallel = 0.7; % hyperfine coupling constant on the z-axis
NV_hf_coupling = [A_perp, A_perp, A_parallel];

% Parameters for the spin-bath

sb_spins = []; %ones(1, 3)*1/2;
sb_zfs = []; %zeros(1, 3);
sb_zeeman = []; %zeros(1,3);

% The coupling for the spin-spin interactions of the whole system

num_of_particles = 2 + length(sb_spins);

kx = ones(num_of_particles);
ky = ones(num_of_particles);
kz = ones(num_of_particles);

kx(1,2) = A_perp;
kx(2,1) = A_perp;
ky(1,2) = A_perp;
ky(2,1) = A_perp;
kz(1,2) = A_parallel;
kz(2,1) = A_parallel;

zfs = [NV_zfs sb_zfs]; % Vector containing the zero-field splittings of the whole system
spins = [NV_spins sb_spins]; % Vector containing all the spins of the whole system
zeeman = [NV_zeeman sb_zeeman]; % Vector containing all the zeeman constants for the whole system

H0 = n_spin_spectrum2(spins, zfs, kx, ky, kz, zeeman); % Define unperturbed Hamiltonian

% Simulation parameters

delta_H_amp = 1;
delta_H_op = spin_matrix_x(spin_e_triplet);
delta_H_op = kron(delta_H_op, kron_id_chain( 2*spins(2:length(spins)) + 1));
delta_H_op = delta_H_amp*delta_H_op;

start_pulse_freq = 0;
end_pulse_freq = 5e6;
delta_pulse_freq = 1e4;
pulse_phase = 0;

tspan = [0, 100];

num_indep_rho_components = (length(H0)^2 + length(H0))/2; % Number of linearly independent components in the density matrix.

initial_state = get_spin_state(spins, [0, 0]);
rho0 = initial_state * initial_state';

range = start_pulse_freq:delta_pulse_freq:end_pulse_freq;
probs = zeros(length(range), 1);
i = 1;

for pulse_freq = range
    
    delta_H = pulse_func(delta_H_op, pulse_freq, pulse_phase);
    [t, rho] = solve_tdse(tspan, rho0, @(t) H0 + delta_H(t));
    
    avg_prob = sum(abs(rho(:, 1)).^2) / (tspan(2) - tspan(1));
    probs(i) = avg_prob;
    
    i = i + 1;
end
