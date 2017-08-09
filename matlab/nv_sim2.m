digits(32)

% Physical constants

hbar = 1;
B0 = 0*513e-4; % [T] Applied DC field
gamma_e = 1.7609e11; % [rad s^-1 T^-1] Electron gyromagnetic ratio.
gamma_N14 = 19.331e16; % [rad s^-1 T^-1 N-14 gyromagnetic ratio.
omega_e = - gamma_e * B0 / (2 * pi); % [Hz] Zeeman splitting of electron triplet
omega_N14 = - gamma_N14 * B0 / (2 * pi); % [Hz] Zeeman splitting of N14
Dzfs = 2870e6; % [Hz] Zero field splitting of electron triplet.
P = - 4.95e6; % [Hz] Zero field splitting of N14.

A_perp = 0; % [Hz] Hyperfine coupling constant, perpendicular to [1, 1, 1] axis. Can be neglected under secular approximation.
A_parallel = -2.16e6; % [Hz] hyperfine coupling constant on the z-axis

% Calculation parameters

phi = pi/4; % Phase of the perturbation.
%f_MW_start = 500e6; % [Hz]
%f_MW_end = 530e7; % [Hz]
f_MW_start = 2.8e9; % [Hz]
f_MW_end = 2.955e9; % [Hz]
calculations = 20000;
f_MW_step = (f_MW_end - f_MW_start) / calculations;

tstart = 0;
time_steps = 1000;
%tend = 1/(f_MW_end) * 10;
%delta_t = (tend - tstart) / time_steps;

% Parameters fodr the NV-centre in contact with a spin-bath

spin_e_triplet = 1;
spin_N = 1;

NV_spins = [spin_e_triplet, spin_N]; % Spin of the electron triplet and the N-atom
NV_zfs = [Dzfs, -omega_N14]; % Zero-field splitting of triplet and N respectively.
NV_zeeman = [omega_e, -omega_N14]; % Zeeman interactions

% Parameters for the spin-bath

sb_spins = []; %ones(1, 3)*1/2;
sb_zfs = []; %zeros(1, 3);
sb_zeeman = []; %zeros(1,3);

% The coupling for the spin-spin interactions of the whole system

num_of_particles = 2 + length(sb_spins);

kx = zeros(num_of_particles);
ky = zeros(num_of_particles);
kz = zeros(num_of_particles);

kx(1,2) = A_perp;
kx(2,1) = A_perp;
ky(1,2) = A_perp;
ky(2,1) = A_perp;
kz(1,2) = A_parallel;
kz(2,1) = A_parallel;

zfs = [NV_zfs sb_zfs]; % Vector containing the zero-field splittings of the whole system
spins = [NV_spins sb_spins]; % Vector containing all the spins of the whole system
zeeman = [NV_zeeman sb_zeeman]; % Vector containing all the zeeman constants for the whole system

spin_multiplicities = 2*spins+1;

H0 = (2*pi)*n_spin_spectrum2(spins, zfs, kx, ky, kz, zeeman); % Define unperturbed Hamiltonian

% Create density matrix

[V, D] = eig(H0);
E = diag(D);

[r, c] = find(E==min(E), 1);
vec = V(:,r);

rho0 = vec * vec';

% Electron triplet basis states

e_basis = { [1;0;0], [0;1;0], [0;0;1] };

% Initialize calculation

omega_1 = 7e6;

zfs_e_op = kron(spin_matrix_z(1) * spin_matrix_z(1), kron_id_chain(spin_multiplicities(2:length(spin_multiplicities))));
sz_op = kron(spin_matrix_z(1), kron_id_chain(spin_multiplicities(2:length(spin_multiplicities))));

Hx = (1 / sqrt(2))* [ 0 1 0 ; 1 0 1 ; 0 1 0];
Hx = kron(Hx, kron_id_chain(spin_multiplicities(2:length(spin_multiplicities))));

Hy = (1 / sqrt(2))* [ 0 -1i 0 ; 1i 0 1i ; 0 -1i 0];
Hy = kron(Hy, kron_id_chain(spin_multiplicities(2:length(spin_multiplicities))));

avg_probs = zeros(calculations, 1);

for i = 1:calculations
    
    f_MW = f_MW_start + f_MW_step * i;
    %f_MW = 1.431e9;
    % Add perturbation

    H_e = - f_MW * zfs_e_op;

    H = H0 + (2*pi)*( H_e + (omega_1 / sqrt(2))*(cos(phi) * Hx + sin(phi) * Hy) );
    
    rho = rho0;
    probs = zeros(time_steps, 1);
    
    tend = 2e3/f_MW;
    delta_t = (tend - tstart) / time_steps;
    
    for t = 0:(time_steps-1)
        U = expm( - 1i * H * delta_t * t);
        rho = U*rho0*U';
        
        % Partial trace
        
        rho_A = zeros(3);
        
        for j = 1:3
            for k = 1:3
                for l=1:3
                    rho_A = rho_A + rho( 3*(j-1) + l , 3*(k - 1) + l ) * ( e_basis{j} * e_basis{k}' );
                end
            end
        end
        
        p = e_basis{2}' * rho_A * e_basis{2};
        
        probs(t+1) = p;
    end
    
    %for t = 2:time_steps
    %    rho = U*rho*U';
    %    
    %    p = trace( sz_op * rho );
    %    probs(t, 1) = p;
    %end

    avg_probs(i) = real(mean(probs));
    
    fprintf("Calculation " + i + " of " + calculations + "\n");
end

%ts = 0:delta_t:(tend-delta_t);
%plot(ts*1e6, real(probs))

fs = (f_MW_start+f_MW_step):f_MW_step:f_MW_end;
plot(fs, avg_probs)
save_dat = [fs' avg_probs];
filename = "./results/matlab_nv_sim" + datestr(now,'mm.dd.yyyy-HH.MM.SS');
save(filename + ".dat", "save_dat");

title('Graph of NV-center ESR. Applied DC field: ' + B0 + "T");
xlabel('Pulse frequency [Hz]');
ylabel('Prob(T_0)');
pl = plot(fs, avg_probs);
saveas(pl,char(filename + ".pdf"));
saveas(pl,char(filename + ".fig"));
