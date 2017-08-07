function [H, V, E]  = n_spin_spectrum2(spins, zfs, kx, ky, kz,  zeeman)
    % N_SPIN_SPECTRUM Constructs the Hamiltonian of interacting spins.
    %   [H, V, D] = N_SPIN_SPECTRUM(spins, ks, H_mag) returns the Hamiltonian
    %   H, the eigenvectors V and the eigenvalues E. spins is a 1D array
    %   containing the spins of each interacting particle, ks is a 2D array
    %   containing the coupling between each spin. H_mag is the applied
    %   magnetic field on the particles.

    spin_multiplicities = 2*spins + 1;
    H_dim = prod(2*spins + 1);
    H = zeros( H_dim ); % Initial empty Hamiltonian

    % Create Hamiltonian
    % Adds each k_nm s_n*s_m term in the Hamiltonian

    for n = 1:length(spins)

        % Add S_X operators

        s1 = spin_matrix_x(spins(n));

        for m = (n+1):length(spins)
            s2 = spin_matrix_x(spins(m));

            if n ~= 1
                h = kron(kron_id_chain( spin_multiplicities(1:(n-1)) ), s1);
            else
                h = s1;
            end

            if m == n+1
                h = kron(h, s2);
            else
                h = kron(h, kron_id_chain( spin_multiplicities(n+1:m-1) ));
                h = kron(h, s2);
            end

            if m ~= length(spins)
                h = kron(h, kron_id_chain( spin_multiplicities(m+1:length(spins))));
            end

            H = H + kx(n, m)*h;
        end

        % Add S_Y operators

        s1 = spin_matrix_y(spins(n));

        for m = (n+1):length(spins)
            s2 = spin_matrix_y(spins(m));

            if n ~= 1
                h = kron(kron_id_chain( spin_multiplicities(1:(n-1)) ), s1);
            else
                h = s1;
            end

            if m == n+1
                h = kron(h, s2);
            else
                h = kron(h, kron_id_chain( spin_multiplicities(n+1:m-1) ));
                h = kron(h, s2);
            end

            if m ~= length(spins)
                h = kron(h, kron_id_chain( spin_multiplicities(m+1:length(spins))));
            end

            H = H + ky(n, m)*h;
        end

        % Add S_Z operators

        s1 = spin_matrix_z(spins(n));

        for m = (n+1):length(spins)
            s2 = spin_matrix_z(spins(m));

            if n ~= 1
                h = kron(kron_id_chain( spin_multiplicities(1:(n-1)) ), s1);
            else
                h = s1;
            end

            if m == n+1
                h = kron(h, s2);
            else
                h = kron(h, kron_id_chain( spin_multiplicities(n+1:m-1) ));
                h = kron(h, s2);
            end

            if m ~= length(spins)
                h = kron(h, kron_id_chain( spin_multiplicities(m+1:length(spins))));
            end

            H = H + kz(n, m)*h;
        end
        
        % Zero-field splitting
        
        s_zf = zfs(n)*s1*s1;
      
        if n == 1
            h = s_zf;
        else
            h = kron(kron_id_chain( spin_multiplicities(1:(n-1)) ), s_zf);
        end
        
        if n ~= length(spins)
            h = kron(h, kron_id_chain( spin_multiplicities(n+1:length(spins))));
        end
        
        H = H + h;

        % Add Zeeman Hamiltonian

        s = spin_matrix_z(spins(n))*zeeman(n);

        if n == 1
            h = s;
        else
            h = kron( kron_id_chain( spin_multiplicities(1:(n-1)) ), s);
        end

        if n ~= length(spins)
            h = kron(h, kron_id_chain( spin_multiplicities(n+1:length(spins)) ) );
        end

        H = H + h;

    end

    % Diagonalise the Hamiltonian

    [V, D] = eig(H);
    
    E = diag(D);

end
