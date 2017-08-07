function [ Us ] = pert_teo_defunct( H0, delta_H, t, delta_t )
% PERT_TEO Returns a set of the time-evolution operators at discrete
% time-steps, with a Hamiltonian under a time-dependent perturbation.

U = eye(size(H0));
n = round(t / delta_t);
delta_t = t / n;
hbar = 1;

Us = cell(n, 1);

for k = 1:n
    
    H = H0 + delta_H(k*delta_t*k);
    dU = expm( (- 1i / hbar ) * H * delta_t );
    U = U * dU;
    
    Us{k} = U;
    
end

end