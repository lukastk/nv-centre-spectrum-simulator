function SZ = spin_matrix_z( s )
%SPIN_MATRIX_Z Constructs a spin matrix for the z-axis, for a specific
%total spin.

hbar = 1;%1.055e-34;

dim = 2*s + 1;
SZ = zeros(dim, dim);

for n = 1:dim
    SZ(n, n) = s - n + 1;
end

SZ = hbar*SZ;

end