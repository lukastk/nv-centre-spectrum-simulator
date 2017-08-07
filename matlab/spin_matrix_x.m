function SX = spin_matrix_x( s )
%SPIN_MATRIX_X Constructs a spin matrix for the x-axis, for a specific
%total spin.

hbar = 1;%1.055e-34;

dim = 2*s + 1;
SX = zeros(dim, dim);

for n = 1:(dim-1)
    m = (-s + n - 1);
    
    SX(n, n+1) = (1/2)*sqrt(s*(s+1) - m*(m+1));
    SX(n+1, n) = (1/2)*sqrt(s*(s+1) - m*(m+1));
end

SX = hbar*SX;

end