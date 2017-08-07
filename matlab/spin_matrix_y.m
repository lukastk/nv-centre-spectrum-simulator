function SY = spin_matrix_y( s )
%SPIN_MATRIX_Y Constructs a spin matrix for the y-axis, for a specific
%total spin.

hbar = 1;%1.055e-34;

dim = 2*s + 1;
SY = zeros(dim, dim);

for n = 1:(dim-1)
    m = (-s + n - 1);
    
    SY(n, n+1) = (1 / ( 2 * 1i ) )*sqrt(s*(s+1) - m*(m+1));
    SY(n+1, n) = - (1 / ( 2 * 1i ) )*sqrt(s*(s+1) - m*(m+1));
end

SY = hbar*SY;

end