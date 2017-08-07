function state = get_spin_state( spins, spin_projections)
%GET_SPIN_STATE Summary of this function goes here
%   Detailed explanation goes here

vals = (spin_projections + spins);
state = circshift(eye(spins(1)*2 + 1, 1), vals(1));

for i = 2:length(vals)
    state = kron(state, circshift(eye(spins(i)*2 + 1, 1), vals(i)));
end

end

