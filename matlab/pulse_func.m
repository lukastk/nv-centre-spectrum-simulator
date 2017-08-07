function [ p_f ] = pulse_func( A, freq, phase )
%PULSE_FUNC Returns a pulse function

    p_f = @(t) A * sin(freq*t + phase);

end

