function [ts, rhos] = solve_tdse(tspan, rho0, H)
    
    hbar = 1;%1.055e-34;

    rho_mat = zeros(length(H(0)));
    rho_len = length(rho_mat);

    function drho_dt = odefunc(t, rho)
        drho_dt = zeros(length(rho), 1);
        
        i = 1;
        for r = 1:rho_len
            for c = r:rho_len
                rho_mat(r, c) = rho(i);
                rho_mat(c, r) = rho(i)';
                i = i + 1;
            end
        end
        
        Ht = H(t);
        rho_mat = (- 1i / hbar) * (Ht*rho_mat - rho_mat*Ht);

        i = 1;
        for r = 1:rho_len
            for c = r:rho_len
                drho_dt(i) = rho_mat(r, c);
                i = i + 1;
            end
        end
    end

    %[ts, rhos] = ode45(@(t, rho) rho, tspan, rho0);
    
    [ts, rhos] = ode45(@odefunc, tspan, rho0);

end
