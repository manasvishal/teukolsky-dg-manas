function [dg_globals]=general_potential(dg_globals)


dg_globals.pot_sch=(1-2*dg_globals.M./dg_globals.r_sch).*(-2*dg_globals.M./(dg_globals.r_sch.^3) - dg_globals.ell_mode_max*(dg_globals.ell_mode_max+1)./dg_globals.r_sch.^2);;


dg_globals.potential_eff_general_s = -(1-(2*dg_globals.M./dg_globals.r_sch)).*(2*dg_globals.M.*(1+dg_globals.spin_field)./(dg_globals.r_sch.^3) ...
                    + (dg_globals.ell_mode_max-dg_globals.spin_field)*(dg_globals.ell_mode_max+dg_globals.spin_field+1)./dg_globals.r_sch.^2);

dg_globals.pot_inf=((dg_globals.omega_double-dg_globals.x.*dg_globals.omegaPrime_double)./(1+dg_globals.capH_double))...
        .*((dg_globals.x-2*dg_globals.M*dg_globals.omega_double)./dg_globals.x)...
        .*(-(2*dg_globals.M.*(1+dg_globals.spin_field).*dg_globals.omega_double)./dg_globals.x.^3 ...
           - (dg_globals.ell_mode_max-dg_globals.spin_field).*(dg_globals.ell_mode_max+dg_globals.spin_field+1)./dg_globals.x.^2);

% dg_globals.potential_eff_general_s(Np,K)=double(dg_globals.pot_inf(Np,K));

dg_globals.potential=dg_globals.potential_eff_general_s;

end