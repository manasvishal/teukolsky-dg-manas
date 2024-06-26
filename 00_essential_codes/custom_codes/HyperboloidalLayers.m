function [struct_hyperboloidal] = HyperboloidalLayers(struct_hyperboloidal,symmetric_layers,num_digits)







digits(num_digits)
struct_hyperboloidal.locR_symbolic= vpa(struct_hyperboloidal.locR);
struct_hyperboloidal.s_symbolic = vpa(struct_hyperboloidal.xR); 
struct_hyperboloidal.P_symbolic = vpa(struct_hyperboloidal.P);
struct_hyperboloidal.x_symbolic = vpa(struct_hyperboloidal.rho);

C=struct_hyperboloidal.hyperboloidal_switch;

if C==0
 %%% no hyperboloidal layers
    struct_hyperboloidal.capH_double=0.*struct_hyperboloidal.rho;
    % struct_hyperboloidal.capHPrime_double=double(struct_hyperboloidal.capHPrime_symbolic);
    struct_hyperboloidal.omega_double=ones(size(struct_hyperboloidal.rho));
    struct_hyperboloidal.omegaPrime_double = 0.*struct_hyperboloidal.rho;
    % struct_hyperboloidal.omegaPrime_double=double(struct_hyperboloidal.omegaPrime_symbolic);
    % struct_hyperboloidal.ohm_double=double(struct_hyperboloidal.ohm_symbolic);
    % struct_hyperboloidal.ohmPrime_double=double(struct_hyperboloidal.ohmPrime_symbolic);
    struct_hyperboloidal.rstarCoord_double=struct_hyperboloidal.rho;

else 
    %%% hyperboloidal layers
    if symmetric_layers==1
    
    
        struct_hyperboloidal.omega_symbolic = 1+(-1).*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^( ...
                                      -1)).^struct_hyperboloidal.P_symbolic.*C.*abs((-1).*struct_hyperboloidal.locR_symbolic+ ...
                                      struct_hyperboloidal.rp_rstar+struct_hyperboloidal.x_symbolic).^struct_hyperboloidal.P_symbolic.* ...
                                      homeHVSD(abs((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.rp_rstar+ ...
                                      struct_hyperboloidal.x_symbolic));
        
        struct_hyperboloidal.omegaPrime_symbolic = (-1).*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^struct_hyperboloidal.P_symbolic.*0 ...
                        +(-1).*struct_hyperboloidal.P_symbolic.*((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic)^(-1).*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1) ...
                        *((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^((-1)+struct_hyperboloidal.P_symbolic).*homeHVSD((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic));
        
        struct_hyperboloidal.ohm_symbolic = struct_hyperboloidal.omega_symbolic.^2./(struct_hyperboloidal.omega_symbolic - abs(struct_hyperboloidal.x_symbolic).*struct_hyperboloidal.omegaPrime_symbolic);
        
        struct_hyperboloidal.ohmPrime_symbolic =abs(struct_hyperboloidal.x_symbolic).*((-2).*C.*struct_hyperboloidal.P_symbolic.*((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1) ...
          .*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^((-1)+struct_hyperboloidal.P_symbolic).*0+(-1).*C.*((-1)+struct_hyperboloidal.P_symbolic).* ...
          struct_hyperboloidal.P_symbolic.*((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-2).*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^( ...
          (-2)+struct_hyperboloidal.P_symbolic).*homeHVSD((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).*(1+(-1).*C.*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^( ...
          -1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^2.*(1+(-1).*C.*( ...
          ((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))+( ...
          -1).*abs(struct_hyperboloidal.x_symbolic).*((-1).*C.*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^struct_hyperboloidal.P_symbolic.* ...
          0+(-1).*C.*struct_hyperboloidal.P_symbolic.*((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*(((-1).* ...
          struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^((-1)+struct_hyperboloidal.P_symbolic).*homeHVSD((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic)))) ...
          .^(-2)+2.*((-1).*C.*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^struct_hyperboloidal.P_symbolic.* ...
          0+(-1).*C.*struct_hyperboloidal.P_symbolic.*((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*(((-1).* ...
          struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^((-1)+struct_hyperboloidal.P_symbolic).*homeHVSD((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))) ...
          .*(1+(-1).*C.*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^struct_hyperboloidal.P_symbolic.*homeHVSD( ...
          (-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).*(1+(-1).*C.*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))) ...
          .^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))+(-1).*abs(struct_hyperboloidal.x_symbolic).*((-1).*C.*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^( ...
          -1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^struct_hyperboloidal.P_symbolic.*0+(-1).*C.*struct_hyperboloidal.P_symbolic.*((-1) ...
          .*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*(((-1).*struct_hyperboloidal.locR_symbolic+struct_hyperboloidal.s_symbolic).^(-1).*((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic))).^((-1)+struct_hyperboloidal.P_symbolic) ...
          .*homeHVSD((-1).*struct_hyperboloidal.locR_symbolic+abs(struct_hyperboloidal.x_symbolic)))).^(-1);
        
        struct_hyperboloidal.capH_symbolic   = 1-struct_hyperboloidal.ohm_symbolic;
        
        struct_hyperboloidal.capHPrime_symbolic=-struct_hyperboloidal.ohmPrime_symbolic;
        
        struct_hyperboloidal.rstarCoord_symbolic=(1./struct_hyperboloidal.omega_symbolic).*(abs(struct_hyperboloidal.x_symbolic));
        
        struct_hyperboloidal.capH_double=double(struct_hyperboloidal.capH_symbolic);
        struct_hyperboloidal.capHPrime_double=double(struct_hyperboloidal.capHPrime_symbolic);
        struct_hyperboloidal.omega_double=double(struct_hyperboloidal.omega_symbolic);
        struct_hyperboloidal.omegaPrime_double=double(struct_hyperboloidal.omegaPrime_symbolic);
        struct_hyperboloidal.ohm_double=double(struct_hyperboloidal.ohm_symbolic);
        struct_hyperboloidal.ohmPrime_double=double(struct_hyperboloidal.ohmPrime_symbolic);
        struct_hyperboloidal.rstarCoord_double=double(struct_hyperboloidal.rstarCoord_symbolic);
    else
        
    
        struct_hyperboloidal.omega_symbolic = 1+(-1).*struct_hyperboloidal.hyperboloidal_switch.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)) ...
                                                  .^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic);
    
        struct_hyperboloidal.omegaPrime_symbolic = (-1).*struct_hyperboloidal.hyperboloidal_switch.* ...
                                                      struct_hyperboloidal.P_symbolic.*((-1).* ...
                                                      struct_hyperboloidal.locR_symbolic+ ...
                                                      struct_hyperboloidal.s_symbolic).^(-1).*(((-1).* ...
                                                      struct_hyperboloidal.locR_symbolic+ ...
                                                      struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                      struct_hyperboloidal.locR_symbolic+ ...
                                                      struct_hyperboloidal.x_symbolic)).^((-1)+ ...
                                                      struct_hyperboloidal.P_symbolic).*homeHVSD((-1).* ...
                                                      struct_hyperboloidal.locR_symbolic+ ...
                                                      struct_hyperboloidal.x_symbolic);
        struct_hyperboloidal.ohm_symbolic = (1+(-1).*struct_hyperboloidal.hyperboloidal_switch.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)) ...
                                                  .^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^2.*(1+ ...
                                                  struct_hyperboloidal.hyperboloidal_switch.* ...
                                                  struct_hyperboloidal.P_symbolic.*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).* ...
                                                  struct_hyperboloidal.x_symbolic.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^((-1)+ ...
                                                  struct_hyperboloidal.P_symbolic).*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)+(-1).* ...
                                                  struct_hyperboloidal.hyperboloidal_switch.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)) ...
                                                  .^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^(-1);
        struct_hyperboloidal.ohmPrime_symbolic =(-1).*struct_hyperboloidal.hyperboloidal_switch.*((-1)+ ...
                                                  struct_hyperboloidal.P_symbolic).* ...
                                                  struct_hyperboloidal.P_symbolic.*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-2).* ...
                                                  struct_hyperboloidal.x_symbolic.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^((-2)+ ...
                                                  struct_hyperboloidal.P_symbolic).*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic).*(1+(-1).* ...
                                                  struct_hyperboloidal.hyperboloidal_switch.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)) ...
                                                  .^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^2.*(1+ ...
                                                  struct_hyperboloidal.hyperboloidal_switch.* ...
                                                  struct_hyperboloidal.P_symbolic.*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).* ...
                                                  struct_hyperboloidal.x_symbolic.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^((-1)+ ...
                                                  struct_hyperboloidal.P_symbolic).*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)+(-1).* ...
                                                  struct_hyperboloidal.hyperboloidal_switch.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)) ...
                                                  .^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^(-2)+(-2).* ...
                                                  struct_hyperboloidal.hyperboloidal_switch.* ...
                                                  struct_hyperboloidal.P_symbolic.*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^((-1)+ ...
                                                  struct_hyperboloidal.P_symbolic).*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic).*(1+(-1).* ...
                                                  struct_hyperboloidal.hyperboloidal_switch.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)) ...
                                                  .^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).*(1+ ...
                                                  struct_hyperboloidal.hyperboloidal_switch.* ...
                                                  struct_hyperboloidal.P_symbolic.*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).* ...
                                                  struct_hyperboloidal.x_symbolic.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^((-1)+ ...
                                                  struct_hyperboloidal.P_symbolic).*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)+(-1).* ...
                                                  struct_hyperboloidal.hyperboloidal_switch.*(((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.s_symbolic).^(-1).*((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)) ...
                                                  .^struct_hyperboloidal.P_symbolic.*homeHVSD((-1).* ...
                                                  struct_hyperboloidal.locR_symbolic+ ...
                                                  struct_hyperboloidal.x_symbolic)).^(-1);
        struct_hyperboloidal.capH_symbolic   = 1-struct_hyperboloidal.ohm_symbolic;
        
        struct_hyperboloidal.capHPrime_symbolic=-struct_hyperboloidal.ohmPrime_symbolic;
        
        struct_hyperboloidal.rstarCoord_symbolic=(1./struct_hyperboloidal.omega_symbolic).*(struct_hyperboloidal.x_symbolic);
        
        struct_hyperboloidal.capH_double=double(struct_hyperboloidal.capH_symbolic);
        struct_hyperboloidal.capHPrime_double=double(struct_hyperboloidal.capHPrime_symbolic);
        struct_hyperboloidal.omega_double=double(struct_hyperboloidal.omega_symbolic);
        struct_hyperboloidal.omegaPrime_double=double(struct_hyperboloidal.omegaPrime_symbolic);
        struct_hyperboloidal.ohm_double=double(struct_hyperboloidal.ohm_symbolic);
        struct_hyperboloidal.ohmPrime_double=double(struct_hyperboloidal.ohmPrime_symbolic);
        struct_hyperboloidal.rstarCoord_double=double(struct_hyperboloidal.rstarCoord_symbolic);
        
    end

%%% end for switch C
end

end