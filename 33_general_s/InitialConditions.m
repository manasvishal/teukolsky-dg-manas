function [dg_globals]=InitialConditions(dg_globals)

if strcmpi(dg_globals.IC_type,'zero-momentum-gaussian')==1
    dg_globals.IC_struct.psi_in = (1/sqrt(2*pi*dg_globals.sigma^2))*exp(-(dg_globals.x-dg_globals.mu).^2/(2*dg_globals.sigma^2));
    dg_globals.IC_struct.pi_in = zeros(size(dg_globals.x));
    dg_globals.IC_struct.phi_in = (1/sqrt(2*pi*dg_globals.sigma^2))*exp(-(dg_globals.x-dg_globals.mu).^2/(2*dg_globals.sigma^2)).*(-2*(dg_globals.x-dg_globals.mu)./(2*dg_globals.sigma^2)); 

elseif strcmpi(dg_globals.IC_type,'non-zero-momentum-gaussian')==1
    %%%%%%%zero IC
    dg_globals.IC_struct.psi_in=zeros(size(dg_globals.x));
    dg_globals.IC_struct.pi_in=(1/sqrt(2*pi*dg_globals.sigma^2))*exp(-(dg_globals.x-dg_globals.mu).^2/(2*dg_globals.sigma^2));
    dg_globals.IC_struct.phi_in=zeros(size(dg_globals.x));

elseif strcmpi(dg_globals.IC_type,'gaussian')==1    

    dg_globals.IC_struct.psi_in=exp(-(dg_globals.x-dg_globals.mu).^2);
    dg_globals.IC_struct.pi_in=zeros(size(dg_globals.x));
    dg_globals.IC_struct.phi_in=dg_globals.rx.*(dg_globals.Dr*dg_globals.IC_struct.psi_in);


elseif strcmpi(dg_globals.IC_type,'zero')==1
    %%%%%%%zero IC
    dg_globals.IC_struct.psi_in=zeros(size(dg_globals.x));
    dg_globals.IC_struct.pi_in=zeros(size(dg_globals.x));
    dg_globals.IC_struct.phi_in=zeros(size(dg_globals.x));

elseif strcmpi(dg_globals.IC_type,'ell2')==1
    %%%%% SCOTT: exact outgoing initial data for ell=2
    dg_globals.IC.x0=-dg_globals.mu;
    dg_globals.IC.f0=2;
    dg_globals.IC.c=1;
    [Psi,Psi_r,Psi_t,Psi_tt] = OutgoingEll2_ID_SineGaussian(dg_globals.x,dg_globals.IC.f0,dg_globals.IC.c,dg_globals.IC.x0,-1*dg_globals.x);
    dg_globals.IC_struct.psi_in = Psi;
    dg_globals.IC_struct.pi_in = -Psi_t;
    dg_globals.IC_struct.phi_in = Psi_r;

end


return