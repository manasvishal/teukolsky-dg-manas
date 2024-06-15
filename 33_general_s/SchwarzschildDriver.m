
% NOTE: WaveRHS1D impliments two different versions of the wave equation 
% using hlayers.
%
% System 1: uses phi = \partial_{rho} psi
%           uses pi  = \partial_{tau} psi
%
% System 2: uses phi = \partial_{r} psi
%           uses pi  = \partial_{t} psi
%

Globals1D;

% disp(["The mode being solved is ",num2str(input.ell_mode_max),",",num2str(input.m_mode)])

N=dg_globals.N; K=dg_globals.K;
[Nv, VX, K, EToV] = MeshGen1D(dg_globals.xL,dg_globals.xR,dg_globals.K);

StartUp1D;

dg_globals.x=x;           dg_globals.Np=Np;        dg_globals.Nfp=Nfp; 
dg_globals.Nfaces=Nfaces; dg_globals.K=K;          dg_globals.vmapM=vmapM; 
dg_globals.vmapP=vmapP;   dg_globals.nx=nx;        dg_globals.rx=rx;
dg_globals.Dr=Dr;         dg_globals.LIFT=LIFT;    dg_globals.Fscale=Fscale;
dg_globals.rk4a=rk4a;     dg_globals.rk4b=rk4b;    dg_globals.rk4c=rk4c;
dg_globals.M=dg_globals.primary_mass;
dg_globals.a=dg_globals.primary_spin;



dg_globals.rho=x;
digits_VPA= 32;
symmetric_layers=0;
[idx1,idx2]=find_value_arr(x,dg_globals.locR);
dg_globals.locR=x(idx1,idx2);
dg_globals=HyperboloidalLayers(dg_globals,symmetric_layers,digits_VPA);


[dg_globals.r_sch,dg_globals.r_schm2m] = RstarToR(dg_globals.rstarCoord_double,dg_globals.M);



dg_globals=general_potential(dg_globals);




dg_globals.matIAB=cell(2,2);
dg_globals.matIAC=cell(2,2);
dg_globals.matIAB=cell(2,1);

dg_globals.matIAB{1,1}=dg_globals.capH_double./(1+dg_globals.capH_double);
dg_globals.matIAB{1,2}=1./(1+dg_globals.capH_double);
dg_globals.matIAB{2,1}=1./(1+dg_globals.capH_double);
dg_globals.matIAB{2,2}=dg_globals.capH_double./(1+dg_globals.capH_double);


dg_globals.prefactorA = -2.*dg_globals.spin_field.*((dg_globals.r_sch-3.*dg_globals.M)./dg_globals.r_sch.^2);
dg_globals.prefactorB = -2.*dg_globals.spin_field.*(dg_globals.r_sch-dg_globals.M)./dg_globals.r_sch.^2;

dg_globals.matIAC{1,1} = (1./(1-dg_globals.capH_double.^2)).*(-1*dg_globals.prefactorA);
dg_globals.matIAC{1,2}= (1./(1-dg_globals.capH_double.^2)).*(dg_globals.prefactorB);

dg_globals.matIAC{2,1}= dg_globals.capH_double.*dg_globals.matIAC{1,1};
dg_globals.matIAC{2,2}= dg_globals.capH_double.*dg_globals.matIAC{1,2};

dg_globals.matIAD{1,1}=(1./(1-dg_globals.capH_double.^2)).*dg_globals.potential_eff_general_s;
dg_globals.matIAD{2,1}=dg_globals.capH_double.*(1./(1-dg_globals.capH_double.^2)).*dg_globals.potential_eff_general_s;


if dg_globals.hyperboloidal_switch==1
    matIAC_11_inf=(2.*dg_globals.spin_field.*(dg_globals.omega_double-dg_globals.x.*dg_globals.omegaPrime_double)./(1+dg_globals.capH_double))...
            .*((dg_globals.rho-3.*dg_globals.M.*dg_globals.omega_double)./dg_globals.x.^2);
    dg_globals.matIAC{1,1}(end)=matIAC_11_inf(end);

    matIAC_12_inf=-(dg_globals.spin_field.*(dg_globals.omega_double-dg_globals.x.*dg_globals.omegaPrime_double)./(1+dg_globals.capH_double))...
        .*((2.*dg_globals.rho-2.*dg_globals.M)./dg_globals.x.^2);
    dg_globals.matIAC{1,2}(end)=matIAC_12_inf(end);
end


%%%%%%%%%%% Initial Conditions, set in the yml file
dg_globals=InitialConditions(dg_globals);

%%%%%%%%%%% particle at rp


[rp_row_idx,rp_col_idx]=find_value_arr(dg_globals.r_sch,dg_globals.rp);

dg_globals.particle.rp_row_idx=rp_row_idx;
dg_globals.particle.rp_col_idx=rp_col_idx;

%%% next two lines are needed to modify the fluxes near the particle
dg_globals.particle.rp_col_idx_left=rp_col_idx;
dg_globals.particle.rp_col_idx_right=rp_col_idx+1;

dg_globals.particle.rp_rsch=dg_globals.r_sch(rp_row_idx,rp_col_idx);

disp(["Particle is at r = ", num2str(dg_globals.particle.rp_rsch,16) , " and r* = ",num2str(x(rp_row_idx,rp_col_idx),16)])



%%%%%%%%%%%% setup the integrator
% filtering
dg_globals.matF=1;

dg_globals.xmin=min(min(diff(x(:,:))));

% dg_globals.dt=0.01046059193835873271649639093539;

dg_globals.dt = dg_globals.CFL*dg_globals.xmin;
% dg_globals.dt = 0.02;
% dg_globals.dt=round(dg_globals.dt*100)/100; %% round it off to two decimal places
if dg_globals.m_mode~=0
    dg_globals.dt=dg_globals.dt/(0.5*dg_globals.m_mode); %%% scaling for m mode
end
next_snapshot = dg_globals.DT;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SchwarzschildIntegrator;



