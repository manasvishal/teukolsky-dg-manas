% function [output] = SchwarzschildIntegrator(dg_globals) 

x=dg_globals.x;
Np=dg_globals.Np;
K=dg_globals.K;
rk4a=dg_globals.rk4a;
rk4b=dg_globals.rk4b;
rk4c=dg_globals.rk4c;
matF=dg_globals.matF;

time = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runge-Kutta residual storage  
resPsi = zeros(Np,K); resPi = zeros(Np,K); resPhi = zeros(Np,K); 

% compute time step size
%xmin = min(abs(x(1,:)-x(2,:)));
xmin=dg_globals.xmin;




disp(['xmin = ', num2str(xmin,8), ' dt = ', num2str(dg_globals.dt,8)])
% dt=1e-4
Nsteps = ceil(dg_globals.FinalTime/dg_globals.dt);

next_snapshot=dg_globals.DT;

Psi=dg_globals.IC_struct.psi_in;
Pi=dg_globals.IC_struct.pi_in;
Phi=dg_globals.IC_struct.phi_in;


%%% output structure with multi-dimensional arrays
%%%  first is t_step,second and third is the solution N,K
Nrecord=ceil(dg_globals.FinalTime/dg_globals.DT)+1;
output.psi_arr=zeros(Nrecord,Np,K);
output.pi_arr=zeros(Nrecord,Np,K);
output.phi_arr=zeros(Nrecord,Np,K);

output.time_arr=zeros(1,Nrecord);

output.psi_arr(1,:,:)=Psi;
output.pi_arr(1,:,:)=Pi;
output.phi_arr(1,:,:)=Phi;
output.time_arr(1,1)=time;
tic % used to estimte remaining simulation time
counter_snap=2;
for tstep=1:Nsteps
   for INTRK = 1:5
      timelocal = time + rk4c(INTRK)*dg_globals.dt;
      [rhsPsi,rhsPi, rhsPhi] = SchwarzschildRHS(Psi,Pi,Phi,dg_globals,timelocal);
       
      resPsi = rk4a(INTRK)*resPsi + dg_globals.dt*(rhsPsi);     
      resPi  = rk4a(INTRK)*resPi  + dg_globals.dt*(rhsPi);
      resPhi = rk4a(INTRK)*resPhi + dg_globals.dt*(rhsPhi);

      Psi = Psi + rk4b(INTRK)*resPsi;     
      Pi  = Pi  + rk4b(INTRK)*resPi;
      Phi = Phi + rk4b(INTRK)*resPhi;
      
      Psi=matF*Psi;
      Pi=matF*Pi;
      Phi=matF*Phi;
   
   end 
   % Increment time
   time = tstep*dg_globals.dt;
   
  
   % if mod(time,dg_globals.DT)<=dg_globals.dt/2 && tstep~=1
   if (mod(time,dg_globals.DT) <= dg_globals.dt && time ~= 0)
       output.psi_arr(counter_snap,:,:)=Psi(:,:);
       output.pi_arr(counter_snap,:,:)=Pi(:,:);
       output.phi_arr(counter_snap,:,:)=Phi(:,:);
  
       output.time_arr(1,counter_snap)=time;
       next_snapshot = next_snapshot + dg_globals.DT;
       counter_snap=counter_snap+1;
   end
    % estimate remaining simulation time
   percmsg=20;
   printTee=floor(dg_globals.FinalTime/percmsg);
   if mod(timelocal,printTee) <= dg_globals.dt && timelocal > 2*dg_globals.dt
        elapsed = toc/60;
        total_est = dg_globals.FinalTime*elapsed/timelocal;
        remain_est = total_est - elapsed;
        disp([num2str(timelocal*100/dg_globals.FinalTime),'% completed ',' | Program time: ',num2str(round(timelocal)) ,' | Elapsed clock(m) : ',num2str(elapsed),' | Estimated Time Remaining : ', num2str(remain_est) ])

   end
   

end
disp('Process completed succesfully')
% return
