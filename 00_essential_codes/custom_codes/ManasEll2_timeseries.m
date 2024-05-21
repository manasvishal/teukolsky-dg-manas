function [y]=ManasEll2_timeseries(x,f0,c,x0,u) 
    X = u - x0;
    
    f = sin(f0*X).*exp(-c*X.^2);
    f1 = (f0*cos(f0*X) - 2*c*X.*sin(f0*X)).*exp(-c*X.^2);
    f2 = ( (4*(c^2).*X.^2 - 2*c - f0^2).*sin(f0*X) - 4*c*f0.*X.*cos(f0*X) ).*exp(-c*X.^2);
    f3 = ( 2*c.*X.*(-4*(c^2).*X.^2 + 6*c + 3*f0^2).*sin(f0*X) - f0*(-12*(c^2).*(X.^2)+6*c+f0^2).*cos(f0*X)).*exp(-c*X.^2);
    
    %%% Assemble outgoing initial data %%%
    % psi(t,r) = f'' + (3/r)f' + (3/r^2) f and psi_t and psi_r from chain rule
    
%     Psi   = f2 + (3./x).*f1 + (3./x.^2).*f;
%     Psi_t = f3 + (3./x).*f2 + (3./x.^2).*f1;
%     Psi_r = -f3 - (3./x.^2).*f1 - (3./x).*f2 - (3./x.^2).*f1 - (6./x.^3).*f;
% 
%     %%% Diagnostic: check violation of compact support %%%
%     XTemp = x(end-1):.01:2*x(end-1);
%     X = -XTemp - x0;
%     
%     f = sin(f0*X).*exp(-c*X.^2);
%     f1 = (f0*cos(f0*X) - 2*c*X.*sin(f0*X)).*exp(-c*X.^2);
%     f2 = ( (4*(c^2).*X.^2 - 2*c - f0^2).*sin(f0*X) - 4*c*f0.*X.*cos(f0*X) ).*exp(-c*X.^2);
%     f3 = ( 2*c.*X.*(-4*(c^2).*X.^2 + 6*c + 3*f0^2).*sin(f0*X) - f0*(-12*(c^2).*(X.^2)+6*c+f0^2).*cos(f0*X)).*exp(-c*X.^2);
%     
%     PsiTemp   = f2 + (3./XTemp).*f1 + (3./XTemp.^2).*f;
%     Psi_tTemp = f3 + (3./XTemp).*f2 + (3./XTemp.^2).*f1;
%     Psi_rTemp = -f3 - (3./XTemp.^2).*f1 - (3./XTemp).*f2 - (3./XTemp.^2).*f1 - (6./XTemp.^3).*f;
%     
%     temp = max(abs(PsiTemp(:))) + max(abs(Psi_tTemp(:))) + max(abs(Psi_rTemp(:)));
%     disp(['Maximum Violation of compact support = ',num2str(temp)])
%  
%     pause(2)
    

    y=f2;
end