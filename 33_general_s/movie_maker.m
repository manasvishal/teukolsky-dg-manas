frm=1;
% rp=10;
close all;
set(gcf, 'Position', [360 360 720 720]);
set(gcf,'DefaultLineLineWidth',2)
%psiarr_unity=load("psiarr_unity.mat");
% v=VideoWriter("./Animation.mp4",'MPEG-4');
% v.FrameRate=24;
% open(v);
while frm<=length(output.time_arr)
        
        subplot(4,1,1)
        plot(x, real(flatten_array(output.psi_arr(frm,:,:),dg_globals)) ); 
        title("Psi")

        subplot(4,1,2)
        plot(x, real(flatten_array(output.pi_arr(frm,:,:),dg_globals)) ); 
        title("Pi")
        
        subplot(4,1,3)
        plot(x, real(flatten_array(output.phi_arr(frm,:,:),dg_globals)) ); 
        title("Phi")

        subplot(4,1,4)
        plot(x, dg_globals.rx.*(dg_globals.Dr*flatten_array(real(output.psi_arr(frm,:,:)),dg_globals)) ...
                -real(flatten_array(output.phi_arr(frm,:,:),dg_globals)) ); 
        title("Constraint")

        % plot(x, dg_globals.rx.*(dg_globals.Dr*flatten_array(real(sin(dg_globals.x+output.time_arr(1,frm))),dg_globals)) ...
        %         -real(flatten_array(cos(dg_globals.x+output.time_arr(1,frm)),dg_globals)) );
        
        
        % xlim([0,dg_globals.xR])
        % xlim([0,100])

        title1=sprintf("Tau = %1.2f",(output.time_arr(frm)));
        sgtitle(title1)
        hold on
        frame = getframe(gcf);
        size(frame.cdata);
%         writeVideo(v,frame);
        drawnow
        pause(1)
        hold off
        frm=frm+10;
end
% close(v)
close all