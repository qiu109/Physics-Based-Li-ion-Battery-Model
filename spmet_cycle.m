clear all;
close all;
clc;

run LFP_parameters
load('dUdT.mat');
load('SOC_ent.mat');

%% Initial Conditions
p.Np=50;
p.Nn=50;
    Upi=p.c_s_p_max*p.theta_p_min;
    Uni=p.c_s_n_max*p.theta_n_max;
    deltasei0=0;
    eps0=0.58;
    Ti=298.15;
% Fixed current controlled by C-rate
 
        p.C=1;   p.Nc=1;
       
% Dynamic current input

%         load('UDDS_SCOTT.mat');
%         I = -current_exp';
%         C = I./2.3;
%         t = time_exp';
%         NT=length(t);
 
%         load('crnt.mat');
%         I = current';
%         p.C = I./2.3;
%         t =0:1:1369;
%         NT=length(t);

%% Cycle
tic;
Q_ni=2.3278; 

%discharge charge
    t=0:1:7800;
    I=@(t) 2.3.*(t<=3600)-2.3.*(t>=4200);

%pulse
%     t=0:1:1890;
%     I=@(t) 2.3.*(t>=900)-2.3*(t>990);
    
% Discharge only
%     t=0:1:3600;
%     I=@(t) 2.3.*(t<=3600);
       
plot(t,I(t));
NT=length(t);

for i=1:p.Nc
    
      
    [c_p,c_n,eps,delta_sei,V,R_tot_n,theta_n,theta_p,T]= spmet(p,t,I,Uni,Upi,deltasei0,eps0,Ti);
    
    
    for k=1:NT
        
        V_spme(i,k)=V(k);
        soc_n(i,k)=theta_n(k);
        core_T(i,k)=T(k);
        cn_s(i,k)=c_n(k,end)'; %last node.
        
        %% State-of-Charge (Bulk), Li Concentration (Bulk)
        r_vec = (0:1/(p.Nn-2):1)';
        
        SOC_n(i,k) = 3/p.c_s_n_max * trapz(r_vec,r_vec.^2.*c_n(k,:)');
        SOC_p(i,k) = 3/p.c_s_p_max * trapz(r_vec,r_vec.^2.*c_p(k,:)');
        
        Cn_bulk(i,k)= 3*trapz(r_vec,r_vec.^2.*c_n(k,:)');
        
        Norm_Cn(i,k)= max(cn_s(i,k))/(p.c_s_n_max*p.theta_n_max) ;
        
        
    end
    
    Decrease_Cn(i)=max(Norm_Cn(i,:)')*100;
    Q_n(i)= (p.Area_n*p.L_n*p.Faraday*max(eps)*p.c_s_n_max*max(theta_n))./3600;
    
    if   1==mod(i,2)
        
        Uni=p.c_s_n_max*min(theta_n);
        Upi=p.c_s_p_max*max(theta_p);
    else
        Uni=p.c_s_n_max*max(theta_n);
        Upi=p.c_s_p_max*min(theta_p);
    end
    
    Q_ni=Q_n(end);
    %     Uni=Un(end);
    %     Upi=c_p(end);
    deltasei0=delta_sei(end);
    eps0=eps(end);
    Ti=T(end);
    dt=Q_ni/1.0121/2.3*3600;
    %    t=dt;
end




toc;


%% Figures
% figure
% plot(t,V_ocv,t,V_spm,t,V);xlabel('time');ylabel('V');
% % legend('V_{ocv}','V_{spm}','V_{spme}');


%% Figures
% figure
% plot(t,V_ocv,t,V_spm,t,V);xlabel('time');ylabel('V');
% % legend('V_{ocv}','V_{spm}','V_{spme}');
