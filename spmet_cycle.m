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
 
        p.C=1;   p.Nc=100;
        I=2.3;
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
Q_ni=2.3278; dt=3600;
        
for i=1:p.Nc
    
   
    t = 0:1:dt;
    
    [charge_current,discharge_current]=InputCurrent(dt,Q_ni); 
        
    
    if   1==mod(i,2)
        I=discharge_current.*ones(size(t));
    else
        I=charge_current.*ones(size(t));
    end
    
    [U_p,U_n,eps,delta_sei,V,R_tot_n,theta_n,c_ss_n,T]= spmet(p,t,I,Uni,Upi,deltasei0,eps0,Ti);
    
   
    
%     [SOCn(i,:),epss(i,:),delta_seii(i,:),R_neg(i,:),LumpT(i,:),V_spme(i,:)]=getvalues_cycle(i,p,t,eps,delta_sei,V,R_tot_n,theta_n,T);
    
    Q_n(i,:)= (p.Area_n*p.L_n*p.Faraday*max(eps)*p.c_s_n_max*max(theta_n))./3600;
    Q_ni=Q_n(end); 
    
    
    Uni=U_n(end);
    Upi=U_p(end);
    deltasei0=delta_sei(end);
    eps0=eps(end);
    Ti=T(end);
   dt=Q_ni/1.0121/2.3*3600;
   t=dt;
end

toc;

