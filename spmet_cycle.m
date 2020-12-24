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
    Qi=(p.Area_n*p.L_n*p.Faraday*p.eps_s_n*p.c_s_n_max*p.theta_n_max);
% Fixed current controlled by C-rate
 
        p.C=1;   p.Nc=40;
       
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


%discharge charge
%     t=0:1:7200;
%     I=@(t) 2.3.*(t<=3600)-2.3.*(t>=3601);

%pulse
%     t=0:1:1890;
%     current=@(t) 2.3.*(t>=900)-2.3*(t>990);
     
% Discharge only
%     t=0:1:3600;
%     current=@(t) 2.3.*(t<=3600);
       
% plot(t,current(t));

tend=3600;
for i=1:p.Nc
    
    t=0:1:tend; NT=length(t);
    current=@(t) 2.3.*(t<=tend);
    
    if   1==mod(i,2)
        I=current(t);
    else 
        I=-current(t);
    end
    
    [c_p,c_n,eps,delta_sei,V,R_tot_n,theta_n,theta_p,T,Q_final,c_e,V_spm,V_ocv,c_ss_n1,c_ss_n]= spmet(p,t,I,Uni,Upi,deltasei0,eps0,Ti,Qi);

    
    for k=1:NT
        
        V_spme(i,k)=V(k);
        soc_n(i,k)=theta_n(k);
        core_T(i,k)=T(k);
        cn_s(i,k)=c_n(k,end)'; %last node.
        Qdecrease(i,k)=Q_final(k)/3600;
        
        
        %% State-of-Charge (Bulk), Li Concentration (Bulk)
        r_vec = (0:1/(p.Nn-2):1)';
        
        SOC_n(i,k) = 3/p.c_s_n_max * trapz(r_vec,r_vec.^2.*c_n(k,:)');
        SOC_p(i,k) = 3/p.c_s_p_max * trapz(r_vec,r_vec.^2.*c_p(k,:)');
        
        Cn_bulk(i,k)= 3*trapz(r_vec,r_vec.^2.*c_n(k,:)');
        
        Norm_Cn(i,k)= max(cn_s(i,k))/(p.c_s_n_max*p.theta_n_max) ;
        
        Qnorm(i,k)=Qdecrease(i,k)/((p.Area_n*p.L_n*p.Faraday*p.eps_s_n*p.c_s_n_max*p.theta_n_max)/3600);
  
        
         n_Li_s(k) = (3*p.eps_s_p*p.L_p*p.Area_p) * trapz(r_vec,r_vec.^2.*c_p(k,:)') ...
            + (3*p.eps_s_n*p.L_n*p.Area_n) * trapz(r_vec,r_vec.^2.*c_n(k,:)');
    end
    
    Decrease_Cn(i)=max(Norm_Cn(i,:)')*100;
    Qratio(i)=Qdecrease(i,k)/((p.Area_n*p.L_n*p.Faraday*p.eps_s_n*p.c_s_n_max*p.theta_n_max)/3600);
     
    if   1==mod(i,2)
        
        Uni=p.c_s_n_max*p.theta_n_min*Qratio(i);
        Upi=p.c_s_p_max*p.theta_p_max*Qratio(i);
    else
        Uni=p.c_s_n_max*p.theta_n_max*Qratio(i);
        Upi=p.c_s_p_max*p.theta_p_min*Qratio(i);
    end
    
    Qi=Q_final(end);
    deltasei0=delta_sei(end);
    eps0=eps(end);
    Ti=T(end);
    
%     dt(i)=Qdecrease(end)/1.0121/2.3*3600;
%     tend=dt;
end




toc;


