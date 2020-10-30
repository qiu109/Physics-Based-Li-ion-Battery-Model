clear all;
close all;
clc;

run LFP_parameters
load('dUdT.mat');
load('SOC_ent.mat');

%% Initial Conditions
p.Np=50;
p.Nn=50;
    Upi=p.c_s_p_max*p.theta_p_min*ones(p.Np-1,1);
    Uni=p.c_s_n_max*p.theta_n_max*ones(p.Nn-1,1);
    deltasei0=0;
    eps0=0.58;
    Ti=298.15;
% Fixed current controlled by C-rate
 
        p.C=1;   p.Nc=4;
        t_end=3599;
        t = 0:1:t_end;
        NT=length(t);
        
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


Vspme=zeros(p.Nc,NT);
deltasei=zeros(p.Nc,NT);
Res_n=zeros(p.Nc,NT);
thetan=zeros(p.Nc,NT);
css_n=zeros(p.Nc,NT);
Un_s=zeros(p.Nc,NT);
Dsn=zeros(p.Nc,NT);
l_T=zeros(p.Nc,NT);
Ds_n=zeros(p.Nc,NT);
tic;
for k=1:p.Nc
    if   1==mod(k,2)
        I=2.3*ones(size(t));
    else
        I=-2.3*ones(size(t));
    end
    [U_p,U_n,eps,delta_sei,V,R_tot_n,theta_n,c_ss_n,p.Ds_n,D_n,T]= spmet(p,t,I,Uni,Upi,deltasei0,eps0,Ti);
    
    Vspme(k,:)=V;
    deltasei(k,:)=delta_sei;
    Res_n(k,:)= R_tot_n;
    thetan(k,:)=theta_n;
    css_n(k,:)=c_ss_n;
    Un_s(k,:)=U_n(:,end);
    Dsn(k,:)=p.Ds_n;
    Ds_n(k,:)=D_n;
    l_T(k,:)=T;
    
    Q_n(k,:)= (p.Area_n*p.L_n*p.Faraday*eps*p.c_s_n_max*max(theta_n))./3600;
     
    Uni=U_n(end,:)';
    Upi=U_p(end,:)';
    deltasei0=delta_sei(end);
    eps0=eps(end);
    Ti=T(end);
  
end

toc;


%% Figures
% figure
% plot(t,V_ocv,t,V_spm,t,V);xlabel('time');ylabel('V');
% % legend('V_{ocv}','V_{spm}','V_{spme}');
