clear all;
close all;
clc;

run LFP_parameters
load('dUdT.mat');
load('SOC_ent.mat');

%% Initial Conditions

    Upi=p.c_s_p_max*p.theta_p_min;
    Uni=p.c_s_n_max*p.theta_n_max;
    deltasei0=0;
    eps0=0.58;
    
% Fixed current controlled by C-rate
 
        p.C=1;   p.Nc=2;
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
Unn=zeros(p.Nc,NT);
tic;
for k=1:p.Nc
    if   1==mod(k,2)
        I=2.3*ones(size(t));
    else
        I=-2.3*ones(size(t));
    end
    [U_p,U_n,eps,delta_sei,V,R_tot_n,theta_n,c_ss_n]= spmet(p,t,I,Uni,Upi,deltasei0,eps0);
    
    Vspme(k,:)=V;
    deltasei(k,:)=delta_sei;
    Res_n(k,:)= R_tot_n;
    thetan(k,:)=theta_n;
    css_n(k,:)=c_ss_n;
    Unn(k,:)=U_n(:,end);
    
    Uni=U_n(end);
    Upi=U_p(end);
    deltasei0=delta_sei(end);
    eps0=eps(end);
    
  
end

toc;


%% Figures
% figure
% plot(t,V_ocv,t,V_spm,t,V);xlabel('time');ylabel('V');
% % legend('V_{ocv}','V_{spm}','V_{spme}');
