%% Single Particle Model w Electrolyte and Thermal model

% Study based on publications 
% Modeling of a Commercial Graphite / LiFePO 4 Cell 158(5), 562–571. 
% Safari, M., & Delacourt, C. (2011). 
% https://doi.org/10.1149/1.3567007

% Battery State Estimation for a Single Particle Model with Electrolyte Dynamics. 1–16.
% Moura, S. J., Argomedo, F. B., Klein, R., Mirtabatabaei, A., & Krstic, M.
% DOI: 10.1109/TCST.2016.2571663


clear all;
close all;
clc;

run LFP_parameters
 load('dUdT.mat');
 load('SOC_ent.mat');

%% Finite difference for spherical particle and electrolyte
p.Np=50;
p.Nn=50;

p.delta_p =  p.R_p/(p.Np);
p.delta_n =  p.R_n/(p.Nn);  

p.Nxn=10; p.Nxp=10;  
p.Nxs=5; p.Nx=p.Nxn +p.Nxs +p.Nxp;
p.dx_n=1/p.Nxn;
p.dx_s=1/p.Nxs;
p.dx_p=1/p.Nxp;
p.Del_xn = p.L_n * p.dx_n;
p.Del_xs = p.L_s * p.dx_s;
p.Del_xp = p.L_p * p.dx_p;
%% Initial Conditions


% Fixed current controlled by C-rate

        p.C=1; 
        t_end=3600./p.C;
        t = 0:1:t_end;
        NT=length(t);
        I=p.C.*p.cell_Ah.*ones(size(t)); 

% Dynamic current input

%         load('UDDS_SCOTT.mat');
%         I = -current_exp';
%         C = I./2.3;
%         t = time_exp';
%         NT=length(t);
 
%         load('crnt.mat');
%         I = current'.*5;
%         C = I./2.3;
%         t =0:1:1369;
%         NT=length(t);


% Initial concentration  of solid particles and electrolyte 

Up0(1:p.Np-1,1)=p.c_s_p_max*p.theta_p_min;
Un0(1:p.Nn-1,1)=p.c_s_n_max*p.theta_n_max;
Ue0=1500.*ones(p.Nx-3,1);

% Temperature
T10 = p.T_amb;
T20 = p.T_amb;
T0  = p.T_amb;

%% Calculation of concentration of solid particles and temperature distribution

x0 = [Un0; Up0; Ue0;T10;T20;T0];
data.time = t;
data.cur = I;
Opt1= odeset('Events', @ode_cont1);
[t,x] = ode23s(@(t,x) ode_spmet(t,x,data,p),t,x0,Opt1);

U_n = x(:,1:(p.Nn-1));
U_p = x(:,p.Nn : 2*(p.Nn-1));
U_e = x(:,2*p.Nn-1 : 2*p.Nn-1+p.Nx-4);
T1 = x(:,end-2);
T2 = x(:,end-1);
T  = x(:,end);

%% Calculation of SOC of the cell based on concentrations of the electrodes
cn_avg= mean(U_n'); cp_avg= mean(U_p');
socn= ( (cn_avg/p.c_s_n_max)-p.theta_n_min ) / (p.theta_n_max-p.theta_n_min);
socp= ( (cp_avg/p.c_s_p_max)-p.theta_p_max ) / (p.theta_p_min-p.theta_p_max) ;

% Capacity 
capacity_n= (p.Area_n*p.L_n*p.Faraday*p.eps_s_n*p.c_s_n_max*p.theta_n_max)./3600; %Ah
capacity_p= (p.Area_p*p.L_p*p.Faraday*p.eps_s_p*p.c_s_p_max*p.theta_p_max)./3600; %Ah
capacity=((1:length(t)).*p.cell_Ah)./length(t);

%% Function outputs
NT = length(t);
theta_p = zeros(NT,1);
theta_n = zeros(NT,1);
Uref_n  = zeros(NT,1);
V       = zeros(NT,1);
V_ocv   = zeros(NT,1);
V_spm   = zeros(NT,1);

for k=1:NT
    
[~,theta_p(k),theta_n(k),Uref_n(k),V(k),V_ocv(k),V_spm(k),...
    Qgen(k),Qohmic(k),Qentropic(k),Qreaction(k),dUndT(k),dUpdT(k)]=ode_spmet(t(k),x(k,:)',data,p);
end

