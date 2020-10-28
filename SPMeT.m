function [varargout]= spmet(p,t,I,Uni,Upi,deltasei0,eps0)

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

% Initial concentration  of solid particles and electrolyte 

Up0(1:p.Np-1,1)=Upi;             %p.c_s_p_max*p.theta_p_min
Un0(1:p.Nn-1,1)=Uni;             %p.c_s_n_max*p.theta_n_max
Ue0=1500.*ones(p.Nx-3,1);

% Temperature
T10 = p.T_amb;
T20 = p.T_amb;
T0  = p.T_amb;

% SEI
delta_sei0=deltasei0;                    %0
eps=eps0;                                %0.58

%% Calculation of concentration of solid particles and temperature distribution

x0 = [Un0; Up0; Ue0;T10;T20;T0;delta_sei0;eps];
data.time = t;
data.cur = I;
Opt1= odeset('Events', @(t, U_n)ode_cont1(t, U_n,p));
[t,x] = ode23s(@(t,x) ode_spmet(t,x,data,p,dUdt,SOC_ent),t,x0,Opt1);

U_n = x(:,1:(p.Nn-1));
U_p = x(:,p.Nn : 2*(p.Nn-1));
U_e = x(:,2*p.Nn-1 : 2*p.Nn-1+p.Nx-4);
T1 = x(:,end-4);
T2 = x(:,end-3);
T  = x(:,end-2);
delta_sei= x(:,end-1);
eps=x(:,end);

%% Calculation of SOC of the cell based on concentrations of the electrodes
% cn_avg= mean(U_n'); cp_avg= mean(U_p');
% socn= ( (cn_avg/p.c_s_n_max)-p.theta_n_min ) / (p.theta_n_max-p.theta_n_min);
% socp= ( (cp_avg/p.c_s_p_max)-p.theta_p_max ) / (p.theta_p_min-p.theta_p_max) ;
% 
% % Capacity 
% capacity_n= (p.Area_n*p.L_n*p.Faraday*p.eps_s_n*p.c_s_n_max*p.theta_n_max)./3600; %Ah
% capacity_p= (p.Area_p*p.L_p*p.Faraday*p.eps_s_p*p.c_s_p_max*p.theta_p_max)./3600; %Ah
% capacity=((1:length(t)).*p.cell_Ah)./length(t);

%% Function outputs
% NT = length(t);
% theta_p = zeros(NT,1);
% theta_n = zeros(NT,1);
% Uref_n  = zeros(NT,1);
% V       = zeros(NT,1);
% V_ocv   = zeros(NT,1);
% V_spm   = zeros(NT,1);
NT=length(t);
for k=1:NT
    
[~,theta_p(k),theta_n(k),V(k),V_spm(k),V_ocv(k),...
    R_tot_n(k),eps(k),delta_sei(k)...
    ]...
    =ode_spmet(t(k),x(k,:)',data,p,dUdt,SOC_ent);

end

varargout{1} = U_p;
varargout{2} = U_n;
varargout{3} = eps;
varargout{4} = delta_sei;
varargout{5} = V;
varargout{6} = R_tot_n;
varargout{7} = theta_n;
end
