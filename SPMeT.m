function [varargout]= spmet(p,t,I,Uni,Upi,deltasei0,eps0,Ti)

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

Up0=Upi*ones(p.Np-1,1);             %p.c_s_p_max*p.theta_p_min
Un0=Uni*ones(p.Nn-1,1);             %p.c_s_n_max*p.theta_n_max
Ue0=1500.*ones(p.Nx-3,1);

% Temperature
T10 = p.T_amb;
T20 = p.T_amb;
T0  = Ti;

% SEI
delta_sei0=deltasei0;                    %0
eps=eps0;                                %0.58

%% Calculation of concentration of solid particles and temperature distribution

x0 = [Un0; Up0; Ue0;T10;T20;T0;delta_sei0;eps];
data.time = t;
data.cur = I;
Opt1= odeset('Events', @(t, c_n)ode_cont1(t, c_n,p));
[ts,x] = ode23s(@(t,x) ode_spmet_degr_cycle(t,x,data,p,dUdt,SOC_ent),t,x0,Opt1);

c_n = x(:,1:(p.Nn-1));
c_p = x(:,p.Nn : 2*(p.Nn-1));
c_e = x(:,2*p.Nn-1 : 2*p.Nn-1+p.Nx-4);
T1 = x(:,end-4);
T2 = x(:,end-3);
T  = x(:,end-2);
delta_sei= x(:,end-1);
eps=x(:,end);



[theta_p,theta_n,V,V_spm,V_ocv,...
    R_tot_n,eps,delta_sei,c_ss_n]=getvalues_spmet(t,x,data,p,dUdt,SOC_ent);

varargout{1} = c_p;
varargout{2} = c_n;
varargout{3} = eps;
varargout{4} = delta_sei;
varargout{5} = V;
varargout{6} = R_tot_n;
varargout{7} = theta_n;
varargout{8} = c_ss_n;
varargout{9} = T;
end
