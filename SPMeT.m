function [varargout]= spmet(p,t,I,Uni,Upi,deltasei0,eps0,Ti,Qi)

load('dUdT.mat');
load('SOC_ent.mat');

%% Finite difference for spherical particle and electrolyte

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
Ue0=1200.*ones(p.Nx-3,1);

% Temperature
T10 = p.T_amb;
T20 = p.T_amb;
T0  = Ti;

% SEI
delta_sei0=deltasei0;                    %0
eps=eps0;                                %0.58


Qdot0=Qi; 
%% Calculation of concentration of solid particles and temperature distribution

x0 = [Un0; Up0; Ue0;T10;T20;T0;delta_sei0;eps;Qdot0];
data.time = t;
data.cur = I;
Opt1= odeset('Events', @(t,socdot)ode_cont1(t,socdot));
[ts,x] = ode23s(@(t,x) ode_spmet_degr_cycle(t,x,data,p,dUdt,SOC_ent),t,x0,Opt1);

c_n = x(:,1:(p.Nn-1));
c_p = x(:,p.Nn : 2*(p.Nn-1));
c_e = x(:,2*p.Nn-1 : 2*p.Nn-1+p.Nx-4);
T1 = x(:,end-5);
T2 = x(:,end-4);
T  = x(:,end-3);
delta_sei= x(:,end-2);
eps=x(:,end-1);
Q_final=x(:,end);


NT=length(t);



for k=1:NT
    
[~,theta_p(k),theta_n(k),V(k),V_spm(k),V_ocv(k),...
    R_tot_n(k),eps(k),delta_sei(k),c_ss_n1(k),c_ss_n(k)...
    ]...
    =ode_spmet_degr_cycle(t(k),x(k,:)',data,p,dUdt,SOC_ent);


end

varargout{1} = c_p;
varargout{2} = c_n;
varargout{3} = eps;
varargout{4} = delta_sei;
varargout{5} = V;
varargout{6} = R_tot_n;
varargout{7} = theta_n;
varargout{8} = theta_p;
varargout{9} = T;
varargout{10}=Q_final;
varargout{11}=c_e;
varargout{12}=V_spm;
varargout{13}=V_ocv;
varargout{14}=c_ss_n1;
varargout{15}=c_ss_n;
end
