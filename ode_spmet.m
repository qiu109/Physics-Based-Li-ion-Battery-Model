function [xdot,varargout]=ode_spmet(t,x,data,p)

U_n = x(1:(p.Nn-1));
U_p = x(p.Nn : 2*(p.Nn-1));
U_e = x(2*p.Nn-1 : 2*p.Nn-1+p.Nx-4);
T1 = x(end-2);
T2 = x(end-1);
T  = x(end);

cur=interp1(data.time,data.cur,t,[]);

%% Solid phase dynamics

% Molar flux for solid phase
J_p=-(cur./p.Area_p)./(p.Faraday*p.a_p*p.L_p);
J_n=(cur./p.Area_n)./(p.Faraday*p.a_n*p.L_n);

% Solid phase diffusivity temperature dependence
p.Ds_n = p.Ds_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/T1));
p.Ds_p = p.Ds_p0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/T1));

% Matrices for solid-phase Li concentration
[A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p]= matrixs(p);

% Calculation of the surface concentration
c_ss_p= C_p*U_p + D_p.*J_p;
c_ss_n= C_n*U_n + D_n.*J_n;

% Calculation of the solid particle concentration
c_p = A_p*U_p + B_p.*J_p;
c_n = A_n*U_n + B_n.*J_n;

%% Electrolyte phase dynamics

% Molar flux for solid phase
Jn=cur.*(1-p.t_plus )/(p.epsilon_e_n*p.Area_n*p.Faraday*p.L_n);
Jp=-cur.*(1-p.t_plus )/(p.epsilon_e_p*p.Area_p*p.Faraday*p.L_p);

% Matrices for electrolyte-phase Li concentration
[Amat_e]=matrixe(p);
B=ones(p.Nx-3,1);
B(1:p.Nxn-1,1)=Jn;
B(p.Nxn:p.Nxn-1+p.Nxs-1,1)=0;
B(p.Nxn-1+p.Nxs:p.Nx-3,1)=Jp;

% Calculation of the electrolyte phase concentration
c_e = Amat_e*U_e + B;

% Electrolyte concentrations
ce_n=U_e(1:p.Nxn,:);
ce_s=U_e((p.Nxn+1):(p.Nxn+p.Nxs+1),:);
ce_p=U_e((p.Nxn+p.Nxs+2):(p.Nxn+p.Nxs+p.Nxp-3),:);     

%Average electrolyte concentration
cen_bar = mean(ce_n);
ces_bar = mean(ce_s);
cep_bar = mean(ce_p);
c_e_bar = [cen_bar; ces_bar; cep_bar];
 
%% Calculation of potential of the cell

    % SOC of the electrodes  
     theta_p = c_ss_p/p.c_s_p_max;
     theta_n = c_ss_n/p.c_s_n_max;
     
    % OCP of the electrodes
     [Uref_p,dUpdT]=refpotantial_p (theta_p);
     [Uref_n,dUndT]=refpotantial_n (theta_n);
     
    % OCV of the cell
    
    V_ocv = Uref_p - Uref_n;
     
     % Kinetic reaction rate, adjusted for Arrhenius temperature dependence
     
     p.k_n = p.k_n0 * exp(p.E.kn/p.R*(1/p.T_ref - 1/T1)); 
     p.k_p = p.k_p0 * exp(p.E.kp/p.R*(1/p.T_ref - 1/T1)); 
     
    % Exchange current density
     i_0n = p.k_n .* ((p.c_s_n_max - c_ss_n) .* c_ss_n .* p.ce).^p.alph;
     i_0p = p.k_p .* ((p.c_s_p_max - c_ss_p) .* c_ss_p .* p.ce).^p.alph;
   
    % Overpotentials
     RTaF=(p.R*T1)./(p.alph*p.Faraday);
     eta_n  = RTaF .* asinh(cur ./ (2.*p.a_n.*p.Area_n.*p.L_n.*i_0n));
     eta_p  = RTaF .* asinh(-cur ./ (2.*p.a_p.*p.Area_p.*p.L_p.*i_0p));
     
    % SPM Voltage
     V_spm= eta_p - eta_n + Uref_p - Uref_n;
     
    % Overpotentials due to electrolyte subsystem
     kap_n0 = electrolyteCond(cen_bar);
     kap_s0 = electrolyteCond(ces_bar);
     kap_p0 = electrolyteCond(cep_bar);
     
    % Adjustment for Arrhenius temperature dependence
    kap_n = kap_n0 .* exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/T1));
    kap_s = kap_s0 .* exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/T1));
    kap_p = kap_p0 .* exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/T1));
    
    % Bruggeman relationships
     kap_n_eff = kap_n * p.epsilon_e_n.^(p.brug);
     kap_s_eff = kap_s * p.epsilon_e_s.^(p.brug);
     kap_p_eff = kap_p * p.epsilon_e_p.^(p.brug);
     
    % Activity coefficient
     dfca_n = electrolyteAct(cen_bar,T1,p);
     dfca_s = electrolyteAct(ces_bar,T1,p);
     dfca_p = electrolyteAct(cep_bar,T1,p);

    % Overpotential due to electrolyte conductivity
     V_electrolyteCond = (p.L_n./(2.*kap_n_eff) + 2*p.L_s./(2.*kap_s_eff)...
                          + p.L_p./(2.*kap_p_eff)).*cur;
    % Concentration at the boundaries
     ce0n=ce_n(1,:); cens=ce_n(end,:); cesp=ce_s(end,:); ce0p=ce_p(end,:);
     
    % Overpotential due to electrolyte polarization
     V_electrolytePolar = (2*p.R*T1)/(p.Faraday) .* (1-p.t_plus).* ...
        ( (1+dfca_n) .* (log(cens) - log(ce0n)) ...
         +(1+dfca_s) .* (log(cesp) - log(cens)) ...
         +(1+dfca_p) .* (log(ce0p) - log(cesp)));
     
    % Overall cell voltage
       V = V_spm + V_electrolyteCond + V_electrolytePolar;
       
    % Heat generation  
      Qohmic =  -cur.*(V - V_ocv);
      
      Qreaction= cur.*(eta_n - eta_p);
                      
      Qentropic= cur*T1*(dUpdT-dUndT)./1000;
      
      Qgen= Qohmic + Qreaction + Qentropic ;
      
     % Heat remove
     Qremv= p.h*p.A*(T-p.T_amb);
     
     % Temperature calculation
     T1_dot= Qgen./p.Cc + (T2-T1)./(p.Rc*p.Cc); %Tc core temp
     T2_dot= (p.T_amb-T2)./(p.Ru*p.Cs) - (T2-T1)./(p.Rc*p.Cs); %Ts surface tem.

     % Lumped Temperature calculation
     T_dot= (Qgen  - Qremv)./(p.M*p.Cp);
     
xdot = [c_n; c_p; c_e; T1_dot; T2_dot; T_dot]; 
varargout{1} = theta_p;
varargout{2} = theta_n;
varargout{3} = Uref_n;
varargout{4} = V;
varargout{5} = V_ocv;
varargout{6} = V_spm;
varargout{7} = Qgen;
varargout{8} = Qohmic;
varargout{9} = Qentropic;
varargout{10} =Qreaction;
varargout{11} =dUndT;
varargout{12} =dUpdT;



end
