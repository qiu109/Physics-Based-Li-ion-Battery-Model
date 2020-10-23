function [xdot,varargout]=ode_spmet(t,x,data,p,dUdt,SOC_ent)

U_n = x(1:(p.Nn-1));
U_p = x(p.Nn : 2*(p.Nn-1));
U_e = x(2*p.Nn-1 : 2*p.Nn-1+p.Nx-4);
T1 = x(end-2);
T2 = x(end-1);
T  = x(end);

cur=interp1(data.time,data.cur,t,[]);
TEMP=T;

%% Solid phase dynamics

% Molar flux for solid phase
J_p=-(cur./p.Area_p)./(p.Faraday*p.a_p*p.L_p);
J_n=(cur./p.Area_n)./(p.Faraday*p.a_n*p.L_n);

% Solid phase diffusivity temperature dependence
p.Ds_n = p.Ds_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/TEMP));
p.Ds_p = p.Ds_p0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/TEMP));

% Matrices for solid-phase Li concentration
[A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p]= matrixs(p);

% Calculation of the surface concentration
c_ss_p= C_p*U_p + D_p.*J_p;
c_ss_n= C_n*U_n + D_n.*J_n;

%% Electrolyte phase dynamics

% Electrolyte phase diffusivity temperature dependence
p.De_n=p.De0_n*exp(p.E.De/p.R*(1/p.T_ref - 1/TEMP));
p.De_p=p.De0_p*exp(p.E.De/p.R*(1/p.T_ref - 1/TEMP));
p.De_s=p.De0_s*exp(p.E.De/p.R*(1/p.T_ref - 1/TEMP));

p.D_en_eff = p.De_n * p.epsilon_e_n^(p.brug-1);
p.D_es_eff = p.De_s * p.epsilon_e_s^(p.brug-1);
p.D_ep_eff = p.De_p * p.epsilon_e_p^(p.brug-1);

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
     [theta_p,theta_n]=getsoc(c_ss_p,c_ss_n,p);
     
    % OCP of the electrodes
     [Uref_p,dUpdT]=refpotantial_p (theta_p);
     [Uref_n,dUndT]=refpotantial_n (theta_n);
   
    % OCV of the cell
    
    V_ocv = Uref_p - Uref_n;
     
     % Kinetic reaction rate, adjusted for Arrhenius temperature dependence
     
     p.k_n = p.k_n0 * exp(p.E.kn/p.R*(1/p.T_ref - 1/TEMP)); 
     p.k_p = p.k_p0 * exp(p.E.kp/p.R*(1/p.T_ref - 1/TEMP)); 
     
    % Exchange current density
     i_0n = p.k_n .* ((p.c_s_n_max - c_ss_n) .* c_ss_n .* p.ce).^p.alph;
     i_0p = p.k_p .* ((p.c_s_p_max - c_ss_p) .* c_ss_p .* p.ce).^p.alph;
   
    % Overpotentials
     RTaF=(p.R*TEMP)./(p.alph*p.Faraday);
     eta_n  = RTaF .* asinh(cur ./ (2.*p.a_n.*p.Area_n.*p.L_n.*i_0n));
     eta_p  = RTaF .* asinh(-cur ./ (2.*p.a_p.*p.Area_p.*p.L_p.*i_0p));
     
    % SPM Voltage
     V_spm= eta_p - eta_n + Uref_p - Uref_n;
     
    % Overpotentials due to electrolyte subsystem
     kap_n0 = electrolyteCond(cen_bar);
     kap_s0 = electrolyteCond(ces_bar);
     kap_p0 = electrolyteCond(cep_bar);
     
    % Adjustment for Arrhenius temperature dependence
    kap_n = kap_n0 .* exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/TEMP));
    kap_s = kap_s0 .* exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/TEMP));
    kap_p = kap_p0 .* exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/TEMP));
    
    % Bruggeman relationships
     kap_n_eff = kap_n * p.epsilon_e_n.^(p.brug);
     kap_s_eff = kap_s * p.epsilon_e_s.^(p.brug);
     kap_p_eff = kap_p * p.epsilon_e_p.^(p.brug);
     
    % Activity coefficient
     dfca_n = electrolyteAct(cen_bar,TEMP,p);
     dfca_s = electrolyteAct(ces_bar,TEMP,p);
     dfca_p = electrolyteAct(cep_bar,TEMP,p);

    % Overpotential due to electrolyte conductivity
     V_electrolyteCond = (p.L_n./(2.*kap_n_eff) + 2*p.L_s./(2.*kap_s_eff)...
                          + p.L_p./(2.*kap_p_eff)).*cur;
    % Concentration at the boundaries
     ce0n=ce_n(1,:); cens=ce_n(end,:); cesp=ce_s(end,:); ce0p=ce_p(end,:);
     
    % Overpotential due to electrolyte polarization
     V_electrolytePolar = (2*p.R*TEMP)/(p.Faraday) .* (1-p.t_plus).* ...
        ( (1+dfca_n) .* (log(cens) - log(ce0n)) ...
         +(1+dfca_s) .* (log(cesp) - log(cens)) ...
         +(1+dfca_p) .* (log(ce0p) - log(cesp))); 
     
    % Overall cell voltage
     V = V_spm + V_electrolyteCond + V_electrolytePolar + ...
     (p.Rsei_p/(p.a_p*p.L_p) + p.Rsei_n/(p.a_n*p.L_n) )*cur;

    %% Degredation   
     
        % Exchange Current density of SEI 
        i0_sei_n= -p.Faraday*p.ksei*0.01;
        i0_sei_p= p.Faraday*p.ksei_p;
        
        % Overpotential of the SEI
        eta_sei_n= Uref_n + eta_n - p.Us + p.Rsei_n*p.L_sei*cur;
        eta_sei_p= Uref_p + eta_p - p.Us + p.Rsei_n*p.L_sei*cur;
        
        % Current density of the SEI
        Jsei_n= i0_sei_n.*exp(-p.alphasei_n*p.Faraday/(p.R*TEMP) * eta_sei_n);
        Jsei_p= i0_sei_p.*exp(p.alphasei_p*p.Faraday/(p.R*TEMP) *eta_sei_p);

        jn = (abs(J_n) - abs(Jsei_n)) * sign(J_n);
        jp = (abs(J_p)- abs(Jsei_p)) * sign(J_p);
        
        % Growth rate of SEI layer
        delta_sei_dot = -p.Msei/(p.Faraday*p.rhos) * Jsei_n;
        
        % Volume fraction
        eps_dot= 7.5e-7*Jsei_n;
        
        % Total resistance (film + growing SEI layer)
        R_tot_n = p.Rsei_n + delta_sei/p.kappa_s;
           
        %% Solid particle concentration after SEI layer formation
        c_p = A_p*U_p + B_p.*J_p;
        c_n = A_n*U_n + B_n.*jn;        
        
        %% State-of-Charge (Bulk)
        r_vec = (0:1/p.Nn:1)';
        c_s_n = [U_n(1); U_n; c_ss_n];
        c_s_p = [U_p(1); U_p; c_ss_p];
        SOC_n = 3/p.c_s_n_max * trapz(r_vec,r_vec.^2.*c_s_n);
        SOC_p = 3/p.c_s_p_max * trapz(r_vec,r_vec.^2.*c_s_p);
        
       %% Heat generation  
      Qohmic =  -cur.*(V - V_ocv);
      Qreaction=   cur.*(eta_n - eta_p) ;     
      dudT= interp1(SOC_ent,dUdt,theta_n,'linear','extrap');
      Qentropic=  cur*TEMP*(dudT)./1000;       %cur*TEMP*(dUpdT-dUndT)./1000  %cur*TEMP*(dudT)./1000
      Qgen= Qohmic + Qreaction + Qentropic;
      
     % Heat remove
     Qremv= p.h*p.A*(T-p.T_amb);
     
     % Temperature calculation
     T1_dot= Qgen./p.Cc + (T2-T1)./(p.Rc*p.Cc); %Tc core temp
     T2_dot= (p.T_amb-T2)./(p.Ru*p.Cs) - (T2-T1)./(p.Rc*p.Cs); %Ts surface tem. 
     
     % Lumped Temperature calculation
     T_dot= (Qgen  - Qremv)./(p.M*p.Cp);
     
 %% Capacity Calculations  
    
    Q_n= (p.Area_n*p.L_n*p.Faraday*eps*p.c_s_n_max*p.theta_n_max)./3600; %Ah
    Q_p= (p.Area_p*p.L_p*p.Faraday*p.eps_s_p*p.c_s_p_max*p.theta_p_max)./3600; %Ah
    
    % Capacity loss based on temp. and cycle numb. 
    Ah=p.Nc*(1-SOC_n)*Q_n;
    B=10000*(15/p.C)^(1/3);
    Qloss=B*exp( (-31700 + 370.3*p.C)/(p.R*TEMP) )*Ah^0.552;
   
    Qrem=(1-Qloss)*100; %???
    
     %% Outputs
    xdot = [c_n; c_p; c_e; T1_dot; T2_dot; T_dot; delta_sei_dot;eps_dot]; 

    varargout{1} = theta_p;
    varargout{2} = theta_n;
    varargout{3} = V;
    varargout{4} = V_spm;
    varargout{5} = V_ocv;
    varargout{6} = c_ss_p;
    varargout{7} = R_tot_n;
    varargout{8} = Qgen;
    varargout{9} = jn;
    varargout{10} =Q_n;
    varargout{11} =Qohmic;
    varargout{12} =Qreaction;
    varargout{13} =Qloss;
    varargout{14} =dudT;
    varargout{15} =c_ss_n;
    varargout{16} =p.Ds_n;
    varargout{17} =p.De_n;
    varargout{18} =p.k_n;
    varargout{19} =eta_n;
    varargout{20} =kap_n;
    varargout{21} =dfca_n;
    varargout{22} =V_electrolytePolar;
    varargout{23} =V_electrolyteCond;
end
