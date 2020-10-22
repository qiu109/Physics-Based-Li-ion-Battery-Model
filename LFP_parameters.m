%% Capacity

p.cell_Ah=2.3;          %Capacity of the battery

%% Geometric Params
% Radius of particles

p.R_p = 3.65e-8;   % Radius of solid particles in positive electrode [m]
p.R_n = 3.5e-6;    % Radius of solid particles in negative electrode [m]

% Thiknesses of layers

p.L_p = 7e-5;      % Thickness of positive electrode [m]
p.L_n = 3.4e-5;    % Thickness of negative electrode [m]
p.L_s = 2.5e-5;
% Volume fractions of electrodes

p.eps_s_p = 0.428; % Volume fraction in solid for pos. electrode
p.eps_s_n = 0.58;  % Volume fraction in solid for neg. electrode 

p.epsilon_e_n = 0.503;      % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 0.55;      % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.633;      % Volume fraction in electrolyte for pos. electrode

p.epsilon_e_filler_p = 0.0535;      % Filler volume fraction in electrolyte for pos. electrode
p.epsilon_e_filler_n = 0.0326;      % Filler Volume fraction in electrolyte for neg. electrode
% Specific interfacial surface area

p.a_p = (3*p.eps_s_p)/p.R_p;  % Positive electrode [m^2/m^3]
p.a_n = (3*p.eps_s_n)/p.R_n;  % Negative electrode [m^2/m^3]

% Electrode area in C26650

p.Area_p = 0.1694;    % Positive electrode area [cm^2]
p.Area_n = 0.1755;    % Negative electrode area [cm^2]

%% Transport Params

% Diffusion coefficient in particles
p.Ds_p0 = 1.18e-18;  % Diffusion coeff for solid in pos. electrode, [m^2/s]
p.Ds_n0 = 2e-14;    % Diffusion coeff for solid in neg. electrode, [m^2/s]

p.Faraday = 96485.33289;      % Faraday's constant, [Coulumbs/mol]
p.t_plus = 0.5;             % Transference number
p.brug = 1.5;                 % Bruggeman porosity

p.De_s=(9e-15);             % Diffusion coeff for solid in sep, [m^2/s]
p.De_p=(10e-12);             % Diffusion coeff for solid in pos, [m^2/s]
p.De_n=(10e-12);             % Diffusion coeff for solid in neg, [m^2/s]

p.D_en_eff = p.De_n * p.epsilon_e_n^(p.brug-1);
p.D_es_eff = p.De_s * p.epsilon_e_s^(p.brug-1);
p.D_ep_eff = p.De_p * p.epsilon_e_p^(p.brug-1);

%% Kinetic Params

p.R = 8.314472;       % Gas constant, [J/mol-K]
p.alph = 0.5;         % Charge transfer coefficients

% Reaction rates
p.k_n0 = 8.19e-7;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
p.k_p0 = 3e-7; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]

%% Concentrations

p.c_s_p_max = 22806;    % Max concentration in cathode, [mol/m^3]  
p.c_s_n_max = 31370;    % Max concentration in anode, [mol/m^3]
p.ce = 1000;            % Fixed electrolyte concentration for SPM, [mol/m^3]

% Stochiometry
p.theta_n_min=0;              % Min soc of the negative electrode
p.theta_n_max=0.8;            % Max soc of the negative electrode

p.theta_p_min=0.03;           % Min soc of the positive electrode
p.theta_p_max=0.76;           % Max soc of the positive electrode

%% Thermodynamic Properties


p.T_ref = 298.15;       % [K]
p.T_amb= 263.15;

p.Cp=1100; %Heat capacity [Jkg^-1K^-1]   
p.Cc=62.7; %Heat capacity J/K (core) 
p.Cs=4.5; %Heat capacity J/K (surface)   
p.Rc=1.94;  % equiv. conduction resistance between core and surface
p.Ru=3.19;  % equiv. conduction resistance around the cell
p.h=10; % 5-10 for free convection air cooling %10-70 for forced air cooling
p.A=6.34e-3; %Cell area
p.M= 0.07;  % Cell mass [kg]

% Heat transfer parameters
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
p.C1 = 62.7;    % [J/K]
p.C2 = 4.5;     % [J/K]
p.h12 = 1.94; %1.9386; % [W/K]
p.h2a = 3.19;  % [W/K]

% Activation Energies,[J/mol]
% Taken from Prada et al (2012)
% https://doi.org/10.1149/2.064209jes

p.E.Dsn = 35e3;
p.E.Dsp = 35e3;
p.E.De = 26600;
p.E.kappa_e = 11e3;

% Taken from Ye, Y., Shi, Y., & Tay, A. A. O. (2012).
% https://doi.org/10.1016/j.jpowsour.2012.06.055

p.E.kn= 20000;
p.E.kp= 30000;


%% Aging submodel parameters

p.kappa_P = 1;      % [S/m] conductivity of side rxn product
