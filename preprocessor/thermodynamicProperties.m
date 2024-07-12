function [rhoG_phy, rhoG_lu, lamb] = thermodynamicProperties(Pi_phy, Ti_phy)
%% Simulation conditions and critical points for methane

% Component: 1-CH4   2-CO2
Component = 1;          % overall mole fraction of CH4
if Component == 1
    k = 1;
elseif zli ==0
    k = 2;
end

% Bulk critical property of components
Tc_phy = [190.564;304.13];          % bulk critical temperature, K;
Pc_phy = [4.5992;7.377]*10^6;       % bulk critical pressure, Pa;
rhoc_phy = [162.66;467.6];          % bulk critical density, kg/m3;
w_phy = [0.011;0.239];              % http://www.kaylaiacovino.com/Petrology_Tools/Critical_Constants_and_Acentric_Factors.htm
Mw_phy = [16.043;44.0095]/1000;     % kg/mol
d_mole_phy = [0.38, 0.33]*10^(-9);      % diameters of CH4 and CO2 molecules, m

% Critical temperature, pressure and molar volume, and acentric factors    
Tc_component = Tc_phy(k);       % K
Pc_component = Pc_phy(k);       % Pa
w_component = w_phy(k);
Mw = Mw_phy(k);                 % kg/mol
d_mole = d_mole_phy(k); 

% PR EOS parameters for components 
R = 8.31446; % M3*Pa/K/mol
k_component = 0.37464+1.54226*w_component-0.26992*w_component.^2;
alpha_component = (1+k_component.*(1-sqrt(Ti_phy./Tc_component))).^2;
a_component = 0.45724*(R*Tc_component).^2./Pc_component;
a_alpha_component = a_component.*alpha_component;
b_component = 0.07780*R*Tc_component./Pc_component;


aG = a_alpha_component;
bG = b_component;

% cubic equation of compressibility (Z) ← Pv=RT (n=1) & PR EOS
AG = aG*Pi_phy/R^2/Ti_phy^2;
BG = bG*Pi_phy/R/Ti_phy;

% Z3*Z^3 + Z2*Z^2 + Z1*Z + Z0 = 0
ZG3 = 1;
ZG2 = BG-1;
ZG1 = AG-3*BG^2-2*BG;
ZG0 = BG^3+BG^2-AG*BG;
ZG_results = roots([ZG3 ZG2 ZG1 ZG0]);
ZG_real = ZG_results(ZG_results==real(ZG_results));
ZG = max(ZG_real);
vG = ZG*R*Ti_phy/Pi_phy;

rhoG_phy_mol = 1/vG;     % mol/m^3 - molar density
 
Rs = R./Mw;
rhoG_phy = rhoG_phy_mol*Mw;     % kg/m^3


%% mfp - mean free path
b_Rs_component = 0.07780*Rs.*Tc_component./Pc_component;                 % m3/kg
bs_lu = 2/21;
Mw_lu = 1.0;

C_bs = b_Rs_component/bs_lu;  % (m3/kg)/(m_lu3/kg_lu)             = Cl^3/Cm       = [ 3, 0,-1, 0, 0]
Crho = 1/C_bs;

% 5 unknown conversion factors and 5 equations                                          % [Cl,Ct,Cm,CT,Cn]
C_bs = b_Rs_component/bs_lu;  % (m3/kg)/(m_lu3/kg_lu)             = Cl^3/Cm       = [ 3, 0,-1, 0, 0]
C_Mw = Mw/Mw_lu;              % (kg/mol)/(kg_lu/mol_lu)           = Cm/Cn         = [ 0, 0, 1, 0,-1]

Crho_mol = 1 / (C_bs * C_Mw);         % conversion factor of mole density, mol/m3         = [-3, 0, 0, 0, 1]

% Density of components ( mass of component / volume of mixture) -- cell-volume densities,
rhoG_lu = rhoG_phy/Crho;            % kg/m^3

% 
r_mole = d_mole/2;

rho_mol = rhoG_lu/Mw_lu;
n = rho_mol*Crho_mol * 6.02214076*10^23;
lamb = 1./(sqrt(2)*n*pi*(2*r_mole+2*r_mole)^2);
end