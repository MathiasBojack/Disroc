%%-------------------------------------------------------------------------
%
%         This function calculates the poroelastic parameters according to
%         the Biot and Coussy theory
%
%%-------------------------------------------------------------------------
function PoroProperty = PoroElasPara()

%%-------------------------------------------------------------------------
%
%                   INPTUT: Mechanical parameters
%   E:          Young's modulus                         [Pa]
%   nu:         Possoin ratio                           [1]
%   Ks:         Bulk modulus of the skeleton            [Pa]
%   Kf:         Bulk modlus of the fluid                [Pa]
%   phi:        porosity                                [1]
%
%--------------------------------------------------------------------------

E   = 1e10;
nu  = 0.25;
Ks  = 4e10;
Kf  = 2e9; % water
% Kf  = 1e8; % CO2
phi = 0.1;

%%-------------------------------------------------------------------------
%
%                   INPUT: Hydraulic parameters
%   
%   k:          The intrinsic permeability              [m2] = 1e12 Darcy
%               100mD = 1e-13m2 -> average permeability
%   mu:         viscosity                               [Pa.s]
%               1cP = 1e-3 Pa.s; 1P = 0.1Pa.s
%%-------------------------------------------------------------------------

k   = 1e-11;  
mu  = 5e-4; % water
% mu  = 5e-5; % CO2  

%--------------------------------------------------------------------------
%                   OUTPUT mechanical parameters
%
%   G:          Shear modulus                           [Pa]
%   K:          Bulk modulus                            [Pa]
%   CM=1/M:     Storage coefficient                     [Pa^-1]
%   b:          biot coeff                              [1]
%   Ku:         Undrained properties                    [Pa]
%   Kv = K + 4/3*G Vertical incompressibility              [Pa]
%%-------------------------------------------------------------------------

K   =    (E/3/(1-2*nu));
b   =    (1- K/Ks);
CM  =    (phi/Kf +(b-phi)/Ks); 
M   =    (1/CM);
Ku  =    (M*b^2+K);
G   =    (E/2/(1+nu)); 
Kv  =    (K + 4/3*G);

B    = b*M/Ku;
vu   = (3*nu+b*B*(1-2*nu))/(3-b*B*(1-2*nu)); % undrained poisson ratio
Kvu  = Ku*3*(1-vu)/(1+vu);
%--------------------------------------------------------------------------
%                        OUTPUT: Hydraulic parameters
%   Cf = 1/Kf:      Compressibility of the fluid            [Pa^-1]
%   c:              Hydraulic Diffussivity in 1D case       [m2.s^-1]
%   gamma:          Loading coefficient                     [1]
%   
%           dp/dt - c*d2p/dz2 = -gamma * dSzz/dt 
%%-------------------------------------------------------------------------
Cf  = 1/Kf;

c   =  (k/mu * M*K/Ku   * (1- 4/3*G /Kv * M*b^2/Ku)^-1);
gamma =  (K/Kv * M*b/Ku * (1- 4/3*G /Kv * M*b^2/Ku)^-1);

PoroProperty.c      = c;
PoroProperty.gamma  = gamma;


PoroProperty.K      = K;
PoroProperty.b      = b; 
PoroProperty.CM     = CM;
PoroProperty.M      = M;
PoroProperty.Ku     = Ku;
PoroProperty.G      = G; 
PoroProperty.Kv     = Kv;
PoroProperty.k      = k;
PoroProperty.mu     = mu;
PoroProperty.B      = B;
PoroProperty.vu     = vu;
PoroProperty.Kvu   = Kvu;

end


