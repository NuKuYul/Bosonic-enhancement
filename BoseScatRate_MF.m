function [outall, outthermal,outBEC] = BoseScatRate_MF(TTC, kappa,Vcontrol)
% TTC is T/T_c
% kappa is the dimsionless momentum
% N is the total number of atoms
% Using the explicit integral equation for the fugacity

% outall is the total scattering rate, while outthermal includes only the
% thermal-thermal scattering, and outBEC includes only BEC-thermal
% scattering

% Define the physical constants
kB = 1.3806504e-23; %J/K
hbar = 1.054571628e-34;
m = 22.9897692807/6.02214086e26;
g = 2*pi*9.7946e6;
lambdaAtom = 589.1583264e-9; %Sodium D2 line
c = 2.99792458e8; % speed of light [m/s]
N=0.4*10^6;
a0=5.29*10^(-11); % Bohr radius
a=85*a0; % scattering length for F=2


% Calculating the parameters of the ODT
Trans = 0.9936 * 0.999* 0.9926 * 0.995 * 0.9975; %Transmission of two lenses (AC508-500-C-ML and AC508-250-C-ML), one laser line mirror, and two surfaces of the vacuum window (from Aviv's thesis, Appendix C, p. 125).
P = Trans.*(1.7893.*Vcontrol + 0.0014); % taking into account optical transmissions [W]
w0 = 7.9786e-6;
lambdaIR = 1064e-9;
omega0 = 2*pi*c/lambdaAtom;
omega = 2*pi*c/lambdaIR;
zR = pi*w0*w0/lambdaIR * 2684/2019; %the extra factor is needed in order to explain the low measured axial frequency (2020-02-27)
U0 = 3*pi*c^2/(2*omega0^3)*( g/(omega0-omega) + g/(omega0+omega))*2.*P./(pi*w0*w0);

omegar = sqrt(4.*U0./(m*w0*w0));
omegaz = sqrt(2.*U0./(m*zR^2));
g31 = 1.202056903159595;
Omega0 = (omegaz .* omegar .* omegar).^(1/3);


TC = hbar .* (omegaz .* omegar .* omegar .* N / g31).^(1/3)./kB; % Critical temperature for BEC - harmonic [K]
% Experimentally, the value of harmonic Tc is quite close to the Tc we
% observed.
% TC = T_from_N_ODT(N, Vcontrol); % Critical temperature for BEC - anharmonic [K]
T=TTC.*TC;
f =MexBECf(TTC,Vcontrol); % Condensate fraction, N0/N
f(TTC>1) = 0; % Avoid negative values
aho=sqrt(hbar/(m*Omega0)); % oscillator length
% ahoz=sqrt(hbar/(m*omegaz));
% ahor=sqrt(hbar/(m*omegar));
% Mu0=(hbar.*Omega0)/2.*(15.*a./aho.*N.*f).^(2/5); %%%%The chemical potential given by the interacting BEC

% calculate the chemical potential and fugacity
if TTC>1
    Mu0=0;
    zODT = fugacityBose(TTC);
else
    Mu0=Mex(TTC,Vcontrol).*1e-6.*kB;
    zODT= exp(Mu0./(kB.*T));
end

rr=sqrt((2*Mu0)/(m*omegar^2)); % radial zise of the cloud
rz=sqrt((2*Mu0)/(m*omegaz^2)); % axial size of the cloud
r0=sqrt((2*Mu0)/(m*Omega0^2)); % geometric mean of the size
lambdaT=hbar*2*pi/sqrt(2*pi*m*kB*T); % thermal wavelength

% calculating the trap averaged overlap between BEC and thermal cloud.
% the parameter x represents the normalized radius, x^2 comes from the
% integration measure, while 1-x^2 is proportional to the BEC density
funBEC1= @(x,TTC,kappa) x.^2.*(1-x.^2)./(-1+exp(kappa.^2./TTC+Mu0.*(1-x.^2)./(kB.*T)));
funBEC2= @(x,TTC,kappa) x.^2.*(1-x.^2);

% Calculate the trapping potential for the thermal cloud
if f>0
    UODT=@(r) 0.5.*m.*Omega0.^2.*r.^2+ 2.*Mu0.*max((1-(r./r0).^2),0);
else
    UODT=@(r) 0.5.*m.*Omega0.^2.*r.^2;
end

int = zeros(size(TTC));

%calculate the light scattering within the thermal cloud
for i = 1 : length(TTC)
  intergrandODT = @(a, y, r) r.^2 ./ ( zODT.^(-1).*exp( a + y.^2 + 1./(kB*T).*UODT(r)) - 1 )  .* 1./( zODT.^(-1).*exp( a + (y + kappa./sqrt(TTC)).^2 + 1./(kB*T).*UODT(r) ) - 1);
   IntSizeR = 1.39857*10;% it was 14*10
   epsilon=1e-6*r0; % cutoff to prevent singularity
%    IntSizeZ = 7;
   intODT1 = integral3(@(a, y, r) intergrandODT(a, y, r), 0, Inf, -Inf, Inf, 0,r0-epsilon,'AbsTol', 1e-3,'RelTol',1e-10, 'method', 'iterated');
   intODT2 = integral3(@(a, y, r) intergrandODT(a, y, r), 0, Inf, -Inf, Inf, r0+epsilon,IntSizeR*w0,'AbsTol', 1e-3,'RelTol',1e-10, 'method', 'iterated');
%    intODT2 = integral(intODT1, -IntSizeZ*zR,  IntSizeZ*zR, 'ArrayValued', true, 'AbsTol', 1e-3,'RelTol',1e-10);

% the factor (1-8*sqrt(2)*a/lambdaT) took care of the modified pair
% correlation function
    int(i) = 1 + (1-8*sqrt(2)*a/lambdaT).*(2*m*kB*T)^(3/2)*(2*pi)^2./( (2*pi*hbar)^3.*N ) .* (intODT1+intODT2);
%     int(i) = integral2(@(a,y)fun1(a, y, TTC(i), kappa, zODT(i)), 0,inf ,-inf,inf,'AbsTol', 1e-8,'RelTol',1e-8);
end

intBEC = zeros(size(TTC));

%calculate the light scattering between BEC and thermal cloud
for i = 1 : length(TTC)
%     test=integral(@(x)funBEC1(x, TTC(i), kappa), 0,1 ,'AbsTol', 1e-8,'RelTol',1e-8);
%     test2=integral(@(x)funBEC2(x, TTC(i), kappa), 0,1 ,'AbsTol', 1e-8,'RelTol',1e-8);
    intBEC(i) = integral(@(x)funBEC1(x, TTC(i), kappa), 0,1 ,'AbsTol', 1e-8,'RelTol',1e-8)./integral(@(x)funBEC2(x, TTC(i), kappa), 0,1 ,'AbsTol', 1e-8,'RelTol',1e-8);
end


% out = 1 + 4./(3*pi*g31)*int + 2.*f .*intBEC;

% the factor (1-8*a/lambdaT) took care of the modified pair
% correlation function
outthermal = int;
outBEC = (1-8*a/lambdaT).* 2.*f .*intBEC;
outall = int + outBEC;
end

