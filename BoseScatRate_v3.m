function [out] = BoseScatRate_v3(TTC, kappa)
% TTC is T/T_c
% kappa is the dimsionless momentum
% N is the total number of atoms
% Using the explicit integral equation for the fugacity
% fun1 = @(a, y, TTC, kappa) TTC.^3./(1 + exp(a + y.^2 - (1 - pi.^2./3.*TTC.^2)./TTC)).*a.^(3/2)./( 1 + exp(-(1 - pi.^2./3.*TTC.^2)./TTC + a + (y + kappa./sqrt(TTC)).^2));
% v2 - using harmonic density of states for integral, to check if we get the same result
% v3 - including missing number of condensate atoms

% Exact equation for z:
z = fugacityBose(TTC);
fun1 = @(a, y, TTC, kappa, z) TTC.^3 .* a.^(3/2) ./ ( z.^(-1).*exp( a + y.^2) - 1)  .* 1./( z.^(-1).*exp( a + (y + kappa./sqrt(TTC)).^2 ) - 1);

% q = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax)
int = zeros(size(TTC));

for i = 1 : length(TTC)
    int(i) = integral2(@(a,y)fun1(a, y, TTC(i), kappa, z(i)), 0,inf ,-inf,inf,'AbsTol', 1e-8,'RelTol',1e-8);
end

g31 = 1.202056903159595;
hbar = 1.054571628e-34;
N0overN = 1-TTC.^3; % Condensate fraction, N0/N
N0overN(TTC>1) = 0; % Avoid negative values
out = 1 + 4./(3*pi*g31)*int + 2.*N0overN ./ (z.^(-1).*exp(kappa.^2./TTC) - 1);
disp(2.*N0overN ./ (z.^(-1).*exp(kappa.^2./TTC) - 1))
end
