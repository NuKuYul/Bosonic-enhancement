function [out] = BoseScatRateBox(TTC, kappa)
% Same calculation as in BoseScatRate, but here we are calcualting for a box potential (and not harmonic trap).
% For a box, the power of coordiante 'a' changes from a^3/2 to a^0=1
% TTC is T/T_c
% kappa is the dimsionless momentum

z = fugacityBoseBox(TTC); % Exact expression for z

fun1 = @(a, y, TTC, kappa, z) TTC.^(3/2) .* 1 ./ ( z.^(-1).*exp( a + y.^2) - 1)  .* 1./( z.^(-1).*exp( a + (y + kappa./sqrt(TTC)).^2 ) - 1);

% q = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax)
int = zeros(size(TTC));

for i = 1 : size(TTC,1)
    for j = 1 : size(TTC,2)
        int(i,j) = integral2(@(a,y)fun1(a, y, TTC(i,j), kappa, z(i,j)), 0,inf ,-inf,inf,'AbsTol', 1e-10,'RelTol',1e-10, 'method', 'iterated');
    end
end
xi32 = 2.612375348685488; %xi32 = zeta(3/2) = polylog(3/2,1) = BoseFunc(3/2,1); 
out = 1 + int ./ (sqrt(pi) .* xi32);

end
