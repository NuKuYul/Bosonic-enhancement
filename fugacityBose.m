function [z, q] = fugacityBose(TTC)
% Calculate the fugacity z=exp(mu*beta) and the logarithmic fugacity q=mu*beta of a Bose gas as a function of T/TC = TTC, where TC is the BEC critical temperature
% Assumes an harmonic potential. Equations hold only above TC
% For T/T_1<1, we set mu=0 => z=1.
syms z
fugacity = zeros(size(TTC));
for i = 1:length(TTC)
    if TTC(i)<=1 % Below one - set to z=1 (assumes no interactions)
        fugacity(i) = 1;
    elseif TTC(i)>5 % T/T_c>5 - use classical expression (the numeric solution cause noise, difference is smaller than 0.05%).
        g31 = 1.202056903159595;
        fugacity(i) = g31./TTC(i).^3;
    else % numeric function
        fugacity(i) = fminsearch( @(z)fugacityEquation(TTC(i),z), 1);
    end
end
clear z
z = fugacity;
q = log(fugacity);

end

function eqDiff = fugacityEquation(TTC,z) % Harmonic fugacity equation, Eq. 24 from Varenna notes, for bosons
persistent polylog3Spline
%spline of the polylog function to speedup computation, in the range -2017:1
if isempty(polylog3Spline)
    polylog3Spline = load('polylog3Spline.mat', 'polylogSpline');
end
g31 = 1.202056903159595; %g31 = polylog(3,1) = BoseFunc(3,1);
eqDiff = abs(polylog3Spline.polylogSpline(z) - g31./(TTC).^3);
end

%% Test the code:
%% Compare harmonic fugacity calcualtions, Bose gas: exact, Sommerfeld, and high tmperature limit
% clear;
% TTC = [ 0.05 linspace(0.1,25,25)];
% [~, qHarmonic] = fugacityBose(TTC); % exact expression for harmonic potential
% mukBTc = qHarmonic.*TTC; %mu/Ef, exact expression for harmonic potential
% % muEfSommerfeld = 1 - pi.^2./3.*TTC.^2 ; %mu/EF, Sommefeld expansion, harmonic
% mukBTcHighT = -TTC.*log(1./BoseFunc(3,1).*TTC.^3) ; %mu/EF, harmonic, high temperature limit
% 
% figure;
% plot(TTC, mukBTc,'b')
% hold on;
% % plot(TTC, muEfSommerfeld, 'b--')
% plot(TTC, mukBTcHighT, 'b--')
% xlabel('T/T_c')
% ylabel('\mu/K_BT_c')
% 
% legend({'Exact','High temperature limit'}, 'Location', 'Best')