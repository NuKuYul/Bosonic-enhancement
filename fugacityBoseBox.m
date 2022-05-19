function [z, q] = fugacityBoseBox(TTC)
% Fugacity of a Bose gas in an box trap, as a function of TTC = T/T_c
% Equation from Pathria, equation 7.26

syms z
fugacity = zeros(size(TTC));
for i = 1:length(TTC)
    if TTC(i)<=1 % Below one - set to z=1 (assumes no interactions)
        fugacity(i) = 1;
    elseif TTC(i)>40 % Use classical expression (numeric cause noise, difference is smaller than 1e-3).
        xi32 = 2.612375348685488;
        fugacity(i) = xi32./(TTC(i).^(3/2));
    else
        fugacity(i) = fminsearch( @(z)fugacityEquation(TTC(i),z), 1);
    end
end

% for i = 1:length(TTC)
%     if TTC(i)<=1 % Below one - set to z=1 (assumes no interactions)
%         fugacity(i) = 1;
% elseif TTC(i)>5 % T/T_c>5 - use classical expression (numeric cause noise, difference is smaller than 0.05%).
%     g31 = 1.202056903159595;
%     fugacity(i) = g31./TTC(i).^3;
%     else % numeric function
%         fugacity(i) = fminsearch( @(z)fugacityEquation(TTC(i),z), 1);
%     end
% end

clear z
z = fugacity;
q = log(fugacity);

end

function eqDiff = fugacityEquation(TTC,z) %Box fugacity equation, equation 7.26 from Pathria
xi32 = 2.612375348685488; %xi32 = zeta(3/2) = polylog(3/2,1) = BoseFunc(3/2,1);
eqDiff = abs(polylog(3/2, z) - xi32./( TTC.^(3/2) ) );
end


%% Compare box fugacity calcualtions, Bose gas: exact, and high tmperature limit
% clear;
% TTC = [ 0.05 linspace(1,50,40)];
% [~, qHarmonic] = fugacityBoseBox(TTC); % exact expression for box potential
% mukBTc = qHarmonic.*TTC; %mu/Ef, exact expression for harmonic potential
% 
% mukBTcHighT = -TTC.*log(TTC.^(3/2) ./ zeta(3/2) ) ; %mu/EF, box potential, high temperature limit
% 
% figure;
% plot(TTC, mukBTc,'b')
% hold on;
% plot(TTC, mukBTcHighT, 'b--')
% xlabel('T/T_c')
% ylabel('\mu/K_BT_c')
% 
% legend({'Exact','High temperature limit'}, 'Location', 'Best')
