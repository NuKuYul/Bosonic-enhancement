function plot_BoseScatRate_MF

tic

% Define the value for physical constants
kB = 1.3806504e-23; %J/K
hbar = 1.054571628e-34;
m = 22.9897692807/6.02214086e26;
g = 2*pi*9.7946e6;
lambdaAtom = 589.1583264e-9; %Sodium D2 line
c = 2.99792458e8; % speed of light [m/s]
N=0.4*10^6;
% a0=5.29*10^(-11);
% a=85*a0;
Erec=2.39*10^(-6)*kB;

% Vcontrol is the control voltage for the ODT
Vcontrol=[0.1, 0.15, 0.2, 0.25];
% Vcontrol = 0.175;

%calculating the parameters of ODT
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
% Omega0 = (omegaz .* omegar .* omegar).^(1/3);

TC = hbar .* (omegaz .* omegar .* omegar .* N / g31).^(1/3)./kB; % Critical temperature for BEC - harmonic [K]
% Experimentally, the value of harmonic Tc is quite close to the Tc we
% observed. In the semi-ideal gas approximation, the backaction from
% thermal clud to BEC is ignored. Therefore, Tc in this case is exactly the
% same as our conventional result.

% TTC = [linspace(0.02, 0.8, 150) linspace(0.8, 1.2, 200) linspace(1.2, 2, 150)];
kappa = sqrt(Erec./TC./kB); % Dimensionless recoil momentum




TTC = [linspace(0.15, 0.99, 500)];
% TTC = [linspace(0.01, 0.9, 90) linspace(0.9,1.1, 500) linspace(1.1, 1.5, 40)];
% Vcontrol=[0.1, 0.15, 0.2, 0.25];
% specific experimental values:
% N = 0.5e6;
% Vcontrol = 0.3;
% ODT = calcXODTSodium(N, Vcontrol);
% T = TTC.*ODT.Tc;
% Trec = 2.3998e-6; % Sodium recoil temperautre, in uK
% alpha = pi/2; %scattering angle
% kappa = 2.*sqrt(Trec./ODT.Tc).*(sin(alpha./2)); %angle dependent kappa, for alpha=90 degrees its the same as kappa = sqrt(Trec./Tf). Corrected for the typo of squre on the sine function from Thywissen's paper
% folder = [fileparts(mfilename('fullpath')) '\'];
% figure('Filename', [folder 'Bose enhancement vs TTC.fig']);
% hold on;
% BECfrac = 1-(TTC).^3;
% kappa = [0.47 0.5  0.53 0.6];
Gammaall=zeros(size(TTC));
Gammathermal=zeros(size(TTC));
GammaBEC=zeros(size(TTC));

%Calculate and plot the reults
for i0=1:length(kappa)
    parfor i = 1 :length(TTC)
      [Gammaall(i),Gammathermal(i),GammaBEC(i) ] = BoseScatRate_MF(TTC(i),kappa(i0),Vcontrol(i0)); 
      % Use BoseScatRate_MF for interacting theory in harmonic trap, use BoseScatRate_v3 for
      % noninteracting theory in harmonic trap, use BoseScatRateBox for
      % noninteracting theory in box trap.
      
%     Gamma(i) = BoseScatRate_v3(TTC(i),kappa);
%     plot(BECfrac, Gamma*2.2e-4)
        
    end
    figure;
    hold on;
    plot(TTC, Gammaall,'-')
    plot(TTC, Gammathermal,'-')
    plot(TTC, GammaBEC,'-')
    legend(['kappa = ' num2str(kappa(i0))])
    xlabel('T/T_C')
    % xlabel('BEC fraction')
    ylabel('Enhancement')
    hold off;
end

toc

%% Plot vs kappa
% kappa = linspace(0.5, 0.7, 30);
% Gamma = zeros(size(kappa));
% for i = 1 :length(kappa)
%     Gamma(i) = BoseScatRate_MF(0.63, kappa(i),Vcontrol);
% 
% end
% figure; %('Filename', [folder 'Bose enhancement vs kappa.fig']);
% plot(kappa, Gamma)
% xlabel('Kappa = sqrt[T_{rec}/T_c]2Sin(\alpha/2)')
% ylabel('Enhancement, T/T_c=1')
%% Check mirror symmetry
% figure;
% plot(TTC(TTC>1), Gamma(TTC>1))
% hold on;
% plot(1./TTC(TTC<1), Gamma(TTC<1))