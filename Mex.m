function [mu] = Mex(TTC,Vcontrol)
% Return the chemical potential for different T/T_c and ODT control voltage
switch Vcontrol
    case 0.1
        
        load('mu0.1V.mat', 'res')
        mu=res(TTC);
   
    case 0.15
        load('mu0.15V.mat', 'res')
        mu=res(TTC);
   
    case 0.2
       load('mu0.2V.mat', 'res')
        mu=res(TTC);
   
        
    case 0.25
        load('mu0.25V.mat', 'res')
        mu=res(TTC);
   
    otherwise
        msg = 'Control voltage not allowed';
        error(msg)
end

end