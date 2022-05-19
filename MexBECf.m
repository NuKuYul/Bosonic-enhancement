function [f] = MexBECf(TTC,Vcontrol)
% Return the BEC fraction for different T/T_c and ODT control voltage
switch Vcontrol
    case 0.1
        
        load('f0.1V.mat', 'res')
        f=res(TTC);
   
    case 0.15
        load('f0.15V.mat', 'res')
        f=res(TTC);
   
    case 0.2
       load('f0.2V.mat', 'res')
        f=res(TTC);
   
        
    case 0.25
        load('f0.25V.mat', 'res')
        f=res(TTC);
   
    otherwise
        msg = 'Control voltage not allowed';
        error(msg)
end

end