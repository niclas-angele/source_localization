%%%%% Warping tutorial
%%%%% Julien Bonnel
%%%%% January 2020


function rho = hamilton(Cp)
Cp = Cp/1000;
if ((Cp>=1.55)&&(Cp<1.9))
    rho = 1.135*Cp-0.190;
elseif((Cp>=1.90) && (Cp<3.0))
    rho = 2.351-7.497*Cp^-4.656;
elseif((Cp>=3.0)&&(Cp<6.60))
    rho = 1.979*Cp+0.112;
else
    warning('out of hamilton bounds');
    rho = 1.5;
end
