%%%%% Warping tutorial
%%%%% Julien Bonnel
%%%%% January 2020


function [kr, kz, A2]=pek_init(rho1, rho2, c1, c2, D, freq)


Nf=length(freq);

% Looking at the highest frequency for the max number of modes
toto=pek_kr_f(c1,c2,rho1,rho2,D,max(freq));
Nm=length(toto);

% Init variables
kr=zeros(Nf,Nm);
kz=zeros(Nf,Nm);
kzb=zeros(Nf,Nm);
A2=zeros(Nf,Nm);

%%% Loop over frequencies
for ff=1:Nf
    
    toto=pek_kr_f(c1,c2,rho1,rho2,D,freq(ff));
    NN=length(toto);
    w=2*pi*freq(ff);
    
    % Horizontal wavenumber
    kr(ff,1:NN)=toto;
    
    % Vertical wavenumber in water
    kz(ff,1:NN)=sqrt(w^2/c1^2-toto.^2);
    
    % Vertical wavenumber in seabed
    kzb(ff,1:NN)=sqrt(toto.^2-w^2/c2^2);
    
    % Normalization for the modal depth function
    A2(ff,1:NN)=2*rho1*kz(ff,1:NN).*kzb(ff,1:NN)./(kz(ff,1:NN).*kzb(ff,1:NN)*D - 0.5*kzb(ff,1:NN).*sin(2*kz(ff,1:NN)*D) + rho1/rho2*kz(ff,1:NN).*(sin(kz(ff,1:NN)*D).^2));    
end
% Global multiplicative factor
A2=A2*1i*exp(1i*pi/4)/rho1/4;
