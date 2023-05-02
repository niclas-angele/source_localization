%%%%% Warping tutorial
%%%%% Julien Bonnel
%%%%% January 2020


function [g]=pek_green(kr, kz, A2, zs, zr, r)

Nz=length(zr);
[Nf Nm]=size(kr);

g=zeros(Nf,Nz);

for ff=1:Nf
    toto=find(kr(ff,:)~=0);
    g(ff,:)=A2(ff,toto).*sin(kz(ff,toto)*zs).*exp(-1i*kr(ff,toto)*r)./sqrt(kr(ff,toto)*r)*sin(kz(ff,toto).'*zr); 
end
