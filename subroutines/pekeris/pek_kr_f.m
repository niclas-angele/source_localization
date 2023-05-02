%%%%% Warping tutorial
%%%%% Julien Bonnel
%%%%% January 2020


function [ root_ok_pek ] = pek_kr_f(c1,c2,rho1,rho2,D,freq)
% Compute horizontal wavenumbers in a Pekeris waveguide
% Inputs :
    % c1 / rho1 : sound speed / density in water
    % c2 / rho2 : sound speed / density in the seabed
    % D : water depth
    % freq : frequency
% Outputs :
    % root_ok_pek : vector of wavenumbers 
    % (size = number of modes at the considered frequency)


jj=sqrt(-1);

% wavenumbers are between k1 and k2
w=2*pi*freq;
k1=w/c1;
k2=w/c2;

% Search space for wavenumber
n_test=5000;
kr=linspace(k2*(1+10^-15),k1*(1-10^-15),n_test);

% Pekeris waveguide equation
func_pek=@(krm) tan(sqrt(k1^2-krm.^2)*D)+rho2*sqrt(k1^2-krm.^2)/rho1./(jj*sqrt(krm.^2-k2^2)).*jj;

% Let's work in small interval to help solving the equation 
toto_pek=tan(sqrt(k1^2-kr.^2)*D)+rho2*sqrt(k1^2-kr.^2)/rho1./(jj*sqrt(kr.^2-k2^2)).*jj;
j=1;
guess_pek=[];
for i=1:n_test-1
    if toto_pek(i)/toto_pek(i+1)<=0
        guess_pek(j,1:2)=[kr(i) kr(i+1)];
        j=j+1;
    end
end

root_ok_pek=[];

% Equation is solved with fzero
if ~isempty(guess_pek)
    for j=1:size(guess_pek,1)
        root_pek(j)=fzero(func_pek,guess_pek(j,1:2));
    end


    j=1;
    for i=1:length(root_pek)
        if abs(func_pek(root_pek(i)))<1 && isreal(root_pek(i))
            root_ok_pek(j)=root_pek(i);
            j=j+1;
        end
    end
    
end

% We order the wavenumbers
root_ok_pek=root_ok_pek(end:-1:1);
