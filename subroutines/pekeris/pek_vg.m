%%%%% Warping tutorial
%%%%% Julien Bonnel
%%%%% January 2020

function [vg,f_]=pek_vg(fmin,fmax,m_min,m_max,c1,c2,rho1,rho2,D,df)
% Compute group speeds in a Pekeris waveguide
% Inputs :
    % c1 / rho1 : sound speed / density in water
    % c2 / rho2 : sound speed / density in the seabed
    % D : water depth
    % freq : frequency
% Outputs :
    % vg : matrix of group speed (size Nmode * Nfreq)
    % f_ : frequency axis

f_=fmin:df:fmax;
Nf=length(f_);
A=rho2/rho1;

vg=zeros(m_max-m_min+1,Nf);

ff=1;
for f=f_    
    w=2*pi*f;
    % Horizontal wavenumbers
    kr = pek_kr_f(c1,c2,rho1,rho2,D,f);
    
    % No need to start computation if there is no propagating mode
    if ~isempty(kr)
        m_max_f=length(kr);
        if m_min <= m_max_f
            m_max_calc=min(m_max_f,m_max);
            % in a Pekeris waveguide, we have an exact formula to compute
            % group speed from  wavenumber
            ind_m=1;  
            for m=m_min:m_max_calc                
                k=kr(m);
                % vertical wavenumber in water
                g1=sqrt(w.^2/c1^2-k.^2);
                % vertical wavenumber in seabed
                g2=sqrt(k.^2-w.^2/c2^2);
                
                % group speed computation (no need to numerically
                % approximate dk/df)
                a=A./(g2.^2+A^2*g1.^2);
                vg(ind_m,ff)=k./w*c1^2*c2^2.*(g2*D+a.*(g1.^2+g2.^2))./(c2^2*(g2*D+a.*g2.^2)+c1^2*a.*g1.^2);
                ind_m=ind_m+1;
            end
            
        end        
    end

    ff=ff+1;
end

vg(vg==0)=NaN;

