function s_r=iwarp_temp_exa(s_w,fs_w,r,c,fs,N)
%%% iwarp_temp_exa.m
%%% Julien Bonnel, Woods Hole Oceanographic Institution
%%% March 2019

%%% Inverse warping main function
%%% Inputs
%%%   s_w: warped signal (time domain)
%%%   fs_w: sampling frequency of the warped signal
%%%   r,c: range and water sound speed (warping parameters)
%%%   fs: sampling frequency of the original/unwarped signals
%%%   N: number of samples of the original/unwarped signals
%%% Outputs
%%%   s_r: signal after inverse warping

%% Preliminary steps
M=length(s_w);

if iscolumn(s_w)
    s_w=s_w';
end

%% Inverse warping computation

% Time axis, uniform sampling (starts from r/c+dt)
t=(1:N)/fs+r/c;

% Time axis, non-uniform sampling
t_iw=iwarp_t(t,r,c);

% factor for energy conservation
coeff=sqrt(t./t_iw);  

% Start exact interpolation (Shannon)
s_aux=repmat(s_w,N,1);  % initial signal replicated N times (N rows)
aux1=repmat(t_iw',1,M);
aux2=repmat((0:M-1)/fs_w,N,1);
aux=sinc(fs_w*(aux1-aux2));
% end of exact interpolation --> interpolated signal is sum(s_aux.*aux,2)

% Final warped signal
s_r=real(coeff'.*sum(s_aux.*aux,2));
      
