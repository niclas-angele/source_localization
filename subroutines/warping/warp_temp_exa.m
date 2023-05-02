function [s_w, fs_w]=warp_temp_exa(s,fs,r,c)
%%% warp_temp_exa.m
%%% Julien Bonnel, Woods Hole Oceanographic Institution
%%% March 2019

%%% Warping main function
%%% Inputs
%%%   s: original signal (time domain)
%%%   fs: sampling frequency
%%%   r,c: range and water sound speed (warping parameters)
%%% Outputs
%%%   s_w: warped signal (time domain)
%%%   fs_w: sampling frequency of the warped signal

if iscolumn(s)
    s=s.';
end

%% Step 1: preliminary computations
s=real(s);
N=length(s);
dt=1/fs;

tmin=r/c+dt;
tmax=N/fs+r/c;

%% Step 2: new time step
dt_w=iwarp_t(tmax,r,c)-iwarp_t(tmax-dt,r,c);

%% Step 3: new sampling frequency
fs_w=2/dt_w;

%% Step 4: new number of points
t_w_max=iwarp_t(tmax,r,c);
M=ceil((t_w_max)*fs_w);


%% Step 5: warped signal computation

% Warped time axis, uniform sampling
t_w=(0:M-1)/fs_w;

% Warped time axis, non-uniform sampling (starts from r/c)
t_ww=warp_t(t_w,r,c);

% factor for energy conservation
coeff=sqrt(t_w./t_ww);  

% Start exact interpolation (Shannon)
s_aux=repmat(s,M,1);
aux1=repmat(t_ww',1,N); % size=(M,1) -> repmat(1,N)
aux2=repmat(tmin+(0:N-1)/fs,M,1);       % size=(1,N) -> repmat(M,1)
aux=sinc(fs*(aux1-aux2));           % size=(M,N)
% end of exact interpolation --> interpolated signal is sum(s_aux.*aux,2)

% Final warped signal
s_w=coeff'.*sum(s_aux.*aux,2);