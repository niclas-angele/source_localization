%%%%% ﻿Automated approach for source location in shallow waters
%%%%% Angele Niclas
%%%%% 2023


clear; close all;
addpath(genpath('./subroutines'))
load pekeris_model_tool.mat
% s_t: propagated modes in a Pekeris waveguide with parameters
%       c1, c2, rho1, rho2: sound speed / density
%       D: depth
%       r: range
%       zs, zr: source/receiver depth
% mode_t: tabular containing in each column the modal components of the
%       signal s_t
% fs: sampling frequency


% %% Artifical noise on data
%signal length
N=length(s_t);
%parameters of the noise 
Tdelta=0.01; 
delta=0.001; 
% creation of the Gaussian process W
time=(-N:N-1)/fs;
Y=randn(1,2*N); 
C=delta^2*exp(-1/2/Tdelta^2*time.^2);
filter=fft(fftshift(C));
W=real(ifft( sqrt(filter).*fft(Y) ));
W=W(N/2:3*N/2-1).';
s_t=s_t+W.';


%% Filtering of modes 
%time vector
time=(0:N-1)/fs;
%FFT parameters
NFFT=2048;
freq=(0:NFFT-1)*fs/NFFT;
%filtering of each mode (the detailed function can be found below).
%mode_t_app countains in each column the approximation of the modal
%components
mode_t_app=filtering(s_t,time,fs,NFFT,4);


%% Computation of the dispersion curves
%value of sigma
sig=5*fs/100; 
%length of the associated slicing window
NW=unique(2*floor((5*fs./sig+1)/2)+1); 

%spectrogram of mode 1
spectro_1=abs(tfrstft(mode_t_app(:,1),1:N,NFFT,gausswin(NW))).^2;
%diversion curves tm_1 and associated frequencies freq_1
[freq_1,tm_1]=meth_max(spectro_1,NFFT,time,fs);
%amplitude of each frequency 
amp_1=max(spectro_1(1:NFFT/2-1,:).');
%same procedure for modes 2,3,4
spectro_2=abs(tfrstft(mode_t_app(:,2),1:N,NFFT,gausswin(NW))).^2;
[freq_2,tm_2]=meth_max(spectro_2,NFFT,time,fs);
amp_2=max(spectro_2(1:NFFT/2-1,:).');
spectro_3=abs(tfrstft(mode_t_app(:,3),1:N,NFFT,gausswin(NW))).^2;
[freq_3,tm_3]=meth_max(spectro_3,NFFT,time,fs);
amp_3=max(spectro_3(1:NFFT/2-1,:).');
spectro_4=abs(tfrstft(mode_t_app(:,4),1:N,NFFT,gausswin(NW))).^2;
[freq_4,tm_4]=meth_max(spectro_4,NFFT,time,fs);
amp_4=max(spectro_4(1:NFFT/2-1,:).');

%vector containing all the dispersion curves 
T=[tm_1,tm_2, tm_3, tm_4].';
%computation of max(S)
M=max([amp_1,amp_2,amp_3,amp_4]);
%threshold 
p=0.4; 
%vector containing 1 if we keep the part of the curve and 0 otherwise
AMP=[1.*(amp_1>M*p);1.*(amp_2>M*p);1.*(amp_3>M*p);1.*(amp_4>M*p)]; 

%%% plot %%%
figure
hold on 
for i=1:4
    plot(AMP(i,:).*T(i,:),freq_1,'bo')
end
vg=pek_vgb(freq,1,4,1500,1600,1000,1500,100); 
plot(10000./vg, freq,'k')
xlim_plots=[6.6 7];
xlim(xlim_plots)
ylim([0,100])
xlabel('Times (s)')
ylabel('Frequency (Hz)')
title('Values of tnapp and estimated dispersion curves')
pause(.1)
%%% end plot %%%


%% Minization of the penalized problem
%initial condition
x0=[r,1500,1600,1000,1500,100,0]; 
%functional J, detailed below
g=@(x) J(x(1),x(2),x(3),x(4),x(5),x(6),x(7),freq_1,T,AMP); 
%penalization
alpha=g(x0); 
g=@(x) g(x)+(norm((x(6)-100)/100)^2+10*norm((1500-x(2))/1500)^2+10*norm((x(4)-1000)/1000)^2+norm((x(5)-1500)/1500)^2+10*((x(1)<9000)+(x(1)>11000))+10*((x(3)<1550)+(x(3)>2000))+10*((x(7)>2)+(x(7)<-2)))*0.001; 
%minization 
options = optimset('Display','none'); 
x=fminsearch(g,[r,1500,1600,1000,1500,100,0],options); 

fprintf('rapp=%dm\n c1app=%dm/s\n c2app=%dm/s\n rho1app=%dkg/m3\n rho2app=%dkg/m3\n Dapp=%dm\n dtapp=%ds\n',x(1),x(2),x(3),x(4),x(5),x(6),x(7))

%%% plot %%%
vg=pek_vgb(freq,1,4,x(2),x(3),x(4),x(5),x(6)); 
plot(x(1)./vg-x(7), freq,'r')
legend('tnapp','','','','exact disp. curves','','','','approx. disp. curves')
%%% end plot %%%


%% Detailed filtering function 
function mode_t_app=filtering(s_t,time,fs,NFFT,nbmode)
    %list to contain all the extracted modal components
    mode_t_app=[]; 
    %cut the signal so keep only the significative part
    ind1=1; 
    ind2=length(s_t);
    m=max(abs(s_t));
    while s_t(ind1)<0.01*m
        ind1=ind1+1; 
    end
    while s_t(ind2)<0.01*m
        ind2=ind2-1; 
    end
    %copy of the received signal
    str=s_t; 
    for nb=1:(nbmode-1)
        %list with all the quality factors
        Q=[];
        %list of positions in time to test t0
        L=1320:1370;
        %test to find the best t0
        for i=1:length(L)
            ind0=L(i); 
            %signal cut after t0
            stbb=str(ind0:ind2);
            %warping associated to t0
            [s_w, ~]=warp_temp_exa(stbb,fs,time(ind0),1);
            %warped spectrogram using a small value of sigma
            spectro_w=abs(tfrstft(s_w,1:length(s_w),NFFT,gausswin(501))).^2;
            [Qt0,~]=separation_warp(spectro_w(1:350,:));
            Q=[Q,Qt0];
        end
        %selection of the best t0 maximising the quality factor
        [~,temp]=max(Q); 
        ind0=L(temp);
        %warping of the signal using this t0
        stbb=str(ind0:ind2);
        N_ok=length(stbb); 
        [s_w, fs_w]=warp_temp_exa(stbb,fs,time(ind0),1);
        %computation of the spectrogram using a small value of sigma
        wind=gausswin(501);
        wind=wind/norm(wind);
        tfr_w=tfrstft(s_w,1:length(s_w),NFFT,wind); 
        spectro_w=abs(tfr_w).^2;
        %creation of the mask (see details below)
        [~,masquep]=separation_warp(spectro_w(1:350,:));
        %completion with zeros
        masque=[masquep;zeros(NFFT-350,length(s_w))];
        %masked spectrogram
        mode_rtf_warp=masque.*tfr_w;
        %inverse spectrogram
        [mode_temp_warp]=real(sum(mode_rtf_warp,1))/NFFT/max(wind)*2;
        %inverse warping
        stib=iwarp_temp_exa(mode_temp_warp,fs_w,time(ind0),1,fs,N_ok);
        %completion with zeros
        stiapp=[zeros(ind0-1,1);stib;zeros(length(s_t)-ind2,1)];
        %add to mode_t_app the new mode
        mode_t_app=[mode_t_app,stiapp]; 
        %remove it from the signal 
        str=str-stiapp.'; 
        %%%%%%%%%%%%
        %%% plot %%% (comment to avoid the plot)
        %%%%%%%%%%%%
        figure 
        sgtitle(sprintf('Extraction of mode %d', nbmode+1-nb))
        subplot 131
        %first plot: warped spectrogram and bassins
        spectro=smooth2a(spectro_w(1:350,:),4,2);
        I=spectro/max(max(spectro)); 
        I=1-I.*(I>=0.01); 
        WS=watershed(I); 
        N=double(max(max(WS)));
        Ymax=0*(1:N); 
        Xmax=0*(1:N); 
        for i=1:N
            S=spectro.*(WS==i); 
            [~,b]=max(S,[],'all');
            [i1,i2]= ind2sub(size(S),b); 
            Ymax(i)=i1; 
            Xmax(i)=i2; 
        end
        [Ymax,ind]=sort(Ymax,'descend');
        Xmax=Xmax(ind); 
        time_w=0:1/fs_w:(length(s_w)-1)/fs_w;
        freq_w=0:fs_w/NFFT:fs_w-fs_w/NFFT; 
        ylimM=min(freq_w(Ymax(1))+10,freq_w(350)); 
        xlimM=time_w(max(Xmax))+0.6; 
        imagesc(time_w,freq_w(1:350),spectro.*(WS>=1)+max(max(spectro)).*(WS==0))
        hold on 
        plot(time_w(Xmax),freq_w(Ymax),'r.','MarkerSize',14)
        axis xy 
        ylim([0,ylimM])
        xlim([0,xlimM])
        xlabel('Times (s)')
        ylabel('Frequency (Hz)')
        title('Warped spectrogram and bassins')
        %second plot: selected bassins
        subplot 132
        imagesc(time_w,freq_w(1:350),spectro.*masquep)
        axis xy 
        ylim([0,ylimM])
        xlim([0,xlimM])
        xlabel('Times (s)')
        ylabel('Frequency (Hz)')
        title('Masked warped spectrogram')
        %third plot: spectrogram of the modal component 
        subplot 133
        spectro=abs(tfrstft(stiapp.',1:length(time),NFFT,gausswin(51))).^2;
        freq=(0:NFFT-1)*fs/NFFT;
        imagesc(time, freq, spectro)
        axis xy
        ylim([0 fs/2])
        xlim_plots=[6.5 7.5];
        xlim(xlim_plots)
        xlabel('Times (s)')
        ylabel('Frequency (Hz)')
        title('Unwarped spectrogram and dispersion curve')
        vg=pek_vgb(freq,5-nb,5-nb,1500,1600,1000,1500,100);
        hold on 
        plot(10000./vg, freq,'k')
        pause(.1)
        %%%%%%%%%%%%
        %%% end plot %%% 
        %%%%%%%%%%%%
    end
    %the rest of the signal is associated to the last mode
    stiapp=str; 
    %%%%%%%%%%%%
    %%% plot %%% (comment to avoid the plot)
    %%%%%%%%%%%%
    figure
    sgtitle('Extraction of mode 1')
    spectro=abs(tfrstft(stiapp.',1:length(time),NFFT,gausswin(51))).^2;
    freq=(0:NFFT-1)*fs/NFFT;
    imagesc(time, freq, spectro)
    axis xy
    ylim([0 fs/2])
    xlim_plots=[6.5 7.5];
    xlim(xlim_plots)
    xlabel('Times (s)')
    ylabel('Frequency (Hz)')
    title('Remaining spectrogram and dispersion curve')
    vg=pek_vgb(freq,1,1,1500,1600,1000,1500,100);
    hold on 
    plot(10000./vg, freq,'k')
    pause(.1)
    %%%%%%%%%%%%
    %%% end plot %%% 
    %%%%%%%%%%%%
    mode_t_app=[mode_t_app,stiapp.']; 
    mode_t_app=mode_t_app(:,nbmode:-1:1);
end

function [Q,Mask]=separation_warp(spectro)
    %smooth the spectrogram to find the drainage bassins
    spectro=smooth2a(spectro,4,2);
    %convert the spectrogram into a black and white picture
    I=spectro/max(max(spectro)); 
    %remove near zero parts
    I=1-I.*(I>=0.01); 
    %do a watershed transform
    WS=watershed(I); 
    %number of drainage bassins
    N=double(max(max(WS)));
    %list of local max
    Lmax=0*(1:N); 
    %list of their positions
    Ymax=0*(1:N); 
    Xmax=0*(1:N); 
    %completion of the three lists
    for i=1:N
        S=spectro.*(WS==i); 
        [a,b]=max(S,[],'all');
        Lmax(i)=a; 
        [i1,i2]= ind2sub(size(S),b); 
        Ymax(i)=i1; 
        Xmax(i)=i2; 
    end
    %sort bassins by descengin frequencies 
    [Ymax,ind]=sort(Ymax,'descend');
    Xmax=Xmax(ind); 
    Lmax=Lmax(ind);
    %computation of the quality factor Q
    if length(Lmax)==1
        Q=0;
    elseif Lmax(1)<Lmax(2) 
        Q=0; 
    else
        i=2;
        while Lmax(i)<Lmax(i-1)
            if i==length(Lmax)
                break
            end
            i=i+1;    
        end
        Q=Lmax(i)-Lmax(i-1);
    end
    %computation of the mask
    i=1; 
    Mask=spectro*0; 
    %add the first bassin
    Mask=Mask+(WS==ind(1)); 
    %add the following bassin if their loc max is smaller
    while i<N && Lmax(i)>Lmax(i+1)
        i=i+1;
        Mask=Mask+(WS==ind(i)); 
    end
end


function [freq,tm]=meth_max(spectro,NFFT,time,fs)
    %set of computed frequencies 
    freq=(0:NFFT-1)*fs/NFFT;
    %list of associated max
    E=spectro(:,1)*0; 
    %list of the position of the local max
    tm=E*0; 
    for i=1:length(E)
        [e,a]=max(spectro(i,:)); 
        E(i)=e; 
        tm(i)=time(a); 
    end
    %we cut the symmetric par in freq and tm 
    freq=freq(1:NFFT/2-1); 
    tm=tm(1:NFFT/2-1); 
end

function x=J(r,c1,c2,rho1,rho2,D,dt,freq,T,AMP) 
    %theoretical dispersion curves
    temp=(r./pek_vgb(freq,1,4,c1,c2,rho1,rho2,D)); 
    x=0; 
    %remove all points with Nan values 
    temp(isnan(temp))=0; 
    for i=1:4
        x=x+norm(rmmissing(AMP(i,:).*(T(i,:)+dt-temp(i,:)))).^2; 
    end
end
