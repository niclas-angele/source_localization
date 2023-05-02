function [tfr,t,f] = tfrstft(x,t,N,h,trace);
%TFRSTFT Short time Fourier transform.
%	[TFR,T,F]=TFRSTFT(X,T,N,H,TRACE) computes the short-time Fourier 
%	transform of a discrete-time signal X. 
% 
%	X     : signal.
%	T     : time instant(s)          (default : 1:length(X)).
%	N     : number of frequency bins (default : length(X)).
%	H     : frequency smoothing window, H being normalized so as to
%	        be  of unit energy.      (default : Hamming(N/4)). 
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR   : time-frequency decomposition (complex values). The
%	        frequency axis is graduated from -0.5 to 0.5.
%	F     : vector of normalized frequencies.
%
%	Example :
%	 sig=[fmconst(128,0.2);fmconst(128,0.4)]; tfr=tfrstft(sig);
%	 subplot(211); imagesc(abs(tfr));
%	 subplot(212); imagesc(angle(tfr));
% 
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

%	F. Auger, May-August 1994, July 1995.
%	Copyright (c) 1996 by CNRS (France).
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


% This matlab function has been extracted from the Time Frequency Toolbox
% (TFTB) by Julien Bonnel (Woods Hole Oceanographic Institution) on August
% 31, 2019. The original toolbox can be downloaded here:
% http://tftb.nongnu.org/ (last seen on August 31, 2019).
% The matlab function has been slightly modified by Julien Bonnel so that
% it can be used without the full toolbox. Each modification is clearly
% identified using comments.


[xrow,xcol] = size(x);
if (nargin < 1),
 error('At least 1 parameter is required');
elseif (nargin <= 2),
 N=xrow;
end;

%%% Lines added by Bonnel
if ~isvector(x)
    error('X must be a vector');
end
if (xcol~=1)
    x=x.';
    [xrow,xcol] = size(x);
end
%%%% end of Bonnel's addition


hlength=floor(N/4);
hlength=hlength+1-rem(hlength,2);

if (nargin == 1),
%  t=1:xrow; h = tftb_window(hlength); trace=0; %%% line commented by Bonnel
 t=1:xrow; h = hamming(hlength); trace=0;  %%% line added by Bonnel
elseif (nargin == 2) | (nargin == 3),
%  h = tftb_window(hlength); trace = 0;  %%% line commented by Bonnel
 h = hamming(hlength); trace = 0;  %%% line added by Bonnel
elseif (nargin == 4),
 trace = 0;
end;

if (N<0),
 error('N must be greater than zero');
end;
[trow,tcol] = size(t);
if (xcol~=1),
 error('X must have one column');
elseif (trow~=1),
 error('T must only have one row'); 
elseif (2^nextpow2(N)~=N),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

[hrow,hcol]=size(h); Lh=(hrow-1)/2; 
if (hcol~=1)|(rem(hrow,2)==0),
 error('H must be a smoothing window with odd length');
end;

if N<length(h) %%% line added by Bonnel
    error('N must be greater than the window length') %%% line added by Bonnel
end %%% line added by Bonnel

h=h/norm(h);

tfr= zeros (N,tcol) ;  
if trace, disp('Short-time Fourier transform'); end;
for icol=1:tcol,
 ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
 indices= rem(N+tau,N)+1; 
 if trace, disprog(icol,tcol,10); end;
 tfr(indices,icol)=x(ti+tau,1).*conj(h(Lh+1+tau));
end;
tfr=fft(tfr); 
if trace, fprintf('\n'); end;

if (nargout==0),
%  tfrqview(abs(tfr).^2,x,t,'tfrstft',h);  %%% line commented by Bonnel
elseif (nargout==3),
 if rem(N,2)==0, 
  f=[0:N/2-1 -N/2:-1]'/N;
 else
  f=[0:(N-1)/2 -(N-1)/2:-1]'/N;  
 end;
end;
