function gsdata = gausmult(data, width, SW, tbegin )
%============================================================================
% function gsdata = gausmult(data, width, SW, tbegin);
% INPUT:
% data = complex time domain data in columns
% [m,n] = size(data) then m is the number of time points and n is the
% number of signals
% width : line broadening in Hz (default 2 Hz)
% SW is spectral width Hz ( default 1500)
% tbegin, tste begin time and step time on s ( default 0 and 1e-3 )
%==========================================================================
[m,n] = size(data);

if nargin < 2
width = 2;
end

if nargin < 3
SW = 1500;
end

if nargin < 4
tbegin = 0;
end

tstep = 1/SW;

% Fourier transform of a Gaussian
% FTexp( -ax^2) = sqrt( pi/a)* exp( -pi^2k^2/a)
% half height in freq domain is when exp( -pi^2k^2/a) = 0.5
% so ln( 2 ) = pi^2k^2/a, or aln(2)/pi^2 = k^2
% Substittute the desired line width lw for k
% ( lw/2 )^2 = -a*ln(2)/pi^2
% Now a = ( lw/2 pi )^2 /ln(2)
widthsq = (piwidth/2)^2/log(2);
time = (tbegin+(0:m-1)’*tstep);

if n ==1
timesq = time.^2;
gsdata = data.exp(-timesqwidthsq);
else
time = time( :, ones(1,n) );
timesq = time.^2;
gsdata = data.exp(-timesqwidthsq);
end

end
