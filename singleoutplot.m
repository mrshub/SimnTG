function [figh, ploth]  =singleoutplot( out , lb , phasecor)
% [figh, wfh]  =singleoutplot( out , lb , T2, phasecor)
% Create a plot of a simulation output spectrum
% out = 
% 
%   1×Nspectra struct array with fields:
% 
%     fids
%     t
%     ppm
%     specs
%     spectralwidth
%     dwelltime
%     n
%     linewidth
%     Bo
%     txfrq
%     sz
%     date
%     dims
%     averages
%     rawAverages
%     subspecs
%     rawSubspecs
%     flags
%     seq
%     te
%     sim
%%

%  Ronald Ouwerkerk NIDDK/NIH 2020


if nargin < 2    
    lb = 2;
end

if nargin < 3
    % no phase correction
    phi0 = 0.0;    
else
    phi0 = pi*phasecor/180;
end

sw = out(1).spectralwidth;
n = out(1).n;
field = out(1).Bo;
F0 = 42.57*field;
Nspecpts = 2*n;

fids = out(1).fids*exp( 1i*phi0);
gfids = gausmult( fids , lb, sw );
F = fftshift( fft( gfids , Nspecpts), 1);

fhz = linspace( -sw/2, sw/2,  Nspecpts);
fppm = fhz/F0 + 4.7;
figh = figure;
hold on;
ploth = plot( fppm, F , 'k');
setgcf;
setgca;
xlabel( 'ppm');
set( gca, 'Xdir', 'reverse');
set( gca, 'Box', 'off');
set( gca, 'Xlim', [0, 7]);