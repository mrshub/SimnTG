function [figh, wfh]  = stackedoutplots( out , lb , T2, phasecor)
% [figh, wfh]  = stackedoutplots( out , lb , T2, phasecor )
% Create stacked plots of spectra
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

% Ronald Ouwerkerk NIDDK/NIH 2020
    
    
    
Nspectra = length( out );

% Set typical reference ppm for 1H
refppm = 4.7;
GAMMA = 42.58;

SW = out(1).spectralwidth; % in Hz
field = out(1).Bo; % in T
% F0 in MHz
F0 = field*GAMMA;
Nfidpts = out(1).n;

if nargin < 2    
    lb = 2;
end

if nargin < 3
    % no T2 decay as a function of TE
    T2 = 0;    
end

if nargin < 4
    % no phase correction
    phi0 = 0.0;    
else
    phi0 = pi*phasecor/180;
end


for ii = 1:Nspectra
    thisTE = out(ii).te;
    TEs(ii) = thisTE;
    % HACK ! put in a scaling factor to give all signals the same T2 decay
    if T2 >0
        T2scaler = exp(-thisTE/T2);
    else
        T2scaler = 1;
    end
    
    TRs(ii) = 10000;
    thisdata = T2scaler*out(ii).fids(:);
    rdata( :, ii ) = thisdata*exp( 1i*phi0);
    TEs(ii) = out(ii).te;
    TRs(ii) = 10000;
end

[Nfidpts, Nspecs] = size( rdata );
Nspecpts = max( 2048, 2*Nfidpts);

fhz = linspace ( -SW/2, SW/2, Nspecpts);
fppm = fhz/F0 + 4.7;

X = fppm;

Y = 1:Nspecs;
if length( unique( TEs ) ) > 1
    Y = TEs;
    ylabelstr = ('TE [ms]' );
elseif length( unique( TRs ) ) > 1
    Y = TRs; 
     ylabelstr = ('TR [ms]' );
end
       
% Gaussian line broadening 
if nargin < 2
    lb = 21.0;
end

gfids = gausmult( rdata, lb, SW );
F = fftshift( fft( gfids , Nspecpts), 1);

Z  = real(F)';
figh = figure;
wfh = waterfall( X, Y, Z );

set( figh, 'Color', 'w');
set( gca, 'Color', 'w');
ylabel( ylabelstr )
xlabel( 'ppm')
set( wfh, 'EdgeColor', 'k')

zlim = get( gca, 'Zlim');

maxz = max(Z(:));
minz = min(Z(:));
zrange = maxz- minz;
set( gca, 'Zlim', [ (minz-0.01*zrange), (maxz+0.02*zrange) ]);
set( gca, 'Xdir', 'reverse');




