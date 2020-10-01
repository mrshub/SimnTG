% sim_outAdd.m
% Ronald Ouwerkerk 2020 created by modifying sim_dAdd.m from Jamie Near.
%
% USAGE:
% out = sim_dAdd(out1,out2, factor)
% 
% DESCRIPTION:
% Add together two simulation output structures.  
% This function can be used to create mixed spectra. 
% 
% INPUTS:
% out1        = first input density matrix to be added.
% out2        = second input density matrix to be added.
% factor    = 1 for sum (default), -1 for diff
%
% OUTPUTS:
% out = sum of out1 and out2.
%%

%  Ronald Ouwerkerk NIDDK/NIH 2020


function out = sim_outAdd(out1,out2,factor)

if nargin<3
    factor=1;
end

%If out1 is an empty structure (or not a structure at all), make out = out2; otherwise, add them:
%Check whether out1 is 
if ~isstruct(out1) && isempty(out1)  %First input is not a structure and is empty
   out=out2;
elseif ~isstruct(out1) && ~isempty(out1)  %First input is not a structure but is not empty (can't use this). 
    error('ERROR:  out1 is not empty.  ABORTING!!');
elseif isstruct(out1) && isempty(out1(1))  %First input is an empty structure.
    out=out2;
elseif isstruct(out1) && ~isempty(out1(1)) %First input is a non-empty structure (actually doing some addition). 
    if out1(1).n ~= out2(1).n  %Make sure the fids and specs to add are the same size
        error('ERROR:  can only add output structutes with fids of the same size!!  ABORTING!!');
    end
    [~,Mfids1] = size( out1(1).fids);
    [~,Mfids2] = size( out2(1).fids);
    [~,Mspecs1] = size( out1(1).specs);
    [~,Mspecs2] = size( out2(1).specs);
    
    %Make sure the number of fids and specs is the same in both structures or struct arrays
    if (Mfids1 ~= Mfids2) || (Mspecs1 ~=Mspecs2) 
        error('ERROR:  can only add outpu t structutes with the same number of fids and spectra!!  ABORTING!!');
    end
    % Use out1 as the basis structure for output
    out = out1;
    % now sum the fid and spec elements of the structures
    for ii=1:length(out1)
        out(ii).fids =out1(ii).fids +out2(ii).fids  .*factor;
        out(ii).specs=out1(ii).fids +out2(ii).specs .*factor;
    end
end



