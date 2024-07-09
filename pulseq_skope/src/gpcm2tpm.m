function op = gpcm2tpm (ip)
% GPCM2TPM Convert gradients in Gauss/cm to those in Tesla/m.
%
% Usage: op = gpcm2tpm (ip)
%
% Returns
% -------
% op: output grad in T/m
%
% Expects
% -------
% ip: input grad in Gauss/cm
% 
% 
% See also: tpm2gpcm
% 

  op = 0.01*ip;
  
