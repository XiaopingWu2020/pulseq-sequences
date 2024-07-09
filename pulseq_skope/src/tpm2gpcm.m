function op = tpm2gpcm (ip)
% TPM2GPCM Convert gradients in Tesla/m to those in Gauss/cm.
%
% Usage: op = tpm2gpcm (ip)
%
% Returns
% -------
% op: output grad in G/cm
%
% Expects
% -------
% ip: input grad in T/m
%
%
% See also: gpcm2tpm
%
%
% Copyright (C) 2008 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Thu Apr 10 13:00:28 2008
%

op = 100*ip;
