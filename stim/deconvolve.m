function [constrainedFoopsi] = deconvolve(dF)
%DECONVOLVE Deconvolved dF/F traces using constrained FOOPSI
%   Constrained FOOPSI https://github.com/epnev/constrained-foopsi
%	Based on wrapper function calcConstrainedFoopsi from Matthias
%   Variables:
%   y:      raw fluorescence data (vector of length(T))
%   c:      denoised calcium concentration (Tx1 vector)
%   b:      baseline concentration (scalar)
%  c1:      initial concentration (scalar)
%   g:      discrete time constant(s) (scalar or 2x1 vector)
%  sn:      noise standard deviation (scalar)
%  sp:      spike vector (Tx1 vector)

%   USAGE:
%   [c,b,c1,g,sn,sp] = constrained_foopsi(y,b,c1,g,sn,OPTIONS)
%   The parameters b,cin,g,sn can be given or else are estimated from the data

%   OPTIONS: (stuct for specifying options)
%         p: order for AR model, used when g is not given (default 2)
%    method: methods for performing spike inference
%   available methods: 'dual' uses dual ascent
%                       'cvx' uses the cvx package available from cvxr.com (default)
%                      'lars' uses the least regression algorithm 
%                     'spgl1' uses the spgl1 package available from
%                     math.ucdavis.edu/~mpf/spgl1/  (usually fastest)
%   bas_nonneg:   flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y)
%   noise_range:  frequency range over which the noise power is estimated. Default [Fs/4,Fs/2]
%   noise_method: method to average the PSD in order to obtain a robust noise level estimate
%   lags:         number of extra autocovariance lags to be considered when estimating the time constants
%   resparse:     number of times that the solution is resparsened (default 0). Currently available only with methods 'cvx', 'spgl'
%   fudge_factor: scaling constant to reduce bias in the time constant estimation (default 1 - no scaling)

opt.p = 4; % Order of autoregressive model.
opt.method = 'cvx'; % spgl1 takes very long
opt.bas_nonneg = 1; % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y)
opt.noise_range = [1/4,1/2]; % frequency range over which the noise power is estimated. Default [Fs/4,Fs/2]
opt.noise_method = 'logmexp';
opt.lags = 20;
opt.ressparse = 0;
opt.fudge_factor = 0.98;

unitaryDF = 0.1; % DF of one spike, for Selmaan's scaling function.

nTraces = size(dF, 1);

c = zeros(size(dF));
b = zeros(nTraces, 1);
c1 = zeros(nTraces, 1);
g = zeros(nTraces, opt.p);
sn = zeros(nTraces, 1);
deconvInSpikes = zeros(size(dF));

parfor i = 1:size(dF, 1)
    try
        [c(i, :), b(i), c1(i), g(i, :), sn(i), sp] = constrained_foopsi(dF(i, :), [], [], [], [], opt);
        maxPulse = max(impulseAR(g(i, :)));
        deconvInSpikes(i, :) = sp*maxPulse/unitaryDF;
    catch err
        warning('Error in parfor: \n%s\n%s\nline %d', err.identifier, err.stack(1).file, err.stack(1).line);
    end
end

constrainedFoopsi.c = c;
constrainedFoopsi.deconvInSpikes = deconvInSpikes;
constrainedFoopsi.b = b;
constrainedFoopsi.c1 = c1;
constrainedFoopsi.g = g;
constrainedFoopsi.sn = sn;

function impulseResponse = impulseAR(p)

impulseResponse = zeros(1e3,1);
impulseResponse(50) = 1;
p = p(:)';

for ind = 51:1e3
    impulseResponse(ind) = p*impulseResponse(ind-1:-1:ind-length(p));
end