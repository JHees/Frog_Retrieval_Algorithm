function [S, sig_sp, sig, shifting_outer, outer] = TraceGenerate(P, G)
%获得两个pulse的FrogTrace
%   Input:
%       P, G: input pulse
%   Output:
%       S:   Trace
%       sig_sp: complex spectra
%
% AUTHOR:  Huang Xingzhao

if nargin == 1
    G = P;
end
P = P(:);
G = G(:);
N = numel(P);

outer = (P) .* (G.');
% shg signal in time domain
shifting_outer = outer((1:N)' + mod((0:N:N^2 - 1) + (0:N:N^2 - 1)', N^2)); % maybe this line could be faster
sig = circshift(shifting_outer, N / 2, 2);

sig_sp = circshift(fft(sig, [], 1), N / 2, 1);

S = sig_sp .* conj(sig_sp);
S = S ./ max(S, [], "all");
end
