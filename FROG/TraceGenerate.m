function [S, sig_sp, sig, shifting_outer, outer] = TraceGenerate(P, G)
    % TraceGenerate. This function generates the FROG trace of two input pulses.
    % 获得两个pulse的FrogTrace
    % This function computes the complex spectra and time-domain signal of two interacting pulses,
    % which is useful for applications like spectral shearing interferometry and pulse characterization.
    % Usage:
    %     [S, sig_sp, sig, shifting_outer, outer] = TraceGenerate(P, G)
    %
    % Input:
    %     P, G: Input pulses (Pulse P and Pulse G)
    %
    % Output:
    %     S: FROG trace (Intensity profile of the signal)
    %     sig_sp: Complex spectra of the signal
    %     sig: Time-domain signal
    %     shifting_outer: Shifted outer product in time domain
    %     outer: Outer product of P and G
    %
    % Note:
    %     Detailed calculations require thorough derivation of the trace.
    %
    % AUTHOR: Huang Xingzhao, June 30, 2024

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
