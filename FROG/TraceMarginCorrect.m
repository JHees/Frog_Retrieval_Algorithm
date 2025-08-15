function S = TraceMarginCorrect(S, D, F, Fc, F_sp, spectrum)
    % TraceMarginCorrect. This function corrects the spectral margins of a SHG FROG trace.
    % 此功能用于校正SHG FROG trace的光谱边缘。
    % It adjusts the spectral margins to better align with the measured spectrum, ensuring more accurate trace reconstruction.
    % Usage:
    %     S = TraceMarginCorrect(S, D, F, Fc, F_sp, spectrum)
    %
    % Input:
    %     S: Input SHG FROG trace (size: N x N)
    %     D: Delay axis of the FROG trace
    %     F: Frequency axis of the FROG trace
    %     Fc: Central frequency of the trace
    %     F_sp: Frequency spectrum axis
    %     spectrum: Measured spectrum (vector)
    %
    % Output:
    %     S: Corrected FROG trace (size: N x N)
    %
    % Note:
    %     This function is specifically designed for SHG FROG traces. It uses Fourier transforms to compute marginal distributions
    %     and corrects the trace based on the ratio of the convolved measured spectrum to the computed marginal.
    %     For detailed information, see Chapter 10, Page 206 in:
    %     Trebino, R. "Frequency-Resolved Optical Gating: The Measurement of Ultrashort Laser Pulses."
    %     (Springer US, Boston, MA, 2000). DOI: 10.1007/978-1-4615-1181-6.
    %
    % AUTHOR: Huang Xingzhao, June 30, 2024

    spectrum = spectrum(:);
    % [~, maxPsp] = max(spectrum);
    % F_sp = F_sp(:)-F_sp(maxPsp);

    N = size(S, 1);
    mg = marginal(S);

    D_sp = FTconvert(F_sp);
    spectrum_rs = FreqTransfer(F_sp, spectrum, F + Fc);
    spectrum_conv = conv(spectrum_rs, spectrum_rs, "same");
    spectrum_conv = spectrum_conv / max(spectrum_conv);

    ratio = spectrum_conv ./ mg;
    ratio(mg <= 1e-4) = 0;
    S = S .* ratio;
end
