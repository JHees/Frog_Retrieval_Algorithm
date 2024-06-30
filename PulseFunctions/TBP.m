function [tbp, t_width, f_width] = TBP(T, E, method)
    % TBP. This function calculates the time-bandwidth product (TBP) of a signal.
    % 此函数计算信号的时间带宽积 (TBP)。
    % TBP is a measure of the signal's spread in both time and frequency domains, using specified width calculation methods.
    % Usage:
    %     [tbp, t_width, f_width] = TBP(T, E, method)
    %
    % Input:
    %     T: Time domain data
    %     E: Electric field amplitude (time domain)
    %     method: Method for width calculation ('fwhm', 'hw1e', 'rms', 'equ')
    %
    % Output:
    %     tbp: Time-bandwidth product
    %     t_width: Width of the pulse in time domain
    %     f_width: Width of the pulse in frequency domain
    %
    % Note:
    %     The function converts time domain data to frequency domain, calculates time and frequency widths using the specified method,
    %     and then computes the TBP as the product of these widths.
    %
    % AUTHOR: Huang Xingzhao, June 30, 2024

    if nargin == 2
        method = 'fwhm';
    end
    [F, T] = FTconvert(T);
    E_fq = fftshift(fft(E));
    t_width = PulseWidth(T, E .* conj(E), method);
    w_width = PulseWidth(2 * pi * F, E_fq .* conj(E_fq), method);
    f_width = w_width / 2 / pi;
    tbp = t_width * w_width;
end
