function S = TraceDenoise(I, D, F, background_subtraction_line, corner_suppression_width, ft_low_pass_filter)
    % TraceDenoise. This function reduces noise in a FROG trace using multiple methods including background subtraction,
    % corner suppression, and Fourier transform low-pass filtering.
    % 噪声减少是通过多种方法来完成的，包括背景减法、角落压制和傅里叶变换低通滤波。
    % The function is essential for enhancing the quality of the trace, ensuring cleaner and more accurate analysis results.
    % Usage:
    %     S = TraceDenoise(I, D, F, background_subtraction_line, corner_suppression_width, ft_low_pass_filter)
    %
    % Input:
    %     I: Original FROG Trace (size: Fn x Dn)
    %     D: Original FROG Trace Delay Axis
    %     F: Original FROG Trace Frequency Axis
    %     background_subtraction_line: Number of lines for background noise calculation
    %     corner_suppression_width: Width for corner suppression in trace
    %     ft_low_pass_filter: Cutoff frequency for low-pass filtering in Fourier transform
    %
    % Output:
    %     S: Denoised FROG Trace (size: Fn x Dn)E
    %
    % AUTHOR: Huang Xingzhao, June 30, 2024
    [~, maxInd] = max(I, [], 'all');
    [FcInd, DcInd] = ind2sub(size(I), maxInd);
    F = F(:) - F(FcInd);
    D = D(:)' - D(DcInd);
    S = I;

    if ~isempty(background_subtraction_line)
        S = background_subtraction(S, background_subtraction_line);
    end

    if ~isempty(corner_suppression_width)
        S = corner_suppression(S, D, F, corner_suppression_width * 2);
    end
    
    if ~isempty(ft_low_pass_filter)
        S = fourier_low_pass_filter(S, ft_low_pass_filter);
    end
end

function S = background_subtraction(I, n)
    S = I - mean(I(:, [1:n, end - n:end]), 2);
    S(S < 0) = 0;
end

function S = corner_suppression(I, D, F, D_width)
    marginal = sum(I, 2);
    marginal = marginal - median(marginal);
    marginal = marginal / max(marginal);

    autocorrelation = sum(I, 1);
    autocorrelation = autocorrelation - median(autocorrelation);
    autocorrelation = autocorrelation / max(autocorrelation);

    marginal_width = PulseMainWidth(F, marginal, 1e-4);
    ac_width = PulseMainWidth(D, autocorrelation, 1e-4);

    % N_width = sum(abs(D) < D_width);
    % F_width = N_width / D_width / 4;

    F_width = marginal_width * D_width / ac_width;

    super_Gauss = exp(-16 * log(2) * ((F / F_width).^2 + (D / D_width).^2).^2);
    S = I .* super_Gauss;
end

function S = fourier_low_pass_filter(I, rou)
    S_ft = fft2(I);
    S_ft = fftshift(S_ft);
    [NF, ND] = size(I);
    [F, D] = meshgrid((-NF / 2:NF / 2 - 1), (-ND / 2:ND / 2 - 1));
    dis = sqrt((2 * F / NF).^2 + (2 * D / ND).^2);
    filter = dis <= rou;
    S_ft_filter = S_ft .* filter;
    S = real(ifft2(ifftshift(S_ft_filter)));
    S(S < 0) = 0;
end
