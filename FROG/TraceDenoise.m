function S = TraceDenoise(I, D, F, background_subtraction_line, corner_suppression_width, ft_low_pass_filter)
[~, maxInd] = max(I, [], 'all');
[FcInd, DcInd] = ind2sub(size(I), maxInd);
F = F(:) - F(FcInd);
D = D(:)' - D(DcInd);

S = background_subtraction(I, background_subtraction_line);
S = corner_suppression(S, D, F, corner_suppression_width * 2);
S = fourier_low_pass_filter(S, ft_low_pass_filter);
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
