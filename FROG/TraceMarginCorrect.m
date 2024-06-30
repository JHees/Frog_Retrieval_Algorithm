function S = TraceMarginCorrect(S, D, F, Fc, F_sp, spectrum)
% Only Shg frog can do
spectrum = spectrum(:);

N = size(S, 1);
iFT_ma = 1 / N * exp(-2i * pi * D .* F.');
FT_ma = exp(2i * pi * D.' .* (F - Fc));
marginal = sum(S, 2);
marginal = real(FT_ma * iFT_ma * marginal);
marginal = marginal / max(marginal);

D_sp = FTconvert(F_sp);
iFT_sp = 1 / N * exp(-2i * pi * D_sp .* F_sp.');
FT_sp = exp(2i * pi * D_sp.' .* F);
spectrum_rs = real(FT_sp * iFT_sp * spectrum);
spectrum_conv = conv(spectrum_rs, spectrum_rs, "same");
spectrum_conv = spectrum_conv / max(spectrum_conv);

% spectrum = spectrum / max(spectrum);
% spectrum_rs = spectrum_rs / max(spectrum_rs);
% plot(F, marginal, F, spectrum_conv, F_sp, spectrum, F, spectrum_rs)
% legend(["mar", "sp\_conv", "sp\_rc", "sp\_ori"])

ratio = spectrum_conv ./ marginal;
S = S .* ratio;
end
