function [tbp,t_width,f_width] = TBP(T, E, method)
if nargin == 2
    method = 'fwhm';
end
[F, T] = FTconvert(T);
E_fq = fftshift(fft(E));
t_width = PulseWidth(T, E .* conj(E), method);
w_width = PulseWidth(2*pi*F, E_fq .* conj(E_fq), method);
f_width = w_width / 2 / pi;
tbp = t_width * w_width;
end
