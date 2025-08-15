clear
% units: fs,nm,PHz
c = 299.792458;

[F, T] = FTconvert((0:1:2000 - 0.1)');
N = floor(numel(T) / 2) * 2;
dT = T(2) - T(1);
dt = dT;
wavelength = 1e4;
f0 = c / wavelength;

Delay = 0.5/f0*0;
gauss = @(t0, tou)exp(- (T - t0).^2 ./ tou.^2 * 2 * log(2)).*exp(-2i*pi*f0*(-t0));

A = [ones(size(T)),ones(size(T)).*exp(-2i*pi*0.5)];
t0 = [0, Delay];
FWHM = [100, 100];
ps = gauss(t0, FWHM) .* A;
P = sum(ps, 2);
P = (P); 
fm = 0.03;
dm = 600;

figure(6)
clf
plotPulse(T, P, [], dm, fm);
% clf
% plot(T,real(ps(:,1)),T,real(ps(:,2)))

% clear
% load('pulse_set.mat')
% f0 = 0;
% [F, T] = FTconvert(time' / 1e-15);
% dT = T(2)-T(1);
% dt = dT;
% P = pulse_set(3, :).';
% G = pulse_set(4,:).';
% P = P ./ max(abs(P));
% P = flip(P);
% G = G ./ max(abs(G));
% G = flip(G);
% fm = 0.3;
% dm = 100;
%%
[I_withoutNoise_origin] = TraceGenerate(P);
pulse_spectrum = abs(fftshift(fft(P))).^2;
isConvert2Wlen = 0;

if isConvert2Wlen
    ww=50;
    dw = 0.5;
    wlen_I = wavelength/2-ww/2:dw:wavelength/2+ww/2;
    wlen_sp= wavelength-ww:dw:wavelength+ww;
    I_withoutNoise = FreqTransfer(F+2*f0,I_withoutNoise_origin,c./wlen_I);
    pulse_spectrum = FreqTransfer(F+f0,pulse_spectrum,c./wlen_sp);
else
    I_withoutNoise = I_withoutNoise_origin;
end
% figure(1)
% plot(wlen_I,I_withoutNoise(:,end/2+1),c./(F+f0*2),I_withoutNoise_origin(:,end/2+1))
% plotTrace(I_withoutNoise_origin, T, F);

% [P, P_sp, df0] = removeFirstOrderPhase(T, P); % P_sp has been removed f0
%% add Noise
bitDepth = 16;
alpha = 0;

if alpha == 0
    I = I_withoutNoise;
else
    maxLevel = 2^bitDepth -1;
    additiveNoiseLevel = alpha * maxLevel;
    I = I_withoutNoise .* (1 + alpha * randn(size(I_withoutNoise))) + alpha * poissrnd(additiveNoiseLevel, size(I_withoutNoise)) / additiveNoiseLevel;
    I(I < 0) = 0;
    I = round(I .* maxLevel) ./ maxLevel;
    I(I > 1) = 1;
end

if isConvert2Wlen
    [I_withoutNoise_re, Do_wn, Fo_wn, Fc_wn] = TraceResample(I_withoutNoise, dT, c ./ wlen_I, 128, 'eps', 1e-4);
    [I_withoutNoise_orre, Do_wn, Fo_wn, Fc_wn] = TraceResample(I_withoutNoise_origin, dT, c ./ wlen_I, 128, 'eps', 1e-4);
else
    [I_withoutNoise_re, Do_wn, Fo_wn, Fc_wn] = TraceResample(I_withoutNoise, dT, F, 128, 'eps', 1e-4);
end
figure(5)
plotTrace(I_withoutNoise_re, Do_wn, Fo_wn + Fc_wn);
return
%% Trace Correct
if isConvert2Wlen
    I_trans = FreqTransfer(c./wlen_I, I,F+f0*2);
    pulse_spectrum_trans = FreqTransfer(c./wlen_sp, pulse_spectrum,F+f0);
else
    I_trans = I;
    pulse_spectrum_trans = pulse_spectrum;
end
% figure(6)
% plot(c./(F+f0*2),I_trans(:,end/2+1),wlen_I,I(:,end/2+1),c./(F+f0*2),I_withoutNoise_origin(:,end/2+1))
%%
% dt = TraceDelayCorrect(I_trans, 1, F, 50)
if alpha == 0
    S_denoise = I_trans;
else
    S_denoise = TraceDenoise(I_trans, dt, F, 5, 600, 0.5);
end

% resample
[S_resample, Do, Fo, Fc] = TraceResample(S_denoise, dt, F, 128, 'eps', 1e-4);
% Dm = dm;
% N = 128;
% dd = Dm / N*2;
% Do = (-N/2:N/2-1)*dd;
% Fo = FTconvert(Do);
% [Dq, Fq] = meshgrid(Do, Fo);
% S_resample = interp2(T, F, S_denoise, Dq, Fq, 'linear', 0);
% S = S_resample;
S = TraceMarginCorrect(S_resample, Do, Fo, Fc, F, pulse_spectrum_trans);
%%
figure(7)
plotTrace(S, Do, Fo);
title("denoise")
%% Solve
[pcgp_svd] = solveFrog_PCGPA(S, [], [2, 200], 1e-4, 'svd');
[solver, ErrorBar] = TraceSolver(S, [], [5, 200], 3, 5, 1e-4);
%%
% figure(1)
% plot(T, abs(P) / max(abs(P)), Do, [abs(pcgp_svd.P), abs(solver.P)])
% legend(["origin", "pcgp svd", "solve"])

figure(2)
semilogy( ...
    1:sum(pcgp_svd.iter), pcgp_svd.err, ...
    1:sum(solver.iter), solver.err)
title("Retrieved Error")
xlabel("iter")
legend(["pcgp svd", "solve"])

figure(3)
clf
plotPulse(Do, solver.P, ErrorBar, dm, fm);
% plotPulse(Do, pcgp_svd.G, [], dm, fm);
plotPulse(T, P, [], dm, fm);
% subplot(2,1,1)
% legend(["pcgp svd","origin"])

figure(4)
clf
plotPulse(Do, pcgp_svd.P, [], dm, fm);
plotPulse(T, P, [], dm, fm);
subplot(2,1,1)
legend(["pcgp svd","origin"])

figure(8)
retrievedTrace = TraceGenerate(solver.P);
plotTrace(retrievedTrace, Do, Fo);
title("retrieved")
realError = TraceError(I_withoutNoise_re, retrievedTrace)
%% marginals
% figure(4)
% subplot(2,1,1)
% yyaxis left
% plot(Fo,sum(S_cutoff,2))
% ylim([0,Inf])
% yyaxis right
% marginal = conv(abs(P_sp).^2,abs(P_sp).^2,"same");
% plot(F,marginal)
% xlim([-fm,fm])
% ylim([0,Inf])
% % autocorrelation
% subplot(2,1,2)
% yyaxis left
% plot(Do,sum(S_cutoff,1))
% ylim([0,Inf])
% yyaxis right
% It = abs(P).^2;
% Iw = fftshift(fft(It));
% At = ifft(ifftshift(abs(Iw).^2));
% % At = conv(It,flip(It),'same');
% plot(T,fftshift(At))
% xlim([-dm,dm])
% ylim([0,Inf])
%% Solve and plot
% [van] = solveFrog_vanilla(S, [], 2000, 1e-5, "over_correction", "k");
% [gp] = solveFrog_GP(S, [], 1e-2, [50, 2000], 1e-5);
% [scgp] = solveFrog_shortcutGP(S, [], 2e-3, 3, [50, 2000], 1e-5);
% [pcgp] = solveFrog_PCGPA(S, [], [50, 2000], 1e-5); % 不收敛的话看看GP是否收敛，GP收敛的话延长GP初始化
% figure(1)
% plot(T, abs(P) / max(abs(P)), Do, [abs(van.P), abs(gp.P), abs(scgp.P), abs(pcgp.P), abs(pcgp_svd.P),abs(solver.P)])
% legend(["origin", "vanilla", "gp", "scgp", "pcgp", "pcgp svd","solve"])
% figure(2)
% semilogy(...
%     1:sum(van.iter), van.err, ...
%     1:sum(gp.iter), gp.err, ...
%     1:sum(scgp.iter), scgp.err, ...
%     1:sum(pcgp.iter), pcgp.err, ...
%     1:sum(pcgp_svd.iter), pcgp_svd.err, ...
%     1:sum(solver.iter) + sum(solver2.iter), [solver.err, solver2.err])
% legend(["vanilla", "gp", "scgp", "pcgp", "pcgp svd", "solve"])

% figure(3)
% clf
% plotPulse(T, P, [], dm, fm);
% plotPulse(Do, pcgp_svd.P, [], dm, fm);
% plotPulse(Do, pcgp_svd.G, dm, fm);
