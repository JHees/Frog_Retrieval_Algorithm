clear
% units: fs,nm,PHz
c = 299.77;

[F, T] = FTconvert((0:1:5000 - 1)');
N = floor(numel(T) / 2) * 2;
FF = FTconvert(F);
wavelength = 850;
f0 = 0;

Delay = 200;
gauss = @(t0, tou)exp(- (T - t0).^2 ./ tou.^2 * 2 * log(2)) .* exp(1i * 2 * pi * (f0 .* T));

A = [1, 1];
t0 = [0, Delay];
FWHM = [100, 100];
P = sum(gauss(t0, FWHM) .* A, 2);
% P = flip(P);
fm = 0.03;
dm = 600;

% clear
% load('pulse_set.mat')
% f0 = 0;
% [F, T] = FTconvert(time' / 1e-15);
% dT = T(2)-T(1);
% dT = 1;
% P = pulse_set(3, :).';
% P = P ./ max(abs(P));
% P = flip(P);
% fm = 0.3;
% dm = 75;

[I_withoutNoise] = TraceGenerate(P);
pulse_spectrum = abs(fftshift(fft(P))).^2;
[~, maxPsp] = max(pulse_spectrum);
Fc_sp = F(maxPsp);
% [P, P_sp, df0] = removeFirstOrderPhase(T, P); % P_sp has been removed f0
%% add Noise
bitDepth = 16;
maxLevel = 2^bitDepth -1;

alpha = 0.1;
additiveNoiseLevel = alpha * maxLevel;

if alpha == 0
    I = I_withoutNoise;
else
    I = I_withoutNoise .* (1 + alpha * randn(N)) + alpha * poissrnd(additiveNoiseLevel, N) / additiveNoiseLevel;
    I(I < 0) = 0;
    I = round(I .* maxLevel) ./ maxLevel;
    I(I>1)=1;
end
[I_withoutNoise, Do_wn, Fo_wn, Fc_wn] = TraceResample(I_withoutNoise, T, F - 2 * Fc_sp, 128, 'eps', 1e-4);
figure(5)
plotTrace(I_withoutNoise, Do_wn, Fo_wn);
%% Delay Correct
dt = TraceDelayCorrect(I, 0.5, F, 50);
%% denoise
S_denoise = TraceDenoise(I, dt, F, 5, 600, 1);
%% resample
[S_resample, Do, Fo, Fc] = TraceResample(S_denoise, dt, F - 2 * Fc_sp, 128, 'eps', 1e-4);
S = TraceMarginCorrect(S_resample, Do, Fo, Fc, F - Fc_sp, pulse_spectrum);

% figure(6)
% plotTrace(I, T, F);
figure(7)
plotTrace(S_resample, Do, Fo);
%% Solve
[pcgp_svd] = solveFrog_PCGPA(S, [], [2, 200], 1e-4, 'svd');
[solver, ErrorBar] = TraceSolver(S, [], [50, 200], 5, 5, 1e-4);
figure(1)
plot(T, abs(P) / max(abs(P)), Do, [abs(pcgp_svd.P), abs(solver.P)])
legend(["origin", "pcgp svd", "solve"])
figure(2)
semilogy( ...
    1:sum(pcgp_svd.iter), pcgp_svd.err, ...
    1:sum(solver.iter), solver.err)
legend(["pcgp svd", "solve"])

figure(3)
clf
plotPulse(T, P, [], dm, fm);
plotPulse(Do, solver.P, ErrorBar, dm, fm);
plotPulse(Do, pcgp_svd.P, [], dm, fm);
figure(8)
retrievedTrace = TraceGenerate(solver.P, solver.P);
plotTrace(retrievedTrace, Do, Fo);
realError = TraceError(I_withoutNoise, retrievedTrace)
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
