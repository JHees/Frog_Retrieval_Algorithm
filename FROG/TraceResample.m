function [Io, Do, Fo, Fc] = TraceResample(I, D, F, N, isPower2Required, varargin)
    % TraceResample. This function generates an NxN trace by resampling the original trace.
    % 由于坐标需要满足傅里叶变换(FTconvert)，但直接进行插值会导致trace失真. 此方法将重新抽样D坐标，并通过傅里叶变换转换F坐标，得到NxN trace
    % Coordinates must satisfy the Fourier transform (FT convert), but direct interpolation could distort the trace.
    % This method resamples the D (Delay) coordinates and uses Fourier transform to convert the F (Frequency) coordinates,
    % resulting in an NxN trace.
    % Usage:
    %     [I, Do, Fo] = cutoffTrace(I, D, F)
    %     [I, Do, Fo] = cutoffTrace(I, D, F, N)
    %     [I, Do, Fo] = cutoffTrace(I, D, F, [maxN, minN])
    %     [I, Do, Fo] = cutoffTrace(I, D, F, N, true)
    %     [I, Do, Fo] = cutoffTrace(I, D, F, 'Dm', Dm)
    %     [I, Do, Fo] = cutoffTrace(I, D, F, [maxN, minN], 'Dm', Dm)
    %     [I, Do, Fo] = cutoffTrace(I, D, F, 'eps', eps)
    %     [I, Do, Fo] = cutoffTrace(I, D, F, [maxN, minN], 'eps', eps)
    %
    % Input:
    %     I: Original Frog Trace (size: Fn x Dn)
    %     D: Original Frog Trace Delay Axis
    %     F: Original Frog Trace Frequency Axis
    %     N: Output trace size (default: 128)
    %     isPower2Required: Indicates whether the output trace size (N) must be a power of 2.
    %     Dm: The maximum delay time to be reserved
    %     eps: The minimum value that will be included in the trace.
    %
    % Output:
    %     Io: Output Frog Trace (size: Fo x Do)
    %     Do: Output Delay Axis
    %     Fo: Output Frequency Axis
    %     Fc: Frequency center, the maximum value on the frequency axis of the Trace
    %
    % Update June 21, 2024:
    %     The code has been updated to ensure that the output trace maintains the correct sample rate
    %     for accurate frequency and delay axis representation after Fourier transform.
    %     For detailed information, see Chapter 10, Page 212 in:
    %     Trebino, R. "Frequency-Resolved Optical Gating: The Measurement of Ultrashort Laser Pulses."
    %     (Springer US, Boston, MA, 2000). DOI: 10.1007/978-1-4615-1181-6.
    %
    % AUTHOR:  Huang Xingzhao, June 15, 2024

    p = inputParser;
    p.addOptional("N", 128);
    p.addOptional("isPower2Required", []);
    p.addParameter("eps", []);
    p.addParameter("Dm", []);
    if nargin <= 4
        isPower2Required = [];
    end
    p.parse(N, isPower2Required, varargin{:});

    Dm = p.Results.Dm;
    eps = p.Results.eps;
    isPower2Required = ~isempty(p.Results.isPower2Required);
    isInputN = isscalar(p.Results.N);
    if isInputN
        maxN = p.Results.N;
        minN = p.Results.N;
    else
        maxN = p.Results.N(1);
        minN = p.Results.N(2);
    end

    [NF, ND] = size(I);
    if isscalar(D)
        D = (-ND / 2:ND / 2 - 1)' .* D;
    end

    [~, maxTraceIndex] = max(I, [], 'all');
    [FcIndex, DcIndex] = ind2sub(size(I), maxTraceIndex);
    Fc = F(FcIndex);
    Dc = D(DcIndex);
    F = F - Fc;
    D = D - Dc;

    mg = marginal(I);
    ac = autocorrelation(I);

    switch bitshift(int8(isInputN), 2) ...
                + bitshift(int8(~isempty(Dm)), 1) ...
                + bitshift(int8(~isempty(eps)), 0)
        case {0, 4} % N
            [~, Do] = FTconvert(D, maxN);
            Dm = -Do(1);
        case {1, 3, 5} %eps
            sp_eps = PulseMainWidth(F, mg, eps);
            it_eps = PulseMainWidth(D, ac, eps);

            if ~isempty(Dm)
                Dm = max(Dm, it_eps / 2);
            elseif isempty(it_eps)
                Dm = -D(1);
                if isempty(sp_eps)
                    sp_eps = -2 * F(1);
                end
            else
                Dm = it_eps / 2;
            end
            if ~isInputN
                minN = max(minN, 2 * sp_eps * Dm); % keep Fm = N/(4*Dm)> sp_eps /2
                maxN = max(minN, maxN);
            end
            Do = FTcutoff(D, Dm);
        case 2 % Dm
            Do = FTcutoff(D, Dm);
            N = maxN(min(length(Do), maxN), minN);
            maxN = N;
            minN = N;

        case 6
            Do = FTcutoff(D, Dm);
            N = maxN;
            % [~, Do] = FTconvert(Do, N);
    end

    spectrum_rms_width = PulseWidth(F, mg, 'rms') / sqrt(2);
    It_rms_width = PulseWidth(D, ac, 'rms') / sqrt(2);

    ac_factor = 1.54;
    spectrum_fwhm = PulseWidth(F, mg, 'fwhm') / ac_factor;
    It_fwhm = PulseWidth(D, ac, 'fwhm') / ac_factor;
    isp = spectrum_fwhm * It_fwhm; %Intensity and Spectrum product
    Mi = It_fwhm / (Do(2) - Do(1));

    dF = 1 / Dm / 2;
    M = spectrum_fwhm / dF;

    N = M^2 / isp;
    k = round(Mi / M);

    if N > maxN
        M = sqrt(isp * maxN);
        k = round(Mi / M);
        N = min(N / k^2,maxN);
    end
    if N < minN
        M = sqrt(isp * minN);
        k = round(Mi / M);
        N = minN;
    end
    if isPower2Required
        N = 2^round(log2(N));
    else
        N = floor(N / 2) * 2;
    end
    [Fo, Do] = FTconvert(Do, N, k);

    M1 = It_fwhm / (Do(2) - Do(1));
    M2 = spectrum_fwhm / (Fo(2) - Fo(1));
    if M > Mi
        warning("Frog Trace sampling rate is not enough!")
    end

    % I = I(:, ismember(D, Do)); %会由于存在数值误差而失效
    I = I(:, ismembertol(D, Do, 1e-8));
    Io = FreqTransfer(F, I, Fo);
    Io = Io ./ max(Io, [], "all");
    Io(Io < 0) = 0;
    if size(Io, 2) < N
        Io = [zeros(N, ceil((N - size(Io, 2)) / 2)), Io, zeros(N, floor((N - size(Io, 2)) / 2))];
    end
end
