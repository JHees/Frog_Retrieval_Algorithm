function plotPulse(T, P, Errorbar, plotRange_Time, plotRange_Freq, phaseEps, varargin)
    if nargin < 6
        phaseEps = 1e-1;
    end
    isPlotErrorbar = ~isempty(Errorbar);

    N = numel(P);
    [P, P_sp, f0] = removeFirstOrderPhase(P);
    % P_sp = fftshift(fft(P));
    % P = ifft(fft(P));

    It = abs(P).^2;
    spectrum = abs(P_sp).^2;
    It = It ./ max(It);
    spectrum = spectrum ./ max(spectrum);

    [F, T] = FTconvert(T);

    [~, max_ind] = max(It);
    [~, max_ind_sp] = max(spectrum);

    subplot(2, 1, 1)
    hold on
    if phaseEps == 0
        region = [T(1), T(end)];
    else
        [~, region] = PulseFullWidth(T, It, phaseEps);
    end
    ind = T >= region(1) & T <= region(2);
    ang = unwrap(angle(P));

    yyaxis left
    if isPlotErrorbar
        errorbar(T, It, Errorbar(:, 1), Errorbar(:, 2), varargin{:})
    else
        plot(T, It, varargin{:})
    end
    xlim([-plotRange_Time, plotRange_Time]);
    ylabel("Intensity (a.u.)")
    yyaxis right
    if isPlotErrorbar
        errorbar(T(ind), (ang(ind) - ang(max_ind)) / pi, Errorbar(ind, 3) / pi, Errorbar(ind, 4) / pi, varargin{:})
    else
        plot(T(ind), (ang(ind) - ang(max_ind)) / pi, varargin{:})
    end
    ylabel("Angle (pi rad)")
    xlabel("time (fs)")

    subplot(2, 1, 2)
    hold on
    ind_sp = spectrum > phaseEps * max(spectrum);
    ang_sp = unwrap(angle(P_sp));

    yyaxis left
    plot(F, spectrum)
    xlim([-plotRange_Freq, plotRange_Freq]);
    ylabel("Intensity (a.u.)")
    yyaxis right
    plot(F(ind_sp), (ang_sp(ind_sp) - ang_sp(max_ind_sp)) / pi)
    ylabel("Angle (pi rad)")
    xlabel("frequency (PHz)")
end
