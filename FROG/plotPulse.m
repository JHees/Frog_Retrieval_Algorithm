function plotPulse(T, P, Errorbar, plotRange_Time, plotRange_Freq, phaseEps)
    if nargin < 6
        phaseEps = 1e-1;
    end
    isPlotErrorbar = ~isempty(Errorbar);
    P=P(:);
    T=T(:);
    N = numel(P);

    [F, T] = FTconvert(T);

    [P,P_sp]=PulseNormalized(T,P);
    It = abs(P).^2;
    [~, max_ind] = max(It);
    spectrum = abs(P_sp).^2;
    [~, max_ind_sp] = max(spectrum);
   
    subplot(2, 1, 1)
    hold("on")
    if phaseEps == 0
        region = [T(1), T(end)];
    else
        [~, region] = PulseFullWidth(T, It, phaseEps);
    end
    ind = T >= region(1) & T <= region(2);
    ang = unwrap(angle(P));

    yyaxis left
    if isPlotErrorbar
        errorbar(T, It, Errorbar(:, 1), Errorbar(:, 2))
    else
        plot(T, It,'.-')
    end
    xlim([-plotRange_Time, plotRange_Time]);
    ylabel("Intensity (a.u.)")
    yyaxis right
    if isPlotErrorbar
        errorbar(T(ind), (ang(ind) - ang(max_ind)) / pi, Errorbar(ind, 3) / pi, Errorbar(ind, 4) / pi)
    else
        plot(T(ind), (ang(ind) - ang(max_ind)) / pi,'.-')
    end
    ylabel("Angle (pi rad)")
    xlabel("time (fs)")

    subplot(2, 1, 2)
    hold("on")
    ind_sp = spectrum > phaseEps * max(spectrum);
    ang_sp = unwrap(angle(P_sp));

    yyaxis left
    if isPlotErrorbar
        errorbar(F, spectrum, Errorbar(:, 5), Errorbar(:, 6))
    else
        plot(F, spectrum,'.-')
    end
    xlim([-plotRange_Freq, plotRange_Freq]);
    ylabel("Intensity (a.u.)")
    yyaxis right
    if isPlotErrorbar
        errorbar(F(ind_sp), (ang_sp(ind_sp) - ang_sp(max_ind_sp)) / pi, Errorbar(ind_sp, 7) / pi, Errorbar(ind_sp, 8) / pi)
    else
        plot(F(ind_sp), (ang_sp(ind_sp) - ang_sp(max_ind_sp)) / pi,'.-')
    end
    ylabel("Angle (pi rad)")
    xlabel("frequency (PHz)")
end
