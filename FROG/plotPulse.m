function plotPulse(T, P, Errorbar, plotRange_Time, plotRange_Freq, phaseEps)
    if nargin < 6
        phaseEps = 1e-1;
    end
    isPlotErrorbar = ~isempty(Errorbar);
    P=P(:);
    T=T(:);
    N = numel(P);

    [F, T] = FTconvert(T);

    P_sp = fftshift(fft(P));
    [~,region]=PulseMainWidth(F,abs(P_sp),0.5);
    if(region(1)>region(2))
        P = P.*exp(2i*pi*F*T(end));
        P_sp = fftshift(fft(P));
    end
    t_center = sum(T.*abs(P).^2)/sum(abs(P).^2);
    P_sp_offset = P_sp.*exp(2i*pi*F*t_center);
    P_offset = ifft(ifftshift(P_sp_offset));
    f_center = sum(F.*abs(P_sp_offset).^2)/sum(abs(P_sp_offset).^2);
    P = P_offset .*exp(-2i*pi*T*f_center);


    It = abs(P).^2;
    [~, max_ind] = max(It);
    P = P.*exp(-1i*angle(P(max_ind)));
    P_sp = fftshift(fft(P));
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
