function Errorbar = TraceErrorbar(T,pulse, RetrievalFunc, BootStrapCount)
    N = length(pulse);
    P_retrieval = zeros(N, BootStrapCount);
    P_retrieval_sp = zeros(N, BootStrapCount);

    [pulse,pulse_sp]=PulseNormalized(T,pulse);
    for i = 1:BootStrapCount
        ind = randperm(N);
        ind = ind(1:N / 2);
        P = pulse;
        P(ind) = 0;
        P = RetrievalFunc(P);
        [P_retrieval(:, i), P_retrieval_sp(:, i)]=PulseNormalized(T,P);
    end
    Errorbar = zeros(N, 8);
    retrieval_intensity = abs(P_retrieval).^2;
    retrieval_phase = angle(P_retrieval) - angle(P_retrieval(N / 2,:));
    retrieval_sp_intensity = abs(P_retrieval_sp).^2;
    retrieval_sp_phase = angle(P_retrieval_sp) - angle(P_retrieval_sp(N / 2,:));

    pulse_intensity = abs(pulse).^2;
    pulse_phase = angle(pulse) - angle(pulse(N / 2));
    pulse_sp_intensity = abs(pulse_sp).^2;
    pulse_sp_phase = angle(pulse_sp) - angle(pulse_sp(N / 2));

    custom_mod = @(a, b) mod(a + b / 2, b) - b / 2;

    % TODO angle差值超+—过pi需要向下取模
    Errorbar(:, 1) = max(retrieval_intensity - pulse_intensity, [], 2);
    Errorbar(:, 2) = min(retrieval_intensity - pulse_intensity, [], 2);
    Errorbar(:, 3) = max(custom_mod(retrieval_phase - pulse_phase, pi), [], 2);
    Errorbar(:, 4) = min(custom_mod(retrieval_phase - pulse_phase, pi), [], 2);
    Errorbar(:, 5) = max(retrieval_sp_intensity - pulse_sp_intensity, [], 2);
    Errorbar(:, 6) = min(retrieval_sp_intensity - pulse_sp_intensity, [], 2);
    Errorbar(:, 7) = max(custom_mod(retrieval_sp_phase - pulse_sp_phase, pi), [], 2);
    Errorbar(:, 8) = min(custom_mod(retrieval_sp_phase - pulse_sp_phase, pi), [], 2);
end
