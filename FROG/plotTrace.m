function plotTrace(S, T, F)
    if isscalar(T)
        ND = size(S, 2);
        T = (-ND / 2:ND / 2 - 1)' .* T;
    end
    imagesc(T, F, S)
    hold on
    % contour(T', F, S, [1e-4, 1e-4], 'g')
    % contour(T', F, S, [0.5, 0.5], 'r')
    xlabel("Delay (fs)")
    ylabel("Frequency (PHz)")
end
