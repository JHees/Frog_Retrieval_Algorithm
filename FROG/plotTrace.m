function plotTrace(S, T, F)
clf
imagesc(T, F, S)
hold on
contour(T', F, S, [1e-4, 1e-4], 'g')
contour(T', F, S, [0.5, 0.5], 'r')
end
