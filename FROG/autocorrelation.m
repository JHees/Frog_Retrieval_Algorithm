function ac = autocorrelation(I)
    n = round(size(I, 1) / 20);
    ac = sum(I, 1);
    ac = ac - mean(ac([1:n, end - n:end]));
    ac = ac / max(ac);
end
