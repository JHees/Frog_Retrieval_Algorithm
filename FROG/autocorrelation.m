function ac = autocorrelation(I)
    ac = sum(I, 1);
    ac = ac - median(ac);
    ac = ac / max(ac);
end