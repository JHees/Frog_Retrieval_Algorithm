function mg = marginal(I)
    n = round(size(I, 2) / 20);
    mg = sum(I, 2);
    mg = mg - mean(mg([1:n, end - n:end]));
    mg = mg / max(mg);
end
