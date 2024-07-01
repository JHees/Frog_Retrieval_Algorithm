function mg = marginal(I)
    mg = sum(I, 2);
    mg = mg - median(mg);
    mg = mg / max(mg);
end
