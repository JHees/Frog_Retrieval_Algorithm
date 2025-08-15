function Io = TraceResampleAgain(I, D, F, Do, Fo, Dc, Fc)
    [NF, ND] = size(I);
    if isscalar(D)
        D = (-ND / 2:ND / 2 - 1)' .* D - Dc;
    end
    I = I(:, ismembertol(D, Do, 1e-8));
    Io = FreqTransfer(F - Fc, I, Fo);
    Io = Io ./ max(Io, [], "all");
    Io(Io < 0) = 0;
    N = size(Io, 1);
    if size(Io, 2) < N
        Io = [zeros(N, ceil((N - size(Io, 2)) / 2)), Io, zeros(N, floor((N - size(Io, 2)) / 2))];
    end
end
