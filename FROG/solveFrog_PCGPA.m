function [res] = solveFrog_PCGPA(I, initialGuess, iterMax, eps, method)

    if nargin > 4 && strcmpi(method, "svd")
        issvd = true;
    else
        issvd = false;
    end
    if isempty(initialGuess)
        if isscalar(iterMax)
            iterMax_init = min(iterMax / 5, 200);
        else
            iterMax_init = iterMax(1:end - 1);
            iterMax = iterMax(end);
        end

        init = solveFrog_GP(I, [], 1e-2, iterMax_init, 0e-2);
        P = init.P;
        G = init.P;
    else
        if isvector(initialGuess)
            P = initialGuess(:);
            G = initialGuess(:);
        else
            if size(initialGuess,2)==2
                P = initialGuess(:,1);
                G = initialGuess(:,2);
            else
                P = initialGuess(1,:).';
                G = initialGuess(2,:).';
            end
        end
        init.err = TraceError(I, TraceGenerate(P,G));
        init.iter = 1;
        iterMax = iterMax(end);
    end

    res.err = zeros(1, iterMax);
    i = 1;
    while i <= iterMax
        if issvd
            [P, G, res.err(i)] = principal_component_generlized_projection_svd(I, P, G, 1);
        else
            [P, G, res.err(i)] = principal_component_generlized_projection(I, P, G, 1);
        end
        if res.err(i) < eps
            break;
        end
        i = i + 1;
        % figure(3)
        % plot(abs(P))
        % hold on
        % plot(abs(G))
        % hold off
        % err(i)
    end
    res.P = P ./ max(abs(P));
    res.G = G ./ max(abs(G));

    % [res.P] = removeFirstOrderPhase(P);
    % [res.G] = removeFirstOrderPhase(G);
    res.iter = [init.iter, i - 1];
    res.err = [init.err, res.err(1:i - 1)];

end

function [sig, sig_sp] = data_constraint(I, sig_sp, b)
    if nargin == 2 || b == 1
        sig_sp = sqrt(I) .* exp(1i * angle(sig_sp));
    else
        sig_sp = sig_sp .* (sqrt(I) ./ abs(sig_sp)).^b;
        sig_sp(isnan(sig_sp)) = 0;
    end
    sig = ifft(ifftshift(sig_sp, 1), [], 1);
end

function [P, G, err] = principal_component_generlized_projection(I, P, G, b)
    N = size(I, 1);
    ind_outer = (1:N)' + mod((0:N:N^2 - 1) + (0:N:N^2 - 1)', N^2);
    ind_outer_invers = (1:N)' + mod((0:N:N^2 - 1) + [0; flip(N:N:N^2 - 1)'], N^2);

    outer = (P) .* (G.');
    % outer = outer + outer.'; % SHG
    shifting_outer = outer(ind_outer);
    sig = circshift(shifting_outer, N / 2, 2);
    sig_sp = circshift(fft(sig, [], 1), N / 2, 1);

    sig = data_constraint(I, sig_sp, b);
    shifting_outer = circshift(sig, N / 2, 2);
    outer = shifting_outer(ind_outer_invers);
    outer = outer + outer.'; % SHG

    P = outer * outer.' * P;
    G = outer.' * outer * G;

    P = P ./ max(abs(P).^2);
    G = G ./ max(abs(G).^2);

    S = TraceGenerate(P, G);
    err = TraceError(I, S);
end

function [P, G, err] = principal_component_generlized_projection_svd(I, P, G, b)
    N = size(I, 1);
    ind_outer = (1:N)' + mod((0:N:N^2 - 1) + (0:N:N^2 - 1)', N^2);
    ind_outer_invers = (1:N)' + mod((0:N:N^2 - 1) + [0; flip(N:N:N^2 - 1)'], N^2);

    outer = (P) .* (G.');
    % outer = outer + outer.'; % SHG
    shifting_outer = outer(ind_outer);
    sig = circshift(shifting_outer, N / 2, 2);
    sig_sp = circshift(fft(sig, [], 1), N / 2, 1);

    sig = data_constraint(I, sig_sp, b);
    shifting_outer = circshift(sig, N / 2, 2);
    outer = shifting_outer(ind_outer_invers);
    outer = outer + outer.'; % SHG
    
    [u, s, v] = svd(outer);
    P = u(:, 1) ./ sqrt(s(1, 1));
    G = conj(v(:, 1)) ./ sqrt(s(1, 1));

    S = TraceGenerate(P, G);
    err = TraceError(I, S);
end
