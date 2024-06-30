function [res, Errorbar] = TraceSolver(I, initialGuess, iters, loop, offspringNum, eps)
    % TraceSolver. This function solves for an optimized trace using evolutionary algorithms, leveraging multiple iteration strategies.
    % 此函数使用多种迭代策略的进化算法来求解优化的trace。
    % The function iteratively improves the trace estimate using a series of solution strategies including vanilla, generalized projection,
    % and principal component methods, each with specific adaptations for enhanced performance.
    % Usage:
    %     [res, Errorbar] = TraceSolver(I, initialGuess, iters, loop, offspringNum, eps)
    %
    % Input:
    %     I: Input trace data (size: Fn x Dn)
    %     initialGuess: Initial guess for the solution (default: computed if empty)
    %     iters: Array with number of iterations per method stage
    %     loop: Number of loops per iteration phase
    %     offspringNum: Number of offspring solutions to generate per iteration
    %     eps: Convergence threshold
    %
    % Output:
    %     res: Struct containing the optimized trace, error metrics, and iteration count
    %     Errorbar: Error metrics for solution robustness analysis (optional)
    %
    % AUTHOR: Huang Xingzhao, June 30, 2024

    N = size(I, 2);

    if isempty(initialGuess)
        x = (-N / 2:N / 2 - 1)';
        x_tou = N / 16;
        P = exp(- x.^2/2 ./ x_tou^2) + 0.5 * exp(- (x + N / 8).^2/2 ./ x_tou^2);
        % init = solveFrog_vanilla(I, [], 50, eps);
        % P = init.P;
    else
        P = initialGuess(:);
    end
    total_iter = iters(1);
    res.err = zeros(1, ceil(total_iter / loop) * loop);
    err_count = 0;
    err = TraceError(I, P);
    sols = { ...
                @vanilla, ...
                @generalized_projection, ...
                @generalized_projection_shortcut, ...
                @principal_component_generlized_projection, ...
                @principal_component_generlized_projection_svd, ...
            };

    P_restore = P;

    for i = 1:ceil(total_iter / loop)
        sols_var = { ...
                        {min(1, (1.1)^(1 + i / 5))}, ...
                        {min(1, (1.1)^(1 + i / 5)), 1e-2}, ...
                        {min(1, (1.1)^(1 + i / 5)), 1e-3, 3}, ...
                        {"G", min(1, (1.1)^(1 + i / 5))}, ...
                        {"G", min(1, (1.1)^(1 + i / 5))}, ...
                    };
        [P_restore, err_restore] = Solver(I, P_restore, loop, sols, sols_var);
        P_restore = P_restore(:, 1:min(size(P_restore, 2), offspringNum));

        if err_restore(end, 1) < err
            P = P_restore(:, 1);
            err_count = err_count + 1;
            res.err((err_count - 1) * loop + 1:err_count * loop) = err_restore(:, 1)';
            err = err_restore(end, end);
            if err < eps
                break;
            end
        end
    end
    P = P ./ max(abs(P));
    res.P = P;
    res.iter = err_count * loop;
    res.err = [res.err(1:err_count * loop)];
    % pcgp_svd
    if numel(iters) >= 2
        pcgp_svd_iter = iters(2);
        err_svd = zeros(1, pcgp_svd_iter);
        G = P;
        for i = 1:pcgp_svd_iter
            [P, G] = principal_component_generlized_projection_svd(I, P, G, 1);
            err_svd(i) = TraceError(I, P);
            if err < eps
                break;
            end
        end
        res.P = P ./ max(abs(P));
        res.iter = err_count * loop + i;
        res.err = [res.err, err_svd(1:i)];
    end
    if nargout == 2 % error bar
        boot_strap = 10;
        P_bs = zeros(N, boot_strap);
        Errorbar = zeros(N, 4);
        for j = 1:boot_strap
            ind = randperm(N);
            ind = ind(1:N / 2);
            P = res.P;
            P(ind) = 0;
            G = P;
            for i = 1:pcgp_svd_iter
                [P, G] = principal_component_generlized_projection_svd(I, P, G, 1);
            end
            [P, ~, ~] = removeFirstOrderPhase(P);
            P_bs(:, j) = P;
        end
        Errorbar(:, 1) = max(abs(P_bs).^2 - abs(P).^2, [], 2);
        Errorbar(:, 2) = min(abs(P_bs).^2 - abs(P).^2, [], 2);
        Errorbar(:, 3) = max(angle(P_bs) - angle(P), [], 2);
        Errorbar(:, 4) = min(angle(P_bs) - angle(P), [], 2);
    end

end
function [P_out, err_out] = Solver(I, P_in, iter, solve_algorithm, varargin)
    solNum = numel(solve_algorithm);
    [N, P_num] = size(P_in);
    P_out = zeros(N, solNum * P_num);
    err_out = zeros(iter, solNum * P_num);

    for k = 1:solNum * P_num
        [i, j] = ind2sub([P_num, solNum], k);
        P = P_in(:, i);
        if strcmpi(varargin{:}{j}{1}, "G")
            G = P_in(:, i);
            for ii = 1:iter
                var = varargin{:}{j};
                var{1} = G;
                [P, G] = solve_algorithm{j}(I, P, var{:});
                err_out(ii, k) = TraceError(I, P);
            end
        else
            for ii = 1:iter
                [P, err_out(ii, k)] = solve_algorithm{j}(I, P, varargin{:}{j}{:});
                % err_out(ii, k) = TraceError(I, P);
            end
        end
        P_out(:, k) = P;
    end
    [~, sort_ind] = sort(err_out(end, :));
    err_out = err_out(:, sort_ind);
    P_out = P_out(:, sort_ind);
end
function [P, err] = vanilla(I, P, b)
    [~, sig_sp] = TraceGenerate(P);
    sig = data_constraint(I, sig_sp, b);
    P = sum(sig, 2);
    P = P ./ max(abs(P));
    err = TraceError(I, P);
end

function [P, err] = generalized_projection(I, P, b, beta)
    persistent ind_plus ind_minus ind_sig
    if isempty(ind_plus) || isempty(ind_minus) || isempty(ind_sig)
        N = size(I, 1);
        ind_plus = mod((-N / 2:N / 2 - 1) + (0:N - 1)', N) + 1;
        ind_minus = flipud(mod(-ind_plus, N + 1));
        ind_sig = ind_minus + (0:N - 1) * N;
        % mask_minus = ones(N) - triu(ones(N), N / 2 + 1) - tril(ones(N), -2);
        % mask_plus = flipud(mask_minus);
    end

    [~, sig_sp] = TraceGenerate(P);
    sig = data_constraint(I, sig_sp, b);

    P_mat_plus = P(ind_plus);
    P_mat_minus = P(ind_minus);
    sig_mat = sig(ind_sig);

    gradient = (-conj(sig) .* P_mat_plus + P_mat_plus .* conj(P_mat_plus) .* conj(P)) ...
        + (-conj(sig_mat) .* P_mat_minus + P_mat_minus .* conj(P_mat_minus) .* conj(P));
    gradient = sum(gradient, 2);
    % Z = norm(sig - P .* P_mat_plus, 'fro');
    P = P - 2 * beta * conj(gradient);
    P = P ./ max(abs(P));
    err = TraceError(I, P);
end

function [P, err] = generalized_projection_shortcut(I, P, b, beta, shortcutStepSize)

    [~, Last.sig_sp] = TraceGenerate(P);
    [~, This.sig_sp] = data_constraint(I, Last.sig_sp, b);
    Last.sig_sp_middle = (Last.sig_sp + This.sig_sp) / 2;

    [P, ~] = generalized_projection(I, P, b, beta);

    [~, Last.sig_sp] = TraceGenerate(P);
    [~, This.sig_sp] = data_constraint(I, Last.sig_sp, b);
    This.sig_sp_middle = (Last.sig_sp + This.sig_sp) / 2;

    sig_sp = shortcutStepSize * (This.sig_sp_middle - Last.sig_sp_middle) + This.sig_sp_middle;
    sig = ifft(ifftshift(sig_sp, 1), [], 1);
    P = sum(sig, 2);
    P = P ./ max(abs(P));
    err = TraceError(I, P);

end

function [P, G, err] = principal_component_generlized_projection(I, P, G, b)
    N = size(I, 1);
    persistent ind_outer ind_outer_invers
    if isempty(ind_outer) || isempty(ind_outer_invers)
        ind_outer = (1:N)' + mod((0:N:N^2 - 1) + (0:N:N^2 - 1)', N^2);
        ind_outer_invers = (1:N)' + mod((0:N:N^2 - 1) + [0; flip(N:N:N^2 - 1)'], N^2);
    end

    outer = (P) .* (G.');
    outer = outer + outer.'; % SHG
    shifting_outer = outer(ind_outer);
    sig = circshift(shifting_outer, N / 2, 2);
    sig_sp = circshift(fft(sig, [], 1), N / 2, 1);

    sig = data_constraint(I, sig_sp, b);
    shifting_outer = circshift(sig, N / 2, 2);
    outer = shifting_outer(ind_outer_invers);

    P = outer * outer.' * P;
    G = outer.' * outer * G;

    P = P ./ max(abs(P));
    G = G ./ max(abs(G));

    S = TraceGenerate(P, G);
    err = TraceError(I, S);
end

function [P, G, err] = principal_component_generlized_projection_svd(I, P, G, b)
    N = size(I, 1);
    persistent ind_outer ind_outer_invers
    if isempty(ind_outer) || isempty(ind_outer_invers)
        ind_outer = (1:N)' + mod((0:N:N^2 - 1) + (0:N:N^2 - 1)', N^2);
        ind_outer_invers = (1:N)' + mod((0:N:N^2 - 1) + [0; flip(N:N:N^2 - 1)'], N^2);
    end

    outer = (P) .* (G.');
    outer = outer + outer.'; % SHG
    shifting_outer = outer(ind_outer);
    sig = circshift(shifting_outer, N / 2, 2);
    sig_sp = circshift(fft(sig, [], 1), N / 2, 1);

    sig = data_constraint(I, sig_sp, b);
    shifting_outer = circshift(sig, N / 2, 2);
    outer = shifting_outer(ind_outer_invers);

    [s, u, v] = svd(outer);
    P = s(:, 1) ./ sqrt(u(1, 1));
    G = conj(v(:, 1)) ./ sqrt(u(1, 1));

    S = TraceGenerate(P, G);
    err = TraceError(I, S);
end

function [sig, sig_sp] = data_constraint(I, sig_sp, b)
    % if nargin == 2 || b == 1
    %     sig_sp = sqrt(I) .* exp(1i * angle(sig_sp));
    % else
    sig_sp = sig_sp .* (sqrt(I) ./ abs(sig_sp)).^b;
    sig_sp(isnan(sig_sp)) = 0;
    % end
    sig = ifft(ifftshift(sig_sp, 1), [], 1);
end
