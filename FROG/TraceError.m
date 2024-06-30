function err = TraceError(I, S)
    % TraceError. This function calculates the Frobenius norm of the error between the measured and simulated traces.
    % 此功能计算测量和模拟trace之间的误差的Frobenius范数。
    % It evaluates the accuracy of a trace simulation relative to actual measured data by calculating a normalized error.
    % Usage:
    %     err = TraceError(I, S)
    %
    % Input:
    %     I: Measured trace (size: N x N)
    %     S: Simulated or reconstructed trace (size: N x N or vector form)
    %
    % Output:
    %     err: Normalized error (scalar)
    %
    % Note:
    %     The function adjusts the scale of the simulated trace S to best fit the measured trace I,
    %     aiming to minimize the error between S and I. It then calculates the normalized Frobenius norm
    %     of the residual error, providing a measure of fit quality.
    %
    % AUTHOR: Huang Xingzhao

    N = size(I, 1);
    if isrow(S) || iscolumn(S)
        S = TraceGenerate(S);
    end
    x = sum(I .* S, 'all') / sum(S .* S, 'all');
    % cvx_begin
    %     variable x
    %     minimize norm(I - x * S,'fro');
    % cvx_end
    err = norm(I - x * S, 'fro') / N;
end
