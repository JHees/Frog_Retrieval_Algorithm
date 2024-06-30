function err = TraceError(I, S)
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
