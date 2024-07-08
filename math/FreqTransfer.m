function [spectrum_out, freq_out] = FreqTransfer(freq_in, spectrum, freq_out)
    freq_in = freq_in(:);
    freq_out = freq_out(:);
    if isrow(spectrum)
        spectrum = spectrum(:);
    end

    [freq_in, ind] = sort(freq_in);
    spectrum = spectrum(ind, :);

    [freq_out, ind] = sort(freq_out);
    inverseIndFreqOut = zeros(size(freq_out));
    inverseIndFreqOut(ind) = 1:length(freq_out);

    % TODO 对重采样大小进行缩放 
    % 自动排除spectrum=0的区域
    n = round(size(spectrum, 1) / 20);
    s_min = mean(spectrum([1:n, end - n:end],:),1);
    spectrum = spectrum - s_min;

    spectrum = spectrum(freq_in >= freq_out(1) & freq_in <= freq_out(end),:);
    freq_in = freq_in(freq_in >= freq_out(1) & freq_in <= freq_out(end));
    

    freq_left = freq_out(freq_out < freq_in(1));
    freq_right = freq_out(freq_out > freq_in(end));
    freq_out = freq_out(freq_out >= freq_in(1) & freq_out <= freq_in(end));

    Fc = min(freq_in);
    freq_in = freq_in - Fc;
    freq_out = freq_out - Fc;

    [N,Ns] = size(spectrum);
    freq_middle = linspace(min(freq_in), max(freq_in), N).';
    t = FTconvert(freq_middle, N);
    
    M_in = freq_transfer_matrix(freq_middle, freq_in, t);
    M_out = freq_transfer_matrix(freq_middle, freq_out, t);

    spectrum_out = real(M_out * pinv(M_in) * spectrum);
    spectrum_out(spectrum_out<0)=0;

    spectrum_out = [zeros(length(freq_left),Ns); spectrum_out; zeros(length(freq_right),Ns)];
    spectrum_out = spectrum_out(inverseIndFreqOut,:);
    
    freq_out = [freq_left; freq_out; freq_right];
    freq_out = freq_out(inverseIndFreqOut) + Fc;

end

function M = freq_transfer_matrix(freq_in, freq_out, t)
    % make sure all input vector is columu vector
    N = length(freq_in);
    iFT = 1 / N * exp(-2i * pi * t .* freq_in.');
    FT = exp(2i * pi * t.' .* freq_out);
    M = FT * iFT;
end