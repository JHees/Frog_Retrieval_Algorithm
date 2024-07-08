function [freq_out, spectrum_out] = wlen2freq(wlen, spectrum)
    % 必须要求单位制为 fs,nm,PHz
    c = 299.792458;
    wlen = wlen(:);
    if isrow(spectrum)
        spectrum = spectrum(:);
    end
    m = max(spectrum, [], 1);
    spectrum = wlen.^2 ./ c .* spectrum;
    spectrum = spectrum ./ max(spectrum) .* m;
    N = size(spectrum, 1);

    freq_in = c ./ wlen;
    freq_out = linspace(min(freq_in), max(freq_in), N).';

    spectrum_out = FreqTransfer(freq_in, spectrum, freq_out);
    spectrum_out(spectrum_out < 0) = 0;
    spectrum_out = spectrum_out ./ max(spectrum_out) .* max(spectrum);
end
