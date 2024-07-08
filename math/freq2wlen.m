function [wlen_out, spectrum_out] = freq2wlen(freq, spectrum, f0)
    % 必须要求单位制为 fs,nm,PHz
    c = 299.792458;
    if nargin <= 2 || isempty(f0)
        f0 = 0;
    end
    if isrow(spectrum)
        spectrum = spectrum(:);
    end
    freq_in = freq(:) + f0;

    m = max(spectrum, [], 1);
    spectrum = freq_in.^2 ./ c .* spectrum;
    spectrum = spectrum ./ max(spectrum) .* m;
    N = size(spectrum, 1);

    wlen_in = c ./ freq_in;
    wlen_out = linspace(min(wlen_in), max(wlen_in), N).';
    freq_out = c ./ wlen_out;

    spectrum_out = FreqTransfer(freq_in, spectrum, freq_out);
    spectrum_out(spectrum_out < 0) = 0;
    spectrum_out = spectrum_out ./ max(spectrum_out) .* max(spectrum);

end
