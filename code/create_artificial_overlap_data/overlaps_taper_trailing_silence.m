function [data1, TSilLen] = overlaps_taper_trailing_silence(VAD1, data1, Fs)
% linearly taper silence at the end of an utterance, to avoid seams
% see also: overlaps_taper_leading_silence
%
% MK, 2018-05
%
% 2018-08-02
%   Moved this into its own file
    
    TSilLen = [];
    for ii = numel(VAD1):-1:1
        if VAD1(ii)
            TSilLen = numel(VAD1) - ii;
            break;
        end
    end
    if isempty(TSilLen)
        warning('No speech found');
        return;
    end

    TSilLen_s = TSilLen * Fs/100;
    if TSilLen > 0
        coeff = (TSilLen_s:-1:1) ./ (TSilLen_s + 1);
        data1((end-TSilLen_s+1):end) = data1((end-TSilLen_s+1):end) .* coeff';
    end
end
