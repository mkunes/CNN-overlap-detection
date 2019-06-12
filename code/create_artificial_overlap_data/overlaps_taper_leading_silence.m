function [data2, LSilLen] = overlaps_taper_leading_silence(VAD2, data2, Fs)
% linearly taper silence at the start of an utterance, to avoid seams
% see also: overlaps_taper_trailing_silence
%
% MK, 2018-05
%
% 2018-08-02
%   Moved this into its own file

    LSilLen = [];
    for ii = 1:numel(VAD2)
        if VAD2(ii)
            LSilLen = ii - 1;
            break;
        end
    end
    if isempty(LSilLen)
        warning('No speech found');
        return;
    end
    
    LSilLen_s = LSilLen * Fs/100;
    if LSilLen > 0
        coeff = (1:LSilLen_s) ./ (LSilLen_s + 1);
        data2(1:LSilLen_s) = data2(1:LSilLen_s) .* coeff';
    end
end

