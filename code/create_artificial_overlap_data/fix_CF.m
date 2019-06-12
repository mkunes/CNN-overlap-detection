function CF = fix_CF(CF, data, Fs, thresh, minLen_t)
% Fix the CF values obtained with (an older version of) Param_KW.exe
%---
% Note: a previous version of Param_KW.exe had an issue with long silence:
%   1. in areas with exact zeroes, there's a constant, non-zero value of CF
%   2. in areas with long, non-zero silence, CF ramps up to ~0.5 and
%       stays around that value (but is not constant)
% Since I still use those files, this fixes it by explicitly assigning
%   CF=0 to longer stretches where wav has low values (dflt. < 1e-3)
% (TODO: recalculate CF with the fixed version, then remove all uses of this)
%
% MK, 2018-05
%
% Changelog:
%
% v1.0.1 2019-06-11
%   - Moved this to its own file (I *really* should just replace the data...)
% v1.0   2018-05
%   - Initial version
    
    if nargin < 4
        thresh = 1e-3;
    end
    if nargin < 5
        minLen_t = 0.25;
    end
    
    minLen_s = minLen_t * Fs;
    
    isSil = data < thresh;
    
    silStart = 0;
    currIsSil = false;
    for ii = 1:numel(data)
        if ~currIsSil && isSil(ii)
            currIsSil = true;
            silStart = ii;
        end
        if currIsSil && ~isSil(ii)
            silLen = ii - silStart;
            if silLen >= minLen_s
                start_v = round(silStart * 100 / Fs);
                end_v = round(ii * 100 / Fs);
                if start_v < 1
                    start_v = 1;
                end
                if end_v > numel(CF)
                    start_v = numel(CF);
                end
                CF(start_v:end_v) = 0;
            end
                
            currIsSil = false;
        end
    end
    if currIsSil
        silLen = numel(data) + 1 - silStart;
        if silLen >= minLen_s
            start_v = round(silStart * 100 / Fs);
            end_v = round(ii * 100 / Fs);
            if start_v < 1
                start_v = 1;
            end
            if end_v > numel(CF)
                start_v = numel(CF);
            end
            CF(start_v:end_v) = 0;
        end
    end

end
