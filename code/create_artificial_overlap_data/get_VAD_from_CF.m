function [VAD, spSegs] = get_VAD_from_CF(CF, minPauseLen, lenTh)
% get speech/nonspeech labels + speech segment times from the CF values
%   given by Param_KW.exe
%
% MK, 2018-05
%
% Changelog:
%
% v1.1   2019-06-11
%   - Moved this to its own file
%   - Fixed incorrect default value of minPauseLen, which caused some very
%       short pauses not to be removed (minPauseLen had effectively been zero)
% v1.0   2018-05
%   Initial version

if nargin == 1
    lenTh = 25; % min. length of silence/speech intervals (no. of frames) 
        % during the first pass of length-checking - longer continuous
        % intervals are left as-is, shorter ones are relabeled based on
        % the ratio of speech/silence in the surrounding area
    minPauseLen = 25; % min. length of pauses between speech segments
        % during the second pass 
        % (if a pause is shorter, the two segments are merged)
end

minTotalSpeech = 0.25;

spTh = 0.25;
d = 25;


nframes = numel(CF);

% --- median filter and thresholding ---
filtCF = medfilt1(CF,3);
P = prctile(filtCF, [5 95]);

CFTh = 0.4;
while CFTh > 0 
    B = filtCF > (P(1) + CFTh * (P(2) - P(1)));
    if sum(B) / nframes < minTotalSpeech
        CFTh = CFTh - 0.05;
    else
        break;
    end
end
% --- end thresh ---


% --- get the starts and ends of all continuous silence/speech regions,
%       regardless of length ----
ns = sum(diff(B) ~= 0) + 1; % total number of the regions
segs = zeros(ns,3);

isSp = B(1);
startF = 1;
jj = 1;
for ii = 2:nframes
    if B(ii) ~= isSp

        endF = ii;
        segs(jj, :) = [isSp, startF, endF-1];
        isSp = B(ii);
        startF = ii;
        jj = jj + 1;
    end
end
segs(jj, :) = [isSp, startF, nframes];
% 

% --- check the length of speech/silence intervals - 
%      - relabel the shorter ones based on neighbourhood ---
C = B;
if B(1)
    ii = 1;
else
    ii = 2; % exclude silence at the beginning
end

if B(end)
    maxI = size(segs,1);
else
    maxI = size(segs,1) - 1; % exclude silence at the end
end

while ii < maxI
    if segs(ii,3) - segs(ii,2) + 1 < lenTh
        L = max(1, segs(ii,2) - d);
        R = min(nframes, segs(ii,3) + d);

        % area to the right or left has 25% or more speech -> label
        %   as speech; otherwise label as silence
        newlabel = (sum(B(L:segs(ii,2))) / d > spTh) || (sum(B(segs(ii,3):R)) / d > spTh);

        if newlabel ~= segs(ii,1)
            segs(ii,1) = newlabel;
            C(segs(ii,2):segs(ii,3)) = segs(ii,1); 
            ii = ii + 1; % skip the next segment (it's connected to this one now)
        end
    end
    ii = ii + 1;
end
% --- end recheck ---


% --- remove all remaining pauses that are shorter than minPauseLen ---
spSegs = zeros(ns,2);

isSp = C(1);
startF = 1;
nsNew = 0;
for ii = 2:nframes
    if C(ii) ~= isSp
        if isSp == true
            endF = ii;
            if nsNew > 0 && (startF - spSegs(nsNew,2) < minPauseLen)
                spSegs(nsNew,2) = endF;
            else
                nsNew = nsNew + 1;
                spSegs(nsNew, :) = [startF, endF-1];
            end
            isSp = false;
        else
            startF = ii;
            isSp = true;
        end
    end
end
if isSp
    nsNew = nsNew + 1;
    spSegs(nsNew, :) = [startF, nframes];
end
spSegs = spSegs(1:nsNew, :);


VAD = false(size(CF));

for ii = 1:size(spSegs,1)
    VAD(spSegs(ii,1):spSegs(ii,2)) = true;
end
% --- end remove pauses ---
