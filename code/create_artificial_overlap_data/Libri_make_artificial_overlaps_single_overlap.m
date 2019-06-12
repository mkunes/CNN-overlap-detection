function Libri_make_artificial_overlaps_single_overlap(wav_parentdir, dir_prm, dir_out_main, targetOLType)
% Libri_make_artificial_overlaps_single_overlap(wav_parentdir, dir_prm, dir_out_main, targetOLType)
%
% creates short artificial overlap data from LibriSpeech:
%   a) two utterances with a short overlap of up to 2s - TYPE_OL_SHORT
%   b) two utterances with a pause of up to 1s - TYPE_PAUSE
%   c) two utterances with a long overlap (third of the total length) - TYPE_OL_LONG
%   d) one utterance with a single word/short phrase (0.25-2s) inserted
%       between words, over speech, or in a random location - TYPE_INSERT
% with added office/meeting/hallway noise and with or without reverberation
%
% MK, 2018-05
%
% Changelog:
% 1.3  2019-06-11
%   Fixed a bug when calling get_VAD_from_CF - due to wrong values of minPauseLen,
%		some very short pauses were not getting removed (minPauseLen was effectively zero)
%   Moved get_VAD_from_CF and fix_CF to their own m-files
% 1.2  2019-05-02
%   Replaced some hardcoded settings with input parameters
%   	(wav_parentdir, dir_prm, dir_out_main, targetOLType)
%   Changed name format of output files to include the type of OL, noise and reverb
%       (originally: 'filename1-filename2'
%        now:        'filename1-filename2-OLType.noisetype[.reverbtype]')
% 1.13 2018-08-02
%   Moved the following functions into separate files:
%       overlaps_taper_trailing_silence
%       overlaps_taper_leading_silence
%       overlaps_add_DEMAND_noise (previously named "add_DEMAND_noise")
%       overlaps_add_AIR_reverb (previously named "add_AIR_reverb")
% 1.12 2018-07-31
%   Modified add_DEMAND_noise:
%    better random start 
%       (originally - concatenate, then choose random start within length constraints; 
%        now - random start, then concatenate)
%    seamless transition between noise files
% 1.11 2018-07-18
%   Added an option to select only one specific type of overlap
% 1.1 2018-05
%   Added TYPE_OL_LONG, TYPE_INSERT subtypes (in-speech, in-pause)
% 1.0 2018-05
%   Original version:
%       TYPE_OL_SHORT, TYPE_PAUSE, TYPE_INSERT - random only
%       noise + reverb
%
% Args:
%   wav_parentdir - main LibriSpeech directory containing individual speakers
%                   (e.g. "LibriSpeech/train-other-500")
%   dir_prm - directory with .prm files from the Param_KW feature extractor
%   dir_out_main - main output directory (script will create subdirs based
%           on overlap type)
%   targetOLType - optional, specifies the type of overlaps to be created;
%           if not given or zero, it will create a mix of all types
%
% -----
% TODO:
%   - Make this compatible with some kind of publicly available VAD, 
%       not just the one in Param_KW
%   - Move all paths to a config file


TYPE_ALL = 0;
TYPE_INSERT = 1;
TYPE_PAUSE = 2;
TYPE_OL_SHORT = 3;
TYPE_OL_LONG = 4;

STRING_OL_LONG = 'overlap_big';
STRING_OL_SHORT = 'overlap_small';
STRING_PAUSE = 'pause';
STRING_INSERT_RANDOM = 'insert_random';
STRING_INSERT_INPAUSE = 'insert_inPause';
STRING_INSERT_INSPEECH = 'insert_inSpeech';
SUFFIX_REVERB = '_reverb';


makePlots = false;


mypath_datasets = getSystemSpecificPath('datasets');

dir_DEMAND_noise = [mypath_datasets '/Demand_DB_noise'];
dir_AIR_reverb = [mypath_datasets '/AIR_1_4/for_AMI'];

if nargin == 0
    wav_parentdir = [mypath_datasets '/LibriSpeech/train-other-500'];
    dir_out_main = [mypath_datasets '/overlaps/LibriSpeech/overlaps_single_overlap_short_other500_v2/'];
    dir_prm = [mypath_datasets '/LibriSpeech/train-other-500-prm'];
end

if nargin < 4
    targetOLType = TYPE_ALL; %TYPE_ALL; % which type of overlap should be created:
    % either a single specific type (TYPE_INSERT, TYPE_PAUSE, TYPE_OL_SHORT, TYPE_OL_LONG) 
    % or a mix of all types (TYPE_ALL)
    % (TYPE_INSERT still creates a mix of "inserted into speech", "inserted
    %   into pause" and "inserted randomly")
end

dir_out_OL_big = [dir_out_main '/' STRING_OL_LONG];
dir_out_OL_small = [dir_out_main '/' STRING_OL_SHORT];
dir_out_pause = [dir_out_main '/' STRING_PAUSE];
dir_out_insert_random = [dir_out_main '/' STRING_INSERT_RANDOM];
dir_out_insert_pause = [dir_out_main '/' STRING_INSERT_INPAUSE];
dir_out_insert_overlap = [dir_out_main '/' STRING_INSERT_INSPEECH];
% reverb dirs are the above plus SUFFIX_REVERB (e.g. ".../pause_reverb")

spkDirs = dir(wav_parentdir);
nspk = numel(spkDirs);


nOutFiles = 0;
maxOutFiles = 10000;

xmin = 0.3;
xmax = 0.7;

num_OL_big = 0;
num_OL_small =0;
num_pause = 0;
num_insert = 0;
num_insert_random = 0;
num_insert_overlap = 0;
num_insert_pause = 0;

% min and max number of seconds on each side of the overlap/pause
minSoloSpkLen = 2; % affects the initial location of the overlapping area
maxSoloSpkLen = 4; % final wav will be cropped to max. this many seconds on each side


% min and max length of an inserted ("backchannel") speech segment 
minInsertLen = 0.25;
maxInsertLen = 2;



nFailedInRow = 0; % number of failed attempts in a row

while nOutFiles < maxOutFiles
    
    if nFailedInRow > 20
        error('Too many failed attempts in a row - something''s probably wrong');
    end
    
    OLType = 0;
    
    
    % -----
    % Select 1st file
    % ----
    ntries = 1;
    spk_id1 = [];
    
    
    while isempty(spk_id1) 
        r = randi(nspk); % select a random folder
        
        % check if it is a valid speaker folder
        if spkDirs(r).isdir && spkDirs(r).name(1) ~= '.'
            spk_id1 = spkDirs(r).name;
        else
            ntries = ntries + 1;
            if ntries > 100 % max 100 attempts - otherwise something is wrong
                error('Could not find a valid speaker');
            end
        end
    end
    
    % select a random file from all of the speaker's books
    flacFiles = dir([wav_parentdir '/' spk_id1 '/*/*.flac']);
    nfiles = numel(flacFiles);
    
    r = randi(nfiles); % select a random file
    
    filepath1 = flacFiles(r).folder;
    filename1 = flacFiles(r).name;
    

    % -----
    % Select 2nd file
    % ----
    ntries = 1;
    spk_id2 = [];
    
    while isempty(spk_id2)
        r = randi(nspk); % select a random folder
        
        % check if it is a valid speaker folder, and different from the 1st
        if spkDirs(r).isdir && spkDirs(r).name(1) ~= '.' && ...
                ~strcmp(spkDirs(r).name, spk_id1)
            spk_id2 = spkDirs(r).name;
        else
            ntries = ntries + 1;
            if ntries > 100 % max 100 attempts - otherwise something is wrong
                error('Could not find a valid speaker');
            end
        end
    end

    % select a random file from all of the speaker's books
    flacFiles = dir([wav_parentdir '/' spk_id2 '/*/*.flac']);
    nfiles = numel(flacFiles);
    
    r = randi(nfiles); % select a random file
    
    filepath2 = flacFiles(r).folder;
    filename2 = flacFiles(r).name;
    
    
    
    % ----
    % Load wavs and CFs, get VAD
    % ----
    
    [~,basename1,~] = fileparts(filename1);
    [~,basename2,~] = fileparts(filename2);

    [data1, Fs] = audioread([filepath1 filesep filename1]);
    [data2, Fs2] = audioread([filepath2 filesep filename2]);
    
    if Fs ~= Fs2
        error('sample rate must be identical for all files');
    end

    [~,CF1] = ReadPRM([dir_prm '/' basename1 '.prm']);
    [~,CF2] = ReadPRM([dir_prm '/' basename2 '.prm']);

    % fixes incorrect CF in long silences
    % (caused by a bug in an older version of ParamKW, 
    %    bug is fixed now, but I still use those data)
    CF1 = fix_CF(CF1, data1, Fs, 1e-3, 0.25);
    CF2 = fix_CF(CF2, data2, Fs, 1e-3, 0.25);
    
    
    VAD1 = get_VAD_from_CF(CF1);
    VAD2 = get_VAD_from_CF(CF2);
    
    nd1 = numel(data1);
    nd2 = numel(data2);
    
    % ---
    % targetOLType:
    %   TYPE_INSERT:
    %     Try to find a short speech segment in the 2nd file that can be inserted into the first file
    %   TYPE_ALL:
    %     Same, unless there are already more such combinations than the other types
    %   other:
    %     Skip this
    % ---
    
    wordIdx = 0;
    
    if targetOLType == TYPE_INSERT || ...
            (targetOLType == TYPE_ALL && num_insert <= num_OL_small && ...
            num_insert <= num_OL_big && num_insert <= num_pause)
            
        [VAD2_b, spSegs] = get_VAD_from_CF(CF2, 10, 10);
        
        if size(spSegs, 1) > 2
            % go through the speech segments in a random order 
            %  and pick the first one with a correct length
            % (the first and the last segment are excluded)
            
            idx = 1 + randperm(size(spSegs, 1)-2); 

            minlen = 100 * minInsertLen;

            % make maxlen such that there will be enough space left around the inserted segment
            maxlen = min(100 * maxInsertLen, numel(CF1) - 2 * minSoloSpkLen * 100);

            for ii = idx
                seglen = spSegs(ii,2) - spSegs(ii,1);
                if seglen >= minlen && seglen <= maxlen
                    wordIdx = ii;
                    wordLen = seglen;
                    break;
                end
            end
        end
    end
    
    if targetOLType == TYPE_INSERT && wordIdx <= 0
        % Failed to find a good location to insert segment
        nFailedInRow = nFailedInRow + 1;
        continue;
    end
    
    
    % ---
    % Calculate overlap/pause position and taper silence at the boundaries
    % ---
    
    if wordIdx > 0 % segment insertion
        
        OLType = TYPE_INSERT;
        
        [~,minIdx] = min(CF2(spSegs(wordIdx-1,2):spSegs(wordIdx,1)));
        segLPauseLen_v = spSegs(wordIdx,1) - spSegs(wordIdx-1,2) - minIdx + 2;
        segStart_v = minIdx + spSegs(wordIdx-1,2) - 1;
        segStart_s = segStart_v * Fs / 100;

        [~,minIdx] = min(CF2(spSegs(wordIdx,2):spSegs(wordIdx+1,1)));
        segEnd_v = minIdx + spSegs(wordIdx,2) - 1;
        segEnd_s = segEnd_v * Fs / 100;

        data2 = data2(segStart_s:segEnd_s);
        VAD2 = VAD2_b(segStart_v:segEnd_v);
        CF2 = CF2(segStart_v:segEnd_v);
        nd2 = numel(data2);
        
        if numel(data2) > numel(data1)
            % first file is way too short
            nFailedInRow = nFailedInRow + 1;
            continue;
        end
            

        % Taper the leading and trailing silence
        [data2, TSilLen] = overlaps_taper_trailing_silence(VAD2, data2, Fs);
        [data2, LSilLen] = overlaps_taper_leading_silence(VAD2, data2, Fs);
        if isempty(TSilLen) || isempty(LSilLen)
            nFailedInRow = nFailedInRow + 1;
            continue;
        end

        num_insert = num_insert + 1;
        
        % ---
        % Try to insert the speech segment into the first file so that it
        %   a) fully overlaps with speech, or
        %   b) is inserted into a pause, with no overlap, or
        %   c) insert it randomly
        % If more than one option is possible, pick the least represented one
        % ---
        
        [~, spSegs] = get_VAD_from_CF(CF1, 10, 25);
        start_OL_v = 0;
        start_pause_v = 0;

        % pick a random speech region that is longer than the segment to be inserted
        %   and that isn't too early in the signal - so that there will 
        %   still be at least minSoloSpkLen seconds BEFORE the overlap
        idx = randperm(size(spSegs, 1));
        for ii = idx
            seglen = spSegs(ii,2) - spSegs(ii,1);
            if seglen >= wordLen && spSegs(ii,1) - segLPauseLen_v >= minSoloSpkLen * 100
                
                % random start point within the speech region (such that it
                %   still fits inside), offset by the leading pause in the inserted segment
                randMax = seglen - wordLen + 1;
                start_OL_v = spSegs(ii,1) + randi(randMax) - segLPauseLen_v;
                end_OL_v = start_OL_v + numel(CF2) - 1;
                if start_OL_v > 0 && end_OL_v <= numel(CF1)
                    break;
                else
                    start_OL_v = 0;
                end
            end
        end
        
        % pick a random pause that is longer than the segment
        %   and that isn't too early in the signal - so that there will 
        %   still be at least minSoloSpkLen seconds BEFORE the overlap
        if size(spSegs, 1) > 1
            idx = randperm(size(spSegs, 1) - 1);
            for ii = idx
                pauselen = spSegs(ii+1,1) - spSegs(ii,2);
                if pauselen >= wordLen && spSegs(ii,2) - segLPauseLen_v >= minSoloSpkLen * 100
                    
                    % random start point within the pause (such that it
                    %   still fits inside), offset by the leading pause in the inserted segment
                    randMax = pauselen - wordLen + 1;
                    start_pause_v = spSegs(ii,2) + randi(randMax) - segLPauseLen_v - 1;
                    end_pause_v = start_pause_v + numel(CF2) - 1;
                    if start_pause_v > 0 && end_pause_v <= numel(CF1)
                        break;
                    else
                        start_pause_v = 0;
                    end
                end
            end
        end

        if (start_OL_v <= 0 && start_pause_v <= 0) || ...
                (num_insert_pause > num_insert_random && ...
                 num_insert_overlap > num_insert_random)
             % no long pauses/speech found, OR there are not enough random
             %      placements -> place it randomly
             %      (but, again, no earlier than minSoloSpkLen seconds after start)
             
             % find a random place to insert the segment
            minSoloSpkLen_s = minSoloSpkLen * Fs;
            randMax = nd1 - (2 * minSoloSpkLen_s + nd2);
            
            if randMax < 0
                % first file is short (but the segment still fits inside, 
                %   that's checked elsewhere) 
                % -> ignore the minlen requirement and allow
                %   placement near the edge
                randMax = nd1 - nd2;
                minSoloSpkLen_s = 0;
            end

            start_s = minSoloSpkLen_s + randi(randMax);
            start_v = round(start_s * 100 / Fs);
            
            num_insert_random = num_insert_random + 1;
            
            dir_out = dir_out_insert_random;
            OLTypeString = STRING_INSERT_RANDOM;
        else
            if start_OL_v <= 0 || (start_pause_v > 0 && num_insert_overlap > num_insert_pause)
                start_v = start_pause_v;
                dir_out = dir_out_insert_pause;
                OLTypeString = STRING_INSERT_INPAUSE;
                num_insert_pause = num_insert_pause + 1;
            else
                start_v = start_OL_v;
                dir_out = dir_out_insert_overlap;
                OLTypeString = STRING_INSERT_INSPEECH;
                num_insert_overlap = num_insert_overlap + 1;
            end
            start_s = round(start_v * Fs / 100);
        end

        end_s = start_s + nd2 - 1;
        end_v = start_v + numel(CF2) - 1;

        data_all = zeros(size(data1));
        
    else  % partial overlap or pause
        
        % ---
        % Find leading/trailing silence and taper it
        % ---

        [data1, TSilLen] = overlaps_taper_trailing_silence(VAD1, data1, Fs);
        [data2, LSilLen] = overlaps_taper_leading_silence(VAD2, data2, Fs);
        if isempty(TSilLen) || isempty(LSilLen)
            nFailedInRow = nFailedInRow + 1;
            continue;
        end
        % offset_v = how much the 2nd was should be shifted from the 1st,
        %   so that there's exactly zero pause/overlap between actual speech
        offset_v = numel(VAD1) - TSilLen - LSilLen;
        
        % ---
        % Generate a random pause/overlap
        % ---
        
        if targetOLType == TYPE_OL_LONG || ...
                (targetOLType == TYPE_ALL && num_OL_big < num_pause && num_OL_big < num_OL_small)
            % make the overlap half the length of 1st file
            pauselen_t = - numel(data1) / (2 * Fs);
            dir_out = dir_out_OL_big;
            num_OL_big = num_OL_big + 1;
            OLType = TYPE_OL_LONG;
            OLTypeString = STRING_OL_LONG;
            
        elseif targetOLType == TYPE_PAUSE || ...
                (targetOLType == TYPE_ALL && num_pause < num_OL_small)
            pauselen_t = randi(1000) / 1000; % random pause of up to 1.0 s
            dir_out = dir_out_pause;
            num_pause = num_pause + 1;
            OLType = TYPE_PAUSE;
            OLTypeString = STRING_PAUSE;
        else % TYPE_OL_SHORT
            pauselen_t = - randi(2000) / 1000; % random overlap of up to 2.0 s
            dir_out = dir_out_OL_small;
            num_OL_small = num_OL_small + 1;
            OLType = TYPE_OL_SHORT;
            OLTypeString = STRING_OL_SHORT;
        end
        
        start_s = round(1 + offset_v * Fs / 100 + pauselen_t * Fs);
        start_v = round(1 + offset_v + pauselen_t * 100);
        
        if start_s < 1
            % first file is probably too short
            nFailedInRow = nFailedInRow + 1;
            continue;
        end

        end_s = start_s + nd2 - 1;
        end_v = start_v + numel(CF2) - 1;
        
        data_all = zeros(nd2 + start_s, 1);
    end
    

    % ---
    % Join wavs
    % ---

    if OLType == TYPE_PAUSE
        r = 0.5 + 0.1 * randn(1);
        if r > xmax, r = xmax; end
        if r < xmin, r = xmin; end
    else
        r = xmin+rand(1)*(xmax-xmin);
    end

    data_all(1:nd1) = data1 * r;
    
    data_all(start_s:end_s) = data_all(start_s:end_s) + data2 * (1-r);
    data_all = data_all(1:max(end_s, nd1));

    CF_both = [];
    CF_both(1:numel(CF1),1) = CF1;
    CF_both(start_v:end_v,2) = CF2;
    
    VAD_both = zeros(size(CF_both));
    VAD_both(1:numel(CF1),1) = VAD1;
    VAD_both(start_v:end_v,2) = VAD2;
    
    
    % ---
    % Limit wav to just 5s on each side of the overlap/pause
    % (except for the large overlaps)
    % ---
    
    if OLType ~= TYPE_OL_LONG
    
        maxSoloSpkLen_s = maxSoloSpkLen * Fs;

        switch OLType
            case TYPE_PAUSE
                boundaryL_v = numel(VAD1) - TSilLen;
                boundaryR_v = start_v + LSilLen;
            case TYPE_INSERT
                boundaryL_v = start_v + LSilLen;
                boundaryR_v = end_v - TSilLen;
            case TYPE_OL_SHORT
                boundaryL_v = start_v + LSilLen;
                boundaryR_v = numel(VAD1) - TSilLen;
            otherwise
                error('Invalid option');
        end

        boundaryL_s = round(boundaryL_v * Fs / 100);
        boundaryR_s = round(boundaryR_v * Fs / 100);
        
        if numel(data_all) - boundaryR_s > maxSoloSpkLen_s
            newEnd_s = boundaryR_s + maxSoloSpkLen_s;
            newEnd_v = boundaryR_v + maxSoloSpkLen * 100 - 1;
            if newEnd_v > end_v
                newEnd_v = end_v;
            end
            data_all = data_all(1:newEnd_s);
            CF_both = CF_both(1:newEnd_v,:); 
            VAD_both = VAD_both(1:newEnd_v,:); 
        end

        if boundaryL_s > maxSoloSpkLen_s
            newStart_s = boundaryL_s - maxSoloSpkLen_s;
            newStart_v = boundaryL_v - maxSoloSpkLen * 100;
            data_all = data_all(newStart_s:end);
            CF_both = CF_both(newStart_v:end,:);
            VAD_both = VAD_both(newStart_v:end,:);
        end

        
    end
        
    
    % ----
    % Add noise and reverb
    % ----
    
    outname_base = [basename1 '-' basename2 '-' OLTypeString];
    strrep(outname_base,'.','_'); % periods will be used as separators ->
        % -> replace existing ones with underscores
    outname = outname_base;

    % tiny level of white noise - just in case
    data_all = data_all + 0.0001 * randn(length(data_all), 1);  
    
    % office noise
    [data_all,noiseString] = overlaps_add_DEMAND_noise(data_all, Fs, dir_DEMAND_noise);
    outname = [outname '.' noiseString]; %#ok<AGROW>
    
    % reverberation
    add_reverb = randi(2) - 1;
    if add_reverb == 1
        [data_all,rirString] = overlaps_add_AIR_reverb(data_all, Fs, dir_AIR_reverb);
        dir_out = [dir_out SUFFIX_REVERB]; %#ok<AGROW>
        outname = [outname '.' rirString]; %#ok<AGROW>
    end
    
    
    
    % ----
    %  Check if this combo of (file1 + file2 + overlap type)
    %   already exists, skip it if it does
    % ----
    
    d = dir([dir_out filesep outname_base '.*']);
    if ~isempty(d)
        nFailedInRow = nFailedInRow + 1;
        continue;
    end
    
    % ---
    % Save it
    % ---
    
    if ~exist(dir_out,'dir')
        mkdir(dir_out)
    end

    % wav is saved with the noise and reverb types as part of the filename,
    %   ref is saved without them
    audiowrite([dir_out '/' outname '.wav'],data_all, Fs);
    save([dir_out '/' outname_base '_VAD.mat'],'VAD_both','CF_both');

        
    %%
    if makePlots
        figure(2018051601);
        ax1 = subplot(2,1,1);
        plot((1:numel(data_all)) / Fs, data_all);
        ax2 = subplot(2,1,2);
        plot((1:size(VAD_both, 1)) / 100, VAD_both);
        hold on;
        plot((1:size(CF_both, 1)) / 100, CF_both);
        linkaxes([ax1 ax2], 'x');
        hold off;
        
        breakhere = true; % good place for a breakpoint
    end

    nOutFiles = nOutFiles + 1;
    nFailedInRow = 0;
end

end

