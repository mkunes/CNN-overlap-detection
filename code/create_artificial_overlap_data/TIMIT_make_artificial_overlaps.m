function TIMIT_make_artificial_overlaps(wav_dir, dir_joined, dir_out, skipFirstStep)
% Create artificial overlap data from (a subset of) the TIMIT train set
%
% MK, 2018-07-19
%
% Changelog:
%   2019-06-12
%       General cleanup
%       Works with the original TIMIT folder structure now
%           (.../TRAIN/dialectRegion/Speaker/sentence.wav)
%   2019-05-03
%       Changed the name format - noise and reverb type are now delimited
%           by periods instead of underscores
%       Made the first step optional (joining the wavs from the same
%           speaker together)
%       Turned some of the hardcoded paths into optional input params
%   2018-08-02
%       Same-speaker utterances are now randomly shuffled
%       Added noise and reverb
%   2018-07-19
%       Original version
%
% TODO:
%   - Get VAD directly from the original TIMIT transcript files,
%       instead of just loading it from a .mat file
%       (technically, all that's needed is the index of the first and last 
%        speech frame, since this script ignores pauses in the middle anyway)
%   - Move all paths to a config file
%   

% ------------
% Paths
% ------------

mypath_datasets = getSystemSpecificPath('datasets');

if nargin == 0
    skipFirstStep = true;
    dir_out = [mypath_datasets '/TIMIT/TRAIN_overlaps/test_v2'];
    dir_joined = [mypath_datasets '/TIMIT/TRAIN_joined_randomised/test'];
    wav_dir = [mypath_datasets '/TIMIT/TRAIN'];
end

dir_DEMAND_noise = [mypath_datasets '/Demand_DB_noise'];
dir_AIR_reverb = [mypath_datasets '/AIR_1_4/for_AMI'];

% "TIMIT_VAD_TRAIN.mat" contains VAD (0/1) for a *subset* of the TRAIN wavs
    %   (some wavs have no VAD, but all VAD entries have a corresponding wav)
    % It's also already alphabetically sorted (by speaker)
file_vad = [mypath_datasets '/TIMIT/TIMIT_VAD_TRAIN.mat'];

% ------------
% Step 1: Join files from the same speaker, with 0-2s pauses in-between
%  (linearly taper silence at the edges, to avoid seams)
% ------------

if ~skipFirstStep

    mkdir(dir_joined);

    DATA = [];
    load(file_vad, 'DATA');

    nfiles = numel(DATA);

    spk = '';
    data_all = [];
    VAD = [];
    wavIds = [];
    
    vps = 100; % VAD values per second

    for ii = 1:nfiles+1

        if ii <= nfiles
            C = strsplit(DATA(ii).name, '_');
            spk_new = C{1};
        end

        % new speaker -> save the previous one
        if ii > nfiles || ~strcmp(spk, spk_new) 
            if ~isempty(wavIds)

                % shuffle the wav order (because some texts are the same)
                wavIds = wavIds(randperm(length(wavIds)));

                for iwav = 1:numel(wavIds)
                    
                    p = [wav_dir '/*/' strrep(DATA(wavIds(iwav)).name,'_','/') '.wav'];
                    d = dir(p);
                    if isempty(d)
                        warning('Found no files matching ''%s''. Skipping.',p);
                        continue;
                    elseif numel(d) > 1
                        warning('Found multiple files matching ''%s''. Using the first one found.',p);
                    end
                    
                    [data,Fs] = audioread([d(1).folder '/' d(1).name]);
                    vad_new = DATA(wavIds(iwav)).vad;

                    % let's assume there are no significant pauses inside the utterances
                    %  -> speech starts at the first *speech* frame and ends at the last
                    % (these are 2-5s long single-sentence utterances, seems reasonable)

                    % taper leading and trailing silence, label everything else as speech
                    [data, TSilLen] = overlaps_taper_trailing_silence(vad_new, data, Fs);
                    [data, LSilLen] = overlaps_taper_leading_silence(vad_new, data, Fs);
                    vad_new(LSilLen:(end-TSilLen)) = true; 

                    if isempty(data_all)
                        pauselen_t = 0;
                    else
                        pauselen_t = randi(2 * vps) / vps; % random pause of up to 2.0 s
                            % (same precision as VAD - to avoid rounding errors)
                    end

                    vad_pause = false(round(pauselen_t * vps), 1);
                    data_pause = zeros(round(pauselen_t * Fs), 1);

                    VAD = [VAD; vad_pause; vad_new];
                    data_all = [data_all; data_pause; data];

                end
                
                if ~isempty(data_all) > 0
                    audiowrite([dir_joined '/' spk '.wav'],data_all, Fs);
                    save([dir_joined '/' spk '_VAD.mat'], 'VAD'); 
                end

                wavIds = [];
            end
            spk = spk_new;
            data_all = [];
            VAD = [];


        end

        wavIds(end+1) = ii;


    end

end    
    
%%

% ------------
% Step 2: Combine pairs of speakers
% ------------

% each call uses every file from step 1 just once -> max nspeakers/2 combinations per cal
% to create more files, run this multiple times with different rng seeds

dir_clean = [dir_out '_clean'];

for rng_seed = 1230:1235
    overlaps_combine_random_pairs(rng_seed, dir_joined, dir_clean, dir_clean, 2, '_', Inf);
end


%%

% ------------
% Step 3: Add more background noise 
%  (slight white noise is already added in step 2)
% ------------


D = dir([dir_clean '/*.wav']);

dir_out_noise = [dir_out '_noise'];
dir_out_noise_reverb = [dir_out '_noise_reverb'];
dir_out_reverb = [dir_out '_reverb'];

mkdir(dir_out_noise);
mkdir(dir_out_noise_reverb);
mkdir(dir_out_reverb);


for ii = 1:numel(D)
    [data,Fs] = audioread([dir_clean '/' D(ii).name]);
    
    [~, basename, ~] = fileparts(D(ii).name);

    % add office noise
    [data, noiseString] = overlaps_add_DEMAND_noise(data, Fs, dir_DEMAND_noise);
    dir_out_noise = [dir_out '_noise'];
    outName = [basename '.' noiseString];
    
    % add reverberation (50% of files)
    add_reverb = randi(2);
    if add_reverb == 1
        [data, rirString] = overlaps_add_AIR_reverb(data, Fs, dir_AIR_reverb);
        
        dir_out_noise = [dir_out '_noise_reverb'];
        outName = [basename '.' noiseString '.' rirString];
    end

    audiowrite([dir_out_noise '/' outName '.wav'],data, Fs);
    
end


