function Libri_make_artificial_overlaps
% Create long artificial overlap data from the LibriSpeech test set
%
% MK, 2019-04-18
%
% Based on "TIMIT_make_artificial_overlaps"
% Also saves multiple versions of the overlap data - different combinations
%   of noise/reverb
%
% Changelog:
%   2019-06-12
%       General cleanup, functionality should be unchanged
%   2019-04-18
%       Initial version (based on "TIMIT_make_artificial_overlaps")
%
% TODO:
%   - Make this compatible with some kind of publicly available VAD, 
%       not just the one in Param_KW
%   - Move all paths to a config file

mypath_datasets = getSystemSpecificPath('datasets');
wav_parentdir = [mypath_datasets '/LibriSpeech/test-other'];

dir_joined = [mypath_datasets '/LibriSpeech/test-other-joined'];
dir_out = [mypath_datasets '/LibriSpeech/test-other-overlaps'];

dir_prm = [mypath_datasets '/LibriSpeech/test-other-prm'];
dir_clean = [dir_out '_clean'];

dir_DEMAND_noise = [mypath_datasets '/Demand_DB_noise'];
dir_AIR_reverb = [mypath_datasets '/AIR_1_4/for_AMI'];


% ------------
% Step 1: Join files from the same speaker, with 5-10s pauses in-between
%  (linearly taper silence at the edges, to avoid seams)
% ------------

% Expects the original LibriSpeech folder structure:
%   mainDir/speakerID/bookID/*.flac

minPause = 5;
maxPause = 10;
maxTotalTime = 600; % make at most ~10 min files

VAD_lenTh = 25; 
VAD_minPauseLen = 50;

mkdir(dir_joined);

spkDirs = dir(wav_parentdir);
nspk = numel(spkDirs);

for ii = 1:nspk
    if ~spkDirs(ii).isdir || spkDirs(ii).name(1) == '.'
        continue;
    end
    
    bkDirs = dir([wav_parentdir '/' spkDirs(ii).name]);
    nbk = numel(bkDirs);
    
    
    for jj = 1:nbk
        if ~bkDirs(jj).isdir || bkDirs(jj).name(1) == '.'
            continue;
        end
        
        filepath = [wav_parentdir '/' spkDirs(ii).name '/' bkDirs(jj).name '/'];
        outname = [spkDirs(ii).name '-' bkDirs(jj).name];

        flacFiles = dir([wav_parentdir '/' spkDirs(ii).name '/' bkDirs(jj).name '/*.flac']);
        nfiles = numel(flacFiles);
        
        data_all = [];
        CF_all = [];
        
        FS_prev = 0;
        
        endTime = 0;
        for k = 1:nfiles
            [data,FS] = audioread([filepath flacFiles(k).name]);
            
            [~,partname,~] = fileparts(flacFiles(k).name);
            
            [~,CF] = ReadPRM([dir_prm '/' partname '.prm']);
            

            if k > 1 && FS_prev ~= FS
                error('Inconsistent sample rate in %s',filepath);
            end
            FS_prev = FS;
            
            if k > 1
                silLen = minPause + (maxPause-minPause) * rand(1);
            else
                silLen = 0;
            end
            
            startTime = endTime + silLen;
            endTime = startTime + numel(data) / FS; 
            
            data_all = [data_all; zeros(round(silLen * FS),1); data]; %#ok<AGROW>
            CF_all = [CF_all; zeros(round(silLen * 100),1); CF]; %#ok<AGROW>

            if endTime > maxTotalTime
                break;
            end
        end
        
        audiowrite([dir_joined '/' outname '.wav'],data_all, FS);
        
        CF = fix_CF(CF_all,data_all, FS);
        VAD = get_VAD_from_CF(CF, VAD_minPauseLen, VAD_lenTh); %#ok<NASGU>
        save([dir_joined '/' outname '_VAD.mat'],'CF','VAD');
    end
    
end


%%

% ------------
% Step 2: Combine pairs of speakers
% ------------

% each call uses every file from step 1 just once -> max nspeakers/2
%   combinations per call
% to create more files, run this multiple times with different rng seeds

for rng_seed = 1230:1235
    overlaps_combine_random_pairs(rng_seed, dir_joined, dir_clean, dir_clean, 2, '-', 0.75);
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
    

    % clean + office noise
    [data_noise, noiseString] = overlaps_add_DEMAND_noise(data, Fs, dir_DEMAND_noise);
    audiowrite([dir_out_noise '/' basename '.' noiseString '.wav'],data_noise, Fs);
    
    rng_seed = randi(1000000); % use the same random seed for both calls of overlaps_add_AIR_reverb
    
    % noise + reverberation
    [data_noise, rirString] = overlaps_add_AIR_reverb(data_noise, Fs, dir_AIR_reverb, rng_seed);
    audiowrite([dir_out_noise_reverb '/' basename '.' noiseString '.' rirString '.wav'],data_noise, Fs);
    clear data_noise

    % clean + reverberation
    [data, rirString] = overlaps_add_AIR_reverb(data, Fs, dir_AIR_reverb, rng_seed);
    audiowrite([dir_out_reverb '/' basename '.' rirString '.wav'],data, Fs);

end

