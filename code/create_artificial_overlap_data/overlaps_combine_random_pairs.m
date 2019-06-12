function overlaps_combine_random_pairs(rng_seed, input_dir, dir_train, dir_test, offset, speakerDelimiter, train_set_size)
% overlaps_combine_random_pairs(rng_seed, input_dir, dir_train, dir_test, offset, speakerDelimiter, train_set_size)
%
% Create artificial overlaps from LibriSpeech or TIMIT data - Step 2
% Combines random pairs of single-speaker files into artifical overlap data
%
% Should be called from "TIMIT_make_artificial_overlaps" or "Libri_make_artificial_overlaps"
%
% Note: A single call creates at most (nfiles / 2) combinations - 
%   each file is used once, same-speaker combinations are discarded
%
% Args:
%   rng_seed - specify the rng seed (to ensure that multiple calls use different ones)
%   input_dir - input dir with single-speaker wav files and VAD
%   dir_train - output dir for the training set
%   dir_test - output dir for the test set
%   offset - start time (in seconds) for the 2nd speaker
%   speakerDelimiter - specifies the delimiter separating the speaker
%       identifier from the rest of the filename in the input files
%       (assumes that the filenames start with speaker id, can be left empty)
%   train_set_size - how many files should be used as the training set; 
%       the rest will be the test set 
%     can be specified as a target number (train_set_size >= 1) of output training files, 
%       or as a ratio of all files (0 <= train_set_size < 1),
%       e.g. 0.5 for a 50:50 train x test split
%       (to put everything in the train set, use train_set_size = Inf)
%     training set = the first X files in the order they are listed by dir();
%      => should only be used if speaker names and order are meaningless 
%         (=> not suitable for TIMIT, but LibriSpeech is ok)
%
% MK, 2018-04-24
%
% 2019-06-12
%   renamed and cleaned up the function
%   removed unneeded support for CF as an alternative to VAD
%   train_set_size can now also be a ratio, not just a specific number
% 2018-07-19, 07-31
%   turned some hardcoded variables into input args, to make this reusable
%       with other corpora

if ~isempty(rng_seed)
    rng(rng_seed);
end
xmin=0.25;
xmax=0.75;

vps = 100; % VAD values per second

if ~exist(dir_train,'dir')
    mkdir(dir_train);
end
if ~exist(dir_test,'dir')
    mkdir(dir_test);
end

D = dir([input_dir '/*.wav']);
nfiles = numel(D);

if ~exist('train_set_size','var')
    % do a 50:50 split by default
    train_set_size = 0.5;
end

if train_set_size < 0
    error('train_set_size must be >= 0');
end
if train_set_size < 1 % it was given as a ratio
    train_set_size = floor(0.5 * nfiles * train_set_size);
end
train_set_size_input = 2 * train_set_size; % no. of output files -> no. of input files

fnames = {D.name}.';

% shuffle the files, but keep train and test sets separate
if train_set_size_input < nfiles
    fnames1 = fnames(1:train_set_size_input);
    fnames2 = fnames((train_set_size_input+1):end);
    
    fnames1 = fnames1(randperm(length(fnames1)));
    fnames2 = fnames2(randperm(length(fnames2)));

    fnames = [fnames1; fnames2];
else
    fnames = fnames(randperm(length(fnames)));
end


dir_out = dir_train;

for ii = 2:2:nfiles
    
    if ii > train_set_size_input
        dir_out = dir_test;
    end
    
    idx1 = ii;
    idx2 = ii-1;
    
    [~,audioname1,~] = fileparts(fnames{idx1});
    [~,audioname2,~] = fileparts(fnames{idx2});
    
    if ~isempty(speakerDelimiter)
        C = strsplit(audioname1, speakerDelimiter);
        spk1 = C{1};

        C = strsplit(audioname2, speakerDelimiter);
        spk2 = C{1};

        if strcmp(spk1,spk2)
            continue;
        end
    end

    outname = [audioname1 '-' audioname2];
    
    % Check if this combination of input files already exists (in either order), 
    %   skip it if yes
    outname_reversed = [audioname2 '-' audioname1];
    if exist([dir_out '/' outname '.wav'],'file') || ...
            exist([dir_out '/' outname_reversed '.wav'],'file')
        continue;
    end
    
    [data1, FS] = audioread([input_dir '/' fnames{idx1}]);
    [data2, FS2] = audioread([input_dir '/' fnames{idx2}]);
    if FS ~= FS2
        error('Different sample rates: %s vs %s', fnames{idx1}, fnames{idx2});
    end
    
    L = min(length(data1), length(data2));
    
    r=xmin+rand(1)*(xmax-xmin);
    data_out = [r * data1(1:L); zeros(offset * FS,1)] + [zeros(offset * FS,1); (1-r) * data2(1:L)];
    
    % add a tiny level of white noise (to fix zero-filled artificial silences)
    data_out = data_out + 0.0001 * randn(length(data_out), 1);  
    
    audiowrite([dir_out '/' outname '.wav'], data_out, FS);

    % save VAD
    VADLen_noOff = round(L * vps / FS);
    VADLen_total = round((L / FS + offset) * vps);
    

    VAD_both = zeros(VADLen_total, 2);
    VAD = [];
    
    load([input_dir '/' audioname1 '_VAD.mat'], 'VAD');

    if numel(VAD) > VADLen_noOff
        VAD_both(1:VADLen_noOff,1) = VAD(1:VADLen_noOff);
    else
        VAD_both(1:numel(VAD),1) = VAD;
    end
    load([input_dir '/' audioname2 '_VAD.mat'], 'VAD');

    VADLen_off = VADLen_total - VADLen_noOff;
    if numel(VAD) > VADLen_noOff
        VAD_both((VADLen_off+1):end,2) = VAD(1:VADLen_noOff); %#ok<NASGU>
    else
        VAD_both((VADLen_off+1):(VADLen_off+numel(VAD)),2) = VAD; %#ok<NASGU>
    end
    
    save([dir_out '/' outname '_VAD.mat'], 'VAD_both');

end
