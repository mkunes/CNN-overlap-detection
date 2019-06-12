function [reverb_data, rirString] = overlaps_add_AIR_reverb(data, Fs, dir_AIR_reverb, rng_seed)
% Reverberate an audio signal using a random room impulse response 
%   from the Aachen Impulse Response (AIR) Database
%
% see also: overlaps_add_DEMAND_noise
%
% MK, 2018-05
%
% 1.1.1 2019-04-20
%   Added the option to specify the RNG seed 
%   (can be used to ensure that two calls use the same RIR)
% 1.1   2018-12-13
%   Corrected how reverb is generated
% 1.0.1 2018-08-02
%   Moved this into its own file
%   Function now also returns a string indicating the room type
% 1.0   2018-05
%   Initial version

    if exist('rng_seed','var')
        rng(rng_seed);
    end

    D = dir([dir_AIR_reverb '/*.mat']);
    
    r = randi(numel(D));
    
    h_air = [];
    load([dir_AIR_reverb filesep D(r).name], 'h_air', 'air_info');

    if size(h_air,2)~=1
        h_air = h_air';
    end

    if air_info.fs ~= Fs
        h_air=resample(h_air,Fs,air_info.fs);
    end
    
    reverb_data = conv(data, h_air);
    reverb_data = reverb_data(1:numel(data));
    scale = mean(abs(data)) ./ mean(abs(reverb_data));
    reverb_data = scale * reverb_data;

    % get the room type (e.g. "office", "hallway") from the filename
    % (Room type is always the 3rd word of the filename,
    %   except for "air_phone_BT_..." files, where it is 4th)
    C = strsplit(D(r).name, '_');
    if strcmp(C{3}, 'BT')
        rirString = C{4};
    else
        rirString = C{3};
    end
    
end
