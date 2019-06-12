function [noisy_data, noiseString] = overlaps_add_DEMAND_noise(data, Fs, dir_DEMAND_noise)
% Add noise from the DEMAND database to an audio signal
% (Concatenate multiple random files of the same noise type, with a random
%  start in the 1st file (others are used whole) and linear transition
%  between files)
%
% see also: overlaps_add_AIR_reverb
%
% MK, 2018-05
%
% v1.01 2018-08-02
%   Moved this into its own file
%   Function now also returns a string indicating the type of added noise
% v1.0  2018-05
%   Initial version

    noiseString = 'noise'; %#ok<NASGU>

    datalen = numel(data);
    
    noiselen = 0;
    noise = [];
    
    r = randi(3);
    switch r
        case 1
            noiseDir = [dir_DEMAND_noise '/OHALLWAY'];
            volume_coeff = 0.5 + 1.5 * rand(1); % change the volume of the noise 0.5-2x
            noiseString = 'OHALLWAY';
        case 2
            noiseDir = [dir_DEMAND_noise '/OMEETING'];
            volume_coeff = 0.1 + 0.4 * rand(1); % change the volume of the noise 0.1-0.5x
            noiseString = 'OMEETING';
        case 3
            noiseDir = [dir_DEMAND_noise '/OOFFICE'];
            volume_coeff = 0.5 + 1.5 * rand(1); % change the volume of the noise 0.5-2x
            noiseString = 'OOFFICE';
        otherwise
            error('Invalid option');
    end
    D = dir([noiseDir '/*.wav']);
    
    nNoiseFiles = 0;
    
    while noiselen < datalen
        r = randi(numel(D));

        [noise_new, Fs2] = audioread([noiseDir '/' D(r).name]);
        
        if Fs2 ~= Fs % sample rate of the noise is different
            % resample the noise to the correct rate
            resample(noise_new,Fs,Fs2);
        end
        
        nNoiseFiles = nNoiseFiles + 1;
        
        if nNoiseFiles == 1
            % first noise file - start at a random point,
            %  no later than 1 sec before end, if possible
            randMax = max(1, numel(noise_new) - Fs);
            start = randi(randMax);
            noise = noise_new(start:end);
        else
            % later noise files - use them whole, 
            %   with a 1 sec linear transition between files
            % (or shorter if the files are very short)
            
            taperLen = min([Fs, noiselen, numel(noise_new)]);

            coeff1 = (taperLen:-1:1) ./ (taperLen + 1);
            coeff2 = (1:taperLen) ./ (taperLen + 1);

            noise((end-taperLen+1):end) = noise((end-taperLen+1):end) .* coeff1' + ...
                noise_new(1:taperLen) .* coeff2';

            noise = [noise; noise_new(taperLen+1:end)]; %#ok<AGROW>

        end
        
        noiselen = numel(noise);
    end

    noisy_data = data + volume_coeff * noise(1:datalen);
    
end
