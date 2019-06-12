MATLAB scripts for generating artificial overlap data

	Libri_make_artificial_overlaps_single_overlap.m 
		- creates several types of short overlaps from LibriSpeech data

	Libri_make_artificial_overlaps.m
		- creates long overlap data from LibriSpeech

	TIMIT_make_artificial_overlaps.m
		- creates long overlap data from TIMIT

Currently, these scripts rely on in-house VAD, so they're not usable out of the box. (Replacing the VAD is on my TODO list.)

Plus there are still some semi-hardcoded paths, but those should be easier to replace.
