clear all;

%get binary mask of same size
[FileNameMask,PathNameMask] = uigetfile('*.mrc','Please select the mask file');
[mask, ~] = ReadMRC(fullfile(PathNameMask, FileNameMask));
mask = double(mask);

%get the shifts
[FileNameShift, PathNameShift] = uigetfile('*.txt','Please select the shift file');

%read shifts and get inverse operation	
fileID = fopen(fullfile(PathNameShift, FileNameShift),'r');
shift = fscanf(fileID, '%f %f', [2,inf] );
shift = shift';
shift = shift;
fclose(fileID);

%get unaligned movie
[FileNameMRC,PathNameMRC] = uigetfile('*.mrc', 'Please select the unaligned movie file');
[~, filename, ~] = fileparts(FileNameMRC);


%ask for directory where output should be saved
outputPath = uigetdir('','Please select a folder for the ouput');

[map, s] = ReadMRC(fullfile(PathNameMRC, FileNameMRC));
map = double(map);

%apply mask for each frame
for frameInd = 1:s.nz
	
	%get frame of movie
	mapIm = map(:, :, frameInd);
	
	%do shift of mask
	maskShift = imtranslate(mask, shift(frameInd, :), 0);

	%apply shifted mask
	mapIm = mapIm .* maskShift;

	%update the map with the masked image
	map(:, :, frameInd) = mapIm;
	
end

%now write the masked movie as mrc
WriteMRC(map, s.pixA, fullfile(outputPath, strcat(filename, '_maskedShifted.mrcs'))); 
