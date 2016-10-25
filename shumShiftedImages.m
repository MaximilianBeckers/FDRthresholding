
%get the shifts
[FileNameShift, PathNameShift] = uigetfile('*.txt','Please select the shift file');

%read shifts	
fileID = fopen(fullfile(PathNameShift, FileNameShift),'r');
shift = fscanf(fileID, '%f %f', [2,inf] );
shift = shift';
shift = -1*shift;
fclose(fileID);

%get the movie stack
[FileNameMRC,PathNameMRC] = uigetfile('Please select the image file');
[~, filename, ~] = fileparts(FileNameMRC);

%ask for directory where output should be saved
outputPath = uigetdir('','Please select a folder for the ouput');

%read first frame of movie
[map, s] = ReadMRC(fullfile(PathNameMRC, FileNameMRC), 1, 1);
map = double(map);

%define summed output Image
summedIm = zeros(s.nx, s.ny);

%normalize map
map = localNormalize(map);

%do first shift of image
map = imtranslate(map, shift(1, :), 0);

%add new slice to the summed Image
summedIm = summedIm + map;	

for frameInd = 2:s.nz

	%read map
	[map, ~] = ReadMRC(fullfile(PathNameMRC, FileNameMRC), frameInd, 1);
	map = double(map);

	%normalize map
	map = localNormalize(map);

	%do shift
	map = imtranslate(map, shift(frameInd, :), 0);

	%add new slice to the summed Image
	summedIm = summedIm + map;			
	
end

pixA = 1;
%now write the summed movie as mrc
WriteMRC(summedIm, pixA, fullfile(outputPath, strcat(filename, '_CorSum.mrc'))); 
