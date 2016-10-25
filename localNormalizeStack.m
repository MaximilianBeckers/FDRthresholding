clear all;

%get the movie stack
[FileNameMRC,PathNameMRC] = uigetfile('*.mrc' ,'Please select the image file');
[~, filename, ~] = fileparts(FileNameMRC);

%ask for directory where output should be saved
outputPath = uigetdir('','Please select a folder for the ouput');

%read first frame of movie
[map, s] = ReadMRC(fullfile(PathNameMRC, FileNameMRC));
map = double(map);
%local normalize complete stack
for frameInd = 1:s.nz
	
	%do the actual normalization
	map(:, :, frameInd) = localNormalize(map(:, :, frameInd));

end

%now write the locNormed movie as mrc
WriteMRC(map, s.pixA, fullfile(outputPath, strcat(filename, '_localNormal.mrcs'))); 
