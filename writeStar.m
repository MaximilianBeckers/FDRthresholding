

%*************************************************************************
%********************* print .star file for relion ***********************
%*************************************************************************

%print header
fileID = fopen('images.star','w');
fprintf(fileID,'data_ \n');
fprintf(fileID,'_loop \n');
fprintf(fileID,'_rlnImageName \n');
fprintf(fileID,'_rlnMicrographName \n');
fprintf(fileID,'_rlnDefocusU \n');
fprintf(fileID,'_rlnDefocusV \n');
fprintf(fileID,'_rlnDefocusAngle \n');
fprintf(fileID,'_rlnVoltage \n');
fprintf(fileID,'_rlnSphericalAberration \n');
fprintf(fileID,'_rlnAmplitudeContrast \n');

%print data
for particleInd = 1:numParticles
	fileName = fullfile(PathName{particleInd}, FileName{particleInd});
	fprintf(fileID, '%6u@%s %s %s %f %f %f %f \n', particleInd, fileName, fileName, def1, def2, defang, volt, sphAbb, AmplCont);
end

fclose(fileID);
