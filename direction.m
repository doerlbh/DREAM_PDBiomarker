% find the direction of the accel
% Baihan Lin
% Columbia University
% October 2017

% function direction(thisFolder, timeFrame, threshold, output_flipped_files)

test_file = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/test_walk_outbound.json';
test_folder = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/test_outbound';
test_output_files = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/test_output2.txt';
test_timeFrame = 20;
test_threshold = 20;
thisFolder = test_folder;
output_flipped_files = test_output_files;
timeFrame = test_timeFrame;
threshold = test_threshold;

% system(['echo >  ' output_flipped_files]);
fprintf('Processing folder %s\n', thisFolder);
filePattern = sprintf('%s/*.tmp', thisFolder);
baseFileNames = dir(filePattern);
numberOfFiles = length(baseFileNames);

for f = 1 : numberOfFiles
    fullFileName = fullfile(thisFolder, baseFileNames(f).name);
    fprintf('Processing file %s\n', fullFileName);
    
    rawData = loadjson(fullFileName,'SimplifyCell',1);
    size = length(rawData);
    time = zeros(1,size);
    acc = zeros(3,size);
    
    for t = 1:size
        time(1,t) = rawData(1,t).timestamp;
        acc(1,t) = rawData(1,t).x;
        acc(2,t) = rawData(1,t).y;
        acc(3,t) = rawData(1,t).z;
    end
    
    for t = 1:size-timeFrame
        if sum(acc(2,t:t+timeFrame-1)) > threshold
            system(['echo ' fullFileName ' >> ' output_flipped_files]);
            
    for t = 1:size
        field1 = 'y';  value1 = -acc(2,t);
            field2 = 'timestamp';  value2 = time(1,t);
            field3 = 'z';  value3 = acc(3,t);
            field4 = 'x';  value4 = acc(1,t);
            
            flipped_Data(t) = struct(field1,value1,field2,value2,field3,value3,field4,value4)
           
    end
    
    break;
        else if sum(acc(2,t:t+timeFrame-1)) < -threshold
                system(['echo ' fullFileName ' >> ' output_flipped_files '.unflipped.txt']);
                break;
            end
        end
    end
    
end




