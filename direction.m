% find the direction of the accel
% Baihan Lin
% Columbia University
% October 2017

% function direction(thisFolder, timeFrame, threshold, output_flipped_files)

test_file = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/test_walk_outbound.json';
test_folder = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/test_outbound';
test_output_files = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/test_output.txt';
test_timeFrame = 20;
test_threshold = 20;
thisFolder = test_folder;
output_flipped_files = test_output_files;
timeFrame = test_timeFrame;
threshold = test_threshold;

system(['echo >  ' output_flipped_files]);
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
            %         side = 1;
            %         sum(acc(2,t:t+timeFrame-1))
            %         t
            system(['echo ' fullFileName ' >> ' output_flipped_files]);
            break;
        else if sum(acc(2,t:t+timeFrame-1)) < -threshold
                system(['echo ' fullFileName ' >> ' output_flipped_files '.unflipped.txt']);
                %             sum(acc(2,t:t+timeFrame-1))
                %             t
                break;
                %         else
                %             disp('bad threshold');
            end
        end
    end
     
end




