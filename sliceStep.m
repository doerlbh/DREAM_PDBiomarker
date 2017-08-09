% Function: slice data into sections
% Baihan Lin
% Columbia University
% July 2017

function [sections, dataNew] = sliceStep(data, slices);

size = length(data);
sections = floor(size/slices);

dataNew = cell(slices,1);

for s = 1 : slices - 1
      dataNew{s} = data(:,((s-1)*sections+1):s*sections);
end
dataNew{slices} = data(:,((slices-1)*sections+1):size);

end
