% Function: slice data into sections
% Baihan Lin
% Columbia University
% July 2017

function [sections, mu, sigma, dataNew] = sliceStep(data, slices)

size = length(data);
sections = floor(size/slices);

dataNew = cell(slices,1);
mu = zeros(slices,1);
sigma = zeros(slices,1);

for s = 1 : slices - 1
      dataNew{s} = data(:,((s-1)*sections+1):s*sections);
      mu(s) = mean(dataNew{s},2);
      sigma(s) = std(dataNew{s}.').';
end
dataNew{slices} = data(:,((slices-1)*sections+1):size);
mu(slices) = mean(dataNew{slices},2);
sigma(slices) = std(dataNew{slices}.').';

end
