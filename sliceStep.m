% Function: slice data into sections
% Baihan Lin
% Columbia University
% July 2017

dataNew = sliceStep(data, slices);

sec = floor(length(data)/slices)

for s = 1 : slices - 1
      dataNew{s} = data(:,((s-1)*sec+1):s*sec);
end
dataNew{slices} = data(:,((slices-1)*sec+1):end);

end
