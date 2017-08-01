% stage characterization with HMM
% Baihan Lin
% Columbia University
% July 2017

clear all; close all;


% path = '/gscratch/stf/sunnylin/Columbia/DREAM_PDBiomarker/';
% path = '/home/sunnylin/Dropbox/Git/DREAM_PDBiomarker/';
path = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/';

pathN = path;

% loc = 'Hyak';
% loc = 'Hyakqsub';
% loc = 'americano';
% loc ='galao';
% loc = 'latte';
% loc = 'espresso';
% loc = 'mocha';
loc = 'doerlbh';

test_file = 'test_walk_outbound.tmp';

addpath(path);
addpath([path 'jsonlab-master']);
% addpath([path 'parameters']);
% addpath([path 'embeddings']);

data = loadjson(test_file,'SimplifyCell',1);
sum = length(data);
time = zeros(1,sum);
acc = zeros(3,sum);

for t = 1:sum
    time(1,t) = data(1,t).timestamp;
    acc(1,t) = data(1,t).x;
    acc(2,t) = data(1,t).y;
    acc(3,t) = data(1,t).z;
end



