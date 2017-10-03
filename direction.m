% find the direction of the accel
% Baihan Lin
% Columbia University
% October 2017

clear all; close all;

% path = '/gscratch/stf/sunnylin/Columbia/DREAM_PDBiomarker/';
% path = '/home/sunnylin/Dropbox/Git/DREAM_PDBiomarker/';
path = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/';

% pathN = path;

% loc = 'Hyak';
% loc = 'Hyakqsub';
% loc = 'americano';
% loc ='galao';
% loc = 'latte';
% loc = 'espresso';
% loc = 'mocha';
loc = 'doerlbh';

addpath(path);
addpath([path 'jsonlab-master']);
% addpath([path 'parameters']);
% addpath([path 'embeddings']);

test_file = '/Users/DoerLBH/Dropbox (Personal)/git/DREAM_PDBiomarker/test_walk_outbound.json';

rawData = loadjson(test_file,'SimplifyCell',1);
size = length(rawData);
time = zeros(1,size);
acc = zeros(3,size);

for t = 1:size
    time(1,t) = rawData(1,t).timestamp;
    acc(1,t) = rawData(1,t).x;
    acc(2,t) = rawData(1,t).y;
    acc(3,t) = rawData(1,t).z;
end



% acc = acc.';
% time = time.';




