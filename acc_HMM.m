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
size = length(data);
time = zeros(1,size);
acc = zeros(3,size);

for t = 1:size
    time(1,t) = data(1,t).timestamp;
    acc(1,t) = data(1,t).x;
    acc(2,t) = data(1,t).y;
    acc(3,t) = data(1,t).z;
end

acc = acc.';
time = time.';

plane1 = var(acc);
s_acc = zeros(size,1);
s_obs = ones(size,1);

parfor n = 1:size
    s_acc(n) = acc(n,:)*plane1.'/norm(plane1,2);
end

num_bins = 4;
s_bins = linspace(min(s_acc), max(s_acc), num_bins+1);

parfor t = 1:size
    for b = 1 : num_bins
        if s_acc(t) > s_bins(b)
            s_obs(t) = b;
        end
    end
end


%% plot obs labels

figure(1);
plot(time, acc);
plot(time, s_acc);

figure(2);
plot(time, s_obs); hold on
plot(time, s_acc); hold on


%% train HMM

trGuess =[0.95 0.05 0; 0 0.95 0.05; 0.05 0 0.95]; 
emitGuess = [0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25]; 
[estTr,estEm] = hmmtrain(s_obs,trGuess,emitGuess);

rng(1);
[seq,states] = hmmgenerate(size,estTr,estEm);

figure(3)
plot(time, s_obs); hold on
plot(time, s_acc); hold on
plot(time, states);

% plot();
