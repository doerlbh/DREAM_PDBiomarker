% stage characterization with HMM
% Baihan Lin
% Columbia University
% July 2017

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

test_file = 'test_walk_outbound.tmp';

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

fig1 = figure(1);
plot(time, acc);

data = acc.';
setSeed(0);
maxIt = 30;
nRndRest = 10;

nstates = 3;

% only choose x and y
data = data(1:2,:);
d = 2;

% test with a bogus prior
if 1
    prior.mu = ones(1, d);
    prior.Sigma = 0.1*eye(d);
    prior.k = d;
    prior.dof = prior.k + 1;
% else 
%     prior.mu = [1 3 5 2 9 7 0 0 0 0 0 0 1];
%     prior.Sigma = randpd(d) + eye(d);
%     prior.k = 12;
%     prior.dof = 15;
end

model = hmmFitEm(data, nstates, 'gauss', 'verbose', true, 'piPrior', ones(1,nstates), ...
    'emissionPrior', prior, 'nRandomRestarts', nRndRest, 'maxIter', maxIt);

T = 100;
[observed, hidden] = hmmSample(model, T, 1);
figure; hold on
[styles, colors, symbols, str] =  plotColors();

for k=1:nstates
  gaussPlot2d(model.emission.mu(:,k), model.emission.Sigma(:,:,k),...
    'color',colors(k),'plotMarker','false');
  ndx=(hidden==k);
  plot(observed(1,ndx), observed(2,ndx), sprintf('%s%s', colors(k), symbols(k)));
end

for t=1:T
  ndx=hidden(t);
  text(observed(1,t), observed(2,t), sprintf('%d', t), ...
    'color', colors(ndx), 'fontsize', 14);
end

plot(observed(1,:),observed(2,:),'k-','linewidth',1);

figure; hold on
for k=1:nstates
  ndx=find(hidden==k);
  plot(ndx, hidden(ndx), 'o', 'color', colors(k));
end
axis_pct


