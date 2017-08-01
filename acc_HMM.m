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

plane1 = var(acc);
s_acc = zeros(size,1);
s_obs = ones(size,1);

parfor n = 1:size
    s_acc(n) = acc(n,:)*plane1.'/norm(plane1,2);
end

%% generate observation labels

num_bins = 4;

s_bins = linspace(min(s_acc), max(s_acc), num_bins+1);

parfor t = 1:size
    for b = 1 : num_bins
        if s_acc(t) > s_bins(b)
            s_obs(t) = b;
        end
    end
end

%% generate default transition matrix and emission matrix

num_stages = 3;
default_home_state = 0.95

trGuess = zeros(num_stages,num_stages);
for j = 1:num_stages
    trGuess(j,j) = default_home_state;
    if j == num_stages
        trGuess(j,1) = 1 - default_home_state;
    else
        trGuess(j,j+1) = 1 - default_home_state;
    end
end

emitGuess = ones(num_stages,num_bins)/num_bins;

% trGuess =[0.95 0.05 0; 0 0.95 0.05; 0.05 0 0.95];
% emitGuess = [0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25];

%% plot obs labels

figure(1);
plot(time, acc); hold on
plot(time, s_acc); hold on
title(['signal at s = ' num2str(num_stages) ' and b =' num2str(num_bins)]);
xlabel('t');
legend('x','y','z','s-acc');
xc = xlim;
xl = xc(1)*0.8+xc(2)*0.2;
yc = ylim;
yl1 = yc(1)*0.83+yc(2)*0.17;
yl2 = yc(1)*0.77+yc(2)*0.23;
% text(xl,yl1,['num stages =' num2str(num_stages)],'FontSize',12);
% text(xl,yl2,['num bins =' num2str(num_bins)],'FontSize',12);
filename = ['signal-s-' num2str(num_stages) '-b-' num2str(num_bins) '.png'];
saveas(gcf, filename,'png');

figure(2);
plot(time, s_obs); hold on
plot(time, s_acc); hold on
legend('s-obs','s-acc');
title(['observation at s = ' num2str(num_stages) ' and b =' num2str(num_bins)]);
xlabel('t');
xc = xlim;
xl = xc(1)*0.8+xc(2)*0.2;
yc = ylim;
yl1 = yc(1)*0.83+yc(2)*0.17;
yl2 = yc(1)*0.77+yc(2)*0.23;
% text(xl,yl1,['num stages =' num2str(num_stages)],'FontSize',12);
% text(xl,yl2,['num bins =' num2str(num_bins)],'FontSize',12);
filename = ['observation-s-' num2str(num_stages) '-b-' num2str(num_bins) '.png'];
saveas(gcf, filename,'png');

%% train HMM

[estTr,estEm] = hmmtrain(s_obs,trGuess,emitGuess);

rng(1);
[seq,states] = hmmgenerate(size,estTr,estEm);

figure(3)
plot(time, s_obs); hold on
plot(time, s_acc); hold on
plot(time, states); hold on
legend('s-obs','s-acc','states');
title(['prediction at s = ' num2str(num_stages) ' and b =' num2str(num_bins)]);
xlabel('t');
xc = xlim;
xl = xc(1)*0.8+xc(2)*0.2;
yc = ylim;
yl1 = yc(1)*0.83+yc(2)*0.17;
yl2 = yc(1)*0.77+yc(2)*0.23;
% text(xl,yl1,['num stages =' num2str(num_stages)],'FontSize',12);
% text(xl,yl2,['num bins =' num2str(num_bins)],'FontSize',12);
filename = ['prediction-s-' num2str(num_stages) '-b-' num2str(num_bins) '.png'];
saveas(gcf, filename,'png');

end

