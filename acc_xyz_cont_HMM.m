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

system(['mkdir ' path 'output/' date]);

slices = 20;
% note = 'nstates-4-slices-20';
note = 'demo';

pathTestData = [path 'test_outbound'];
pathOut = [path 'output/' date '/' note];

system(['mkdir ' pathOut]);

[~,list] = system(['find ' pathTestData ' -type f -name "*.tmp"']);

if strcmp(note, 'demo')
    
    test_file = 'test_walk_outbound.tmp';
    
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
    
    ncases = 1;
    n = 1;
    figN = figure(n);
    plot(time, acc);
    title('signal')
    xlabel(t);
    legend('x','y','z');
    filename = [pathOut  '/sig-demo-slices-' num2str(slices)];
    saveas(gcf, [filename '.png'],'png');
    saveas(gcf, [filename '.fig']);
    close(figN);
    
    % data = acc.';
    data = acc;
    
    [sections, data] = sliceStep(data, slices);
    steps = slices;
    
else
    
    files = strsplit(list);
    ncases = length(files)-1;
    data = [];
    minsec = 10000000;
    for n = 1:ncases
        file = files{n};
        
        rawData = loadjson(file,'SimplifyCell',1);
        size = length(rawData);
        time = zeros(1,size);
        acc = zeros(3,size);
        
        for t = 1:size
            time(1,t) = rawData(1,t).timestamp;
            acc(1,t) = rawData(1,t).x;
            acc(2,t) = rawData(1,t).y;
            acc(3,t) = rawData(1,t).z;
        end
        
        figN = figure(n);
        plot(time, acc);
        title('signal')
        xlabel(t);
        legend('x','y','z');
        filename = [pathOut  '/sig-case-' num2str(n) '-slices-' num2str(slices)];
        saveas(gcf, [filename '.png'],'png');
        saveas(gcf, [filename '.fig']);
        close(figN);
        
        % data = acc.';
        dataUncomb = acc;
        [secN, dataUncomb] = sliceStep(dataUncomb, slices);
        if secN < minsec
            minsec = secN;
        end
        data = [data; dataUncomb];
    end
    sections = minsec;
    steps = length(data);
    
end

%% train HMM

setSeed(0);
maxIt = 30;
nRndRest = 10;

nstates = 4;
d = 3;

% % only choose x and y
% data = data(1:2,:);
% d = 2;

% test with a bogus prior
if 1
    prior.mu = ones(1, d);
    prior.Sigma = 0.1*eye(d);
    prior.k = d;
    prior.dof = prior.k + 1;
end

model = hmmFitEm(data, nstates, 'gauss', 'verbose', true, 'piPrior', ones(1,nstates), ...
    'emissionPrior', prior, 'nRandomRestarts', nRndRest, 'maxIter', maxIt);

T = sections;
stptime = zeros(1,T);
timeObs = zeros(1,T);

[observed, hidden] = hmmSample(model, T, 1);
% figure;


%% plot clustering

fig2 = figure(ncases+1); hold on
title(['stage characterizations of walking outbound acc']);
[styles, colors, symbols, str] =  plotColors();

% for t = 1 : T - 1
%    ndx=hidden(t);
%    rate = model.pi(ndx);
%    tau = log(1/rand()) / rate;
%    stptime(t+1) = tau;
%    timeObs(t+1) = timeObs(t) + tau;
% end

for k=1:nstates
    gaussPlot2d(model.emission.mu(1:2,k), model.emission.Sigma(1:2,1:2,k),...
        'color',colors(k),'plotMarker','false');
    xlabel('x');
    ylabel('y');
    %   ndx=(hidden==k);
    %   plot(observed(1,ndx), observed(2,ndx), sprintf('%s%s', colors(k), symbols(k)));
end


for t=1:T-1
    ndx=hidden(t);
    text(observed(1,t), observed(2,t), sprintf('%d', t), ...
        'color', colors(ndx), 'fontsize', 14);
    
    plot(observed(1,t:t+1),observed(2,t:t+1),'k-','linewidth',1);
    if strcmp(note, 'demo')
        pause(0.1);
    end
end

filename = [pathOut '/stages-slices-' num2str(slices)];
saveas(gcf, [filename '.png'],'png');
saveas(gcf, [filename '.fig']);
close(fig2);

% plot(observed(1,:),observed(2,:),'k-','linewidth',1);

%% plot hidden states

% figure;
fig3 = figure(ncases+2); hold on
title(['hidden states']);
xlabel('t-sample');
ylabel('hidden states');
for k=1:nstates
    ndx=find(hidden==k);
    plot(ndx, hidden(ndx), 'o', 'color', colors(k));
end
axis_pct

filename = [pathOut  '/hidden-slices-' num2str(slices)];
saveas(gcf, [filename '.png'],'png');
saveas(gcf, [filename '.fig']);
close(fig3);


%% plot x, y, z observed

% figure;
fig4 = figure(ncases+3);
title(['sampled observations']);
xlabel('t');
ylabel('observations');
hold on
plot(1:T, observed(1,:));
plot(1:T, observed(2,:));
plot(1:T, observed(3,:));
legend('x','y','z');

filename = [pathOut  '/obs-slices-' num2str(slices)];
saveas(gcf, [filename '.png'],'png');
saveas(gcf, [filename '.fig']);
close(fig4);

%% Comparing observed with original datasets

% figure;
% hold on

for s = 1: steps
    
    figS = figure(ncases+3+s);
    title(['observed vs original datasets']);
    
    data_sec = data{s};
    
    subplot(3,1,1);
    plot(time(1:T), observed(1,1:T)); hold on
    plot(time(1:T), data_sec(1,1:T));
    legend('obs_x','data_x');
    xlabel('t');
    
    subplot(3,1,2);
    plot(time(1:T), observed(2,1:T)); hold on
    plot(time(1:T), data_sec(2,1:T));
    legend('obs_y','data_y');
    xlabel('t');
    
    subplot(3,1,3);
    plot(time(1:T), observed(3,1:T)); hold on
    plot(time(1:T), data_sec(3,1:T));
    legend('obs_z','data_z');
    xlabel('t');
    
    filename = [pathOut '/comp-slices-' num2str(slices) '-s' num2str(s)];
    saveas(gcf, [filename '.png'],'png');
    saveas(gcf, [filename '.fig']);
    close(figS);
    
end



