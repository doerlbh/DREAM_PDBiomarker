% stage characterization with HMM
% Baihan Lin
% Columbia University
% July 2017

function acc_xyz_cont_HMM_resize(slices,)

clear all; close all;

% path = '/gscratch/stf/sunnylin/Columbia/DREAM_PDBiomarker/';
% path = '/home/sunnylin/Dropbox/Git/DREAM_PDBiomarker/';
path = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/';

addpath(path);
addpath([path 'jsonlab-master']);
% addpath([path 'parameters']);
% addpath([path 'embeddings']);

system(['mkdir ' path 'output/' date]);

slices = 1;
note = 'nstates-10-slices-1';
% note = 'demo';

pathTestData = [path 'test_outbound'];
pathOut = [path 'output/' date '/' note];

system(['mkdir ' pathOut]);
system(['mkdir ' pathOut '/sig']);
system(['mkdir ' pathOut '/hidden']);
system(['mkdir ' pathOut '/obs']);
system(['mkdir ' pathOut '/stages']);
system(['mkdir ' pathOut '/comp']);

[~,list] = system(['find ' pathTestData ' -type f -name "*.tmp"']);


resize_length = 2000;
sample_timerange = [71010 71024;
    2795187 2795230;
    225408 225420;
    354178 354191;
    197253 197264;
    379882 379898;
    2244631 2244644;
    46694 46708;
    19299 19318;
    361677 361692;
    66982 66998;
    44843 44857;
    122531 122546;
    239818 239834;
    166240 166253;
    112776 112787;
    406465 406478;
    525549 525561;
    731633 731647];

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
        %         time = zeros(1,size);
        %         acc = zeros(3,size);
        
        isRecording = 0;
        newSize = 0;
        
        for t = 1:size
            
            if ~isRecording && rawData(1,t).timestamp > sample_timerange(n,1)
                isRecording = 1;
                newSize = newSize + 1;
                time(1,newSize) = rawData(1,t).timestamp;
                acc(1,newSize) = rawData(1,t).x;
                acc(2,newSize) = rawData(1,t).y;
                acc(3,newSize) = rawData(1,t).z;
            else if rawData(1,t).timestamp > sample_timerange(n,2)
                    isRecording = 0;
                else if isRecording
                        newSize = newSize + 1;
                        time(1,newSize) = rawData(1,t).timestamp;
                        acc(1,newSize) = rawData(1,t).x;
                        acc(2,newSize) = rawData(1,t).y;
                        acc(3,newSize) = rawData(1,t).z;
                    end
                end
            end    
            
        end
        
        time = imresize(time, [1 resize_length]);
        acc = imresize(acc, [3 resize_length]);
        
        figN = figure(n);
        plot(time, acc);
        title('signal')
        xlabel(t);
        legend('x','y','z');
        filename = [pathOut  '/sig/sig-case-' num2str(n) '-slices-' num2str(slices)];
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

nstates = 10;
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

model = hmmFitEm(data, nstates, 'gauss', 'verbose', true, 'piPrior', [100 ones(1,nstates-1)], ...
    'emissionPrior', prior, 'nRandomRestarts', nRndRest, 'maxIter', maxIt);

T = sections;
% stptime = zeros(1,T);
% timeObs = zeros(1,T);

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
    pause(0.01);
    if strcmp(note, 'demo')
        %         pause(0.1);
    end
end

filename = [pathOut '/stages/stages-slices-' num2str(slices)];
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

filename = [pathOut  '/hidden/hidden-slices-' num2str(slices)];
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

filename = [pathOut  '/obs/obs-slices-' num2str(slices)];
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
    plot(time(((s-1)*T+1):s*T), observed(1,((s-1)*T+1):s*T)); hold on
    plot(time(((s-1)*T+1):s*T), data_sec(1,((s-1)*T+1):s*T));
    legend('obs_x','data_x');
    xlabel('t');
    
    subplot(3,1,2);
    plot(time(((s-1)*T+1):s*T), observed(2,((s-1)*T+1):s*T)); hold on
    plot(time(((s-1)*T+1):s*T), data_sec(2,((s-1)*T+1):s*T));
    legend('obs_y','data_y');
    xlabel('t');
    
    subplot(3,1,3);
    plot(time(((s-1)*T+1):s*T), observed(3,((s-1)*T+1):s*T)); hold on
    plot(time(((s-1)*T+1):s*T), data_sec(3,((s-1)*T+1):s*T));
    legend('obs_z','data_z');
    xlabel('t');
    
    filename = [pathOut '/comp/comp-slices-' num2str(slices) '-s' num2str(s)];
    saveas(gcf, [filename '.png'],'png');
    saveas(gcf, [filename '.fig']);
    close(figS);
    
end



