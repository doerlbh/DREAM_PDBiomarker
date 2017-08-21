% stage characterization with HMM
% acc_xyz_cont_HMM_resize(note, slices, nstates, resize_length)
% Baihan Lin
% Columbia University
% July 2017

% function acc_xyz_cont_HMM_resize_New(note, slices, nstates, resize_length)
function acc_xyz_cont_HMM_resize_New(note, slices, single_nstates)
% acc_xyz_cont_HMM_resize(note, slices, nstates)
disp('single nstate must be even. Abort if odd.')

nstates = single_nstates + 2; % initial state and end state and two cases

d = 3;

% clear all; close all;

% path = '/gscratch/stf/sunnylin/Columbia/DREAM_PDBiomarker/';
% path = '/home/sunnylin/Dropbox/Git/DREAM_PDBiomarker/';
path = '/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/';

addpath(path);
addpath([path 'jsonlab-master']);

system(['mkdir ' path 'output/' date]);

% slices = 1;
% note = 'nstates-10-slices-1';
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

timeAll = [];
if strcmp(note, 'demo')
    
    test_file = 'test_walk_outbound.tmp';
    
    rawData = loadjson(test_file,'SimplifyCell',1);
    datasize = length(rawData);
    time = zeros(1,datasize);
    acc = zeros(3,datasize);
    
    for t = 1:datasize
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
    
    [sections, mu, sigma, data] = sliceStep(data, d, slices);
    steps = slices;
    
else
    
    files = strsplit(list);
    ncases = length(files)-1;
    data = [];
    %     minsec = 10000000;
    for n = 1:ncases
        file = files{n};
        
        rawData = loadjson(file,'SimplifyCell',1);
        datasize = length(rawData);
        %         time = zeros(1,size);
        %         acc = zeros(3,size);
        
        rawData = loadjson(test_file,'SimplifyCell',1);
        datasize = length(rawData);
        time = zeros(1,datasize);
        acc = zeros(3,datasize);
        
        for t = 1:datasize
            time(1,t) = rawData(1,t).timestamp;
            acc(1,t) = rawData(1,t).x;
            acc(2,t) = rawData(1,t).y;
            acc(3,t) = rawData(1,t).z;
        end
        
        %         time = imresize(time, [1 resize_length]);
        %         acc = imresize(acc, [3 resize_length]);
        
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
        timeAll = [timeAll; time];
        [sections, mu, sigma, dataUncomb] = sliceStep(dataUncomb, d, slices);
        
        %         if secN < minsec
        %             minsec = secN;
        %         end
        data = [data; dataUncomb];
    end
    %     sections = minsec;
    steps = length(data);
    
end

%% train HMM

setSeed(0);

pPrior = [1 zeros(1, nstates-1)].';
% pPrior = ones(1, nstates);
% size(pPrior)

% emPrior.mu = ones(1, d);
% emPrior.Sigma = 0.1*eye(d);
% emPrior.k = d;
% emPrior.dof = emPrior.k + 1;

% emPrior.mu = [zeros(1,d); mu; zeros(1,d)];
% emPrior.Sigma = diag([zeros(1,d); sigma; zeros(1,d)]);
% emPrior.k = d;
% emPrior.dof = emPrior.k + 1;

emPrior.mu = ones(1, d);
emPrior.Sigma = 0.1*eye(d);
emPrior.k = d;
emPrior.dof = emPrior.k + 1;

trPrior = zeros(nstates, nstates);
for tit = 1:nstates-1
    trPrior(tit,tit:tit+1) = 1/2;
end
trPrior(nstates,nstates) = 1;
trPrior([single_nstates/2+1 3*single_nstates/2+1],nstates) = 1/2;
trPrior([single_nstates+1 2*single_nstates+1],2) = 1/2;
trPrior([single_nstates/2+1 3*single_nstates/2+1 single_nstates+1 ...
    2*single_nstates+1],:) = 2/3*trPrior([single_nstates/2+1 ...
    3*single_nstates/2+1 single_nstates+1 2*single_nstates+1],:);
trPrior(1,[1 2 single_nstates/2+2 single_nstates+2 3*single_nstates/2+2]) = 1/5;
% size(trPrior)

nRndRest = 10;
maxIt = 30;

model = hmmFitEmGauss(data, nstates, 'gauss', ...
    'piPrior', pPrior, ...
    'emissionPrior', emPrior,...
    'transPrior', trPrior, ...
    'nRandomRestarts', nRndRest,...
    'maxIter', maxIt);

T = int8(resize_length/slices);

[observed, hidden] = hmmSample(model, T, 1);
% figure;


%% plot clustering

fig2 = figure(ncases+1); hold on
title(['stage characterizations of walking outbound acc']);
[~, colors, ~, ~] =  plotColors();

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
    plot(timeAll(((s-1)*T+1):s*T), observed(1,1:T)); hold on
    plot(timeAll(((s-1)*T+1):s*T), data_sec(1,1:T));
    legend('obs_x','data_x');
    xlabel('t');
    
    subplot(3,1,2);
    plot(timeAll(((s-1)*T+1):s*T), observed(2,1:T)); hold on
    plot(timeAll(((s-1)*T+1):s*T), data_sec(2,1:T));
    legend('obs_y','data_y');
    xlabel('t');
    
    subplot(3,1,3);
    plot(timeAll(((s-1)*T+1):s*T), observed(3,1:T)); hold on
    plot(timeAll(((s-1)*T+1):s*T), data_sec(3,1:T));
    legend('obs_z','data_z');
    xlabel('t');
    
    filename = [pathOut '/comp/comp-slices-' num2str(slices) '-s' num2str(s)];
    saveas(gcf, [filename '.png'],'png');
    saveas(gcf, [filename '.fig']);
    close(figS);
    
end

end

