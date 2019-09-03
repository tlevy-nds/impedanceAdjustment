% generate figures for the paper

%% Figure 1
% schematic

%% Plexon/Cortec
% Use entire recordings to reproduce the figures
% filenames = {'mouse1\File1.mat', ...
%     'mouse2\File1.mat', ...
%     'mouse3\File1.mat'};
% specify the names of the variables that represent the two contacts on the cuff electrode
% myvars = {{'y1', 'y2'}, ...
%     {'y1', 'y2'}, ...
%     {'y3', 'y4'}};

% Use the provided smaller segment of data to run the code
filenames = {'ecgArtifacts.mat'};
myvars = {{'y1', 'y2'}};

sr = 40e3;
dt = 1/sr;

% load the data files
for ifile = 1:length(filenames)
    myvars_ = myvars{ifile};
    temp = load(filenames{ifile}, myvars_{:});
    data{ifile} = 20*double(struct2array(temp));  % Our recordings have a gain of 50, so multiplying by 20 converts to ?V
end

% run the impedance matching algorithm
% load('mymatfile.mat', 'adjusted', 'locs');
if ~exist('adjusted', 'var') || ~exist('locs', 'var')
    % Don't rerun the impedance matching algorithm if the results have been loaded from a mat file
    minPeakProminence = 500;  % parameter in findpeaks
    n1 = 300;                 % number of samples before the ECG artifact peak to use
    n2 = 300;                 % number of samples after the ECG artifact peak to use
    stride = 1;               % number of ECG artifacts to skip before adjusting the correction
    nmed = 9;                 % number of ECG artifacts for smoothing the correction
    for ifile = 1:length(filenames)
        duration = size(data{ifile}, 1) / sr;
        
        % impedance matching algorithm
        t0 = clock();
        [adjusted{ifile}, locs{ifile}] = impedance_adjustment( ...
            data{ifile}(:, 1), data{ifile}(:, 2), sr, ...
            minPeakProminence, n1, n2, stride, nmed);
        tf = clock();
        
        % measure the time this takes and print it
        fprintf(1, ['ran impedance matching algorithm on %f seconds of' ...
            ' data sampled at %f Hz in %f seconds\n'], duration, sr, etime(tf, t0));
    end
    % save('mymatfile.mat', 'adjusted', 'locs', '-v7.3');
end

% Power distributions (Figure 6)
INRuni=[];INRbi=[];INRadj=[];
for ifile = 1:length(filenames)
    % ECG interferernce signals for each method to compare
    uni = mean(data{ifile}, 2)/1e6;                   % mean unipolar in V
    bi = diff(data{ifile}, 1, 2)/1e6;                 % bipolar derivation in V
    adj = (data{ifile}(:, 2) - adjusted{ifile})/1e6;  % impedance adjusted bipolar derivation in V
    
    % mean ECG interference power
    Suni = arrayfun(@(x) mean(uni(x-300:x+300).^2), locs{ifile});   % mean power of mean unipolar
    Sbi = arrayfun(@(x) mean(bi(x-300:x+300)'.^2), locs{ifile});    % mean power of bipolar
    Sadj = arrayfun(@(x) mean(adj(x-300:x+300)'.^2), locs{ifile});  % mean power of impedance adjusted bipolar

    % mean noise power (The noise is measured from the regions between the ECG artifact windows)
    Nuni = arrayfun(@(xm1, x, xp1) mean(uni([xm1+300:x-300 x+300:xp1-300])'.^2), ...
        [0; locs{ifile}(1:end-1)], locs{ifile}, [locs{ifile}(2:end); length(uni)]);
    Nbi = arrayfun(@(xm1, x, xp1) mean(bi([xm1+300:x-300 x+300:xp1-300])'.^2), ...
        [0; locs{ifile}(1:end-1)], locs{ifile}, [locs{ifile}(2:end); length(bi)]);
    Nadj = arrayfun(@(xm1, x, xp1) mean(adj([xm1+300:x-300 x+300:xp1-300])'.^2), ...
        [0; locs{ifile}(1:end-1)], locs{ifile}, [locs{ifile}(2:end); length(adj)]);
    
    % interference to noise ratio (INR)
    INRuni = [INRuni; Suni./Nuni];
    INRbi = [INRbi; Sbi./Nbi];
    INRadj = [INRadj; Sadj./Nadj];
end
% INR histograms
[funi, xiuni] = ksdensity(10*log10(INRuni));
[fbi, xibi] = ksdensity(10*log10(INRbi));
[fadj, xiadj] = ksdensity(10*log10(INRadj));

%% Figure 3
figure;ax = gca;
hold on;
plot(xiuni, funi, 'LineWidth', 1);
plot(xibi, fbi, 'LineWidth', 1);
plot(xiadj, fadj, 'LineWidth', 1);
hold off;
set(ax, 'XScale', 'linear');
xlabel(ax, 'INR (dB)');
ylabel(ax, 'Likelihood');
title(ax, 'ECG INR Distributions');
legend(ax, {'Mean Unipolar' 'Simple Subtraction' 'Impedance Adjusted Subtraction'}, ...
    'Location', 'northeast');

% Test for significance
p = signrank(10*log10(INRuni), 10*log10(INRbi), 'Tail', 'right');  % p = 0
p = signrank(10*log10(INRbi), 10*log10(INRadj), 'Tail', 'right');  % p = 0

%% Figure 2A and Figure 2B
ifile = 1;
t = (0:size(data{ifile}, 1)-1)'/sr;
inds = find(t >= 0.90 & t <= 1.44);
% unipolar
figure;plot(t(inds), [data{ifile}(inds, :) diff(data{ifile}(inds, :), 1, 2)], 'LineWidth', 1);
xlabel('Time (s)');ylabel('Voltage (\muV)');legend({'Channel 1', 'Channel 2', 'Simple Subtraction'});
xlim([0.90 1.44]);ylim([-284.5 752.8]);  % 5 ECG artifacts
% xlim([1.36371664579075, 1.3763536985019]);ylim([-253.219552529183, 586.307295719844]);  % zoomed in on last ECG artifact
% bipolar
figure;plot(t(inds), [diff(data{ifile}(inds, :), 1, 2), data{ifile}(inds, 2) - adjusted{ifile}(inds)], 'LineWidth', 1);
xlabel('Time (s)');ylabel('Voltage (\muV)');legend({'Simple Subtraction', 'Impedance Adjusted Subtraction'});
xlim([0.90 1.44]);ylim([-31.2 28.8]);  % 5 ECG artifacts
% xlim([1.36371664579075, 1.3763536985019]);  % zoomed in on last ECG artifact

% Annotations for figure 2A and 2B
nlocs = [0 cumsum(cellfun(@(x) length(x), locs))];
ind = nlocs(ifile) + find(t(locs{ifile}) >= 1.3, 1, 'first');  % Figure 2 shows ECG events 9-13 
disp(10*log10([INRuni(ind) INRbi(ind) INRadj(ind)]));

%% Figure 2C Plexon/Cortec time-frequency (ECG removal)
% unipolar channel 1
[s, f, t] = spectrogram(data{ifile}(:, 1), 5*sr, 0, 5*sr, sr);
figure;imagesc('Xdata', t, 'YData', f, 'CData', 20*log10(abs(s)));
axis image;axis normal;axis ij;colorbar;xlabel('Time (s)');ylabel('Frequency (Hz)');
caxis([0 120]);ylim([0 1000]);

% bipolar
[s, f, t] = spectrogram(data{ifile}(:, 2) - data{ifile}(:, 1), 5*sr, 0, 5*sr, sr);
figure;imagesc('Xdata', t, 'YData', f, 'CData', 20*log10(abs(s)));
axis image;axis normal;axis ij;colorbar;xlabel('Time (s)');ylabel('Frequency (Hz)');
caxis([0 100]);ylim([0 1000]);

% adjusted bipolar
[s, f, t] = spectrogram(data{ifile}(:, 2) - adjusted{ifile}, 5*sr, 0, 5*sr, sr);
figure;imagesc('Xdata', t, 'YData', f, 'CData', 20*log10(abs(s)));
axis image;axis normal;axis ij;colorbar;xlabel('Time (s)');ylabel('Frequency (Hz)');
caxis([0 100]);ylim([0 1000]);

%% Intan/Flex
% clear  % Load new files so previous data can be cleared if desired
% Use the entire recordings to reproduce the figures
% dirnames = {'rat1\2018-12-03-uy-1', ...
%     'rat2\2018-12-06-u-1', ...
%     'rat3\2019-01-04-u-1'};
% % data files to use in each path
% datafiles = {{'exp-2018-12-03-uy-1-1_181203_140142.mat', 'exp-2018-12-03-uy-1-1_181203_143142.mat'}, 'exp-2018-12-03-uy-1-3_181203_151408.mat', 'exp-2018-12-03-uy-1-4_181203_161216.mat'; ...
%         'exp-2018-12-06-u-1-1_181206_142732.mat', 'exp-2018-12-06-u-1-2_181206_150710.mat', 'exp-2018-12-06-u-1-3_181206_162501.mat'; ...
%         'exp-2019-01-04-uy-1-1_190104_131456.mat', 'exp-2019-01-04-uy-1-2_190104_141532.mat', 'exp-2019-01-04-uy-1-3_190104_150854.mat'};
% % stimulation data to load for each path
% matfiles = {'rat1\2018-12-03-uy-1\2018-12-3-14-6-uy_randpulse_SEQ.mat', ...
%     'rat2\2018-12-06-u-1\2018-12-6-14-25-u_randpulse_SEQ.mat', ...
%     'rat3\2019-01-04-u-1\2019-1-4-13-8-u_randpulse_SEQ.mat'};
% chans = {[4 5], [4 5], [4 5]};  % data channels from the files for each experiment directory
% fig4dirInd = 3;
% fig6dirInd = 3;

% Use the provided smaller segment of data to run the code
segnames = {'Before VEC', 'During VEC', 'After VEC'};
sr = 30e3;
dt = 1/sr;

clear('data');
if exist('stimulationArtifacts.mat', 'file')
    % load the sample data
    dirnames = {''};  % need something of length 1
    load('stimulationArtifacts.mat', 'data', 'rseq', 'WID');
    fig4dirInd = 1;
    fig6dirInd = 1;
end
if ~exist('data', 'var')
    % load the entire recordings
    [rseq, WID] = deal(cell(1, length(dirnames)));
    for idir = 1:length(dirnames)
        % hut{idir} = UTExperiment(dirnames{idir}, [], 'Yes', '', true, true);
        data{idir} = struct('trigchan', [], 'neuraldata', []);
        data{idir}.trigchan = cell(1, length(segnames));
        data{idir}.neuraldata = cell(1, length(segnames));
        for iseg = 1:length(segnames)
            if iscell(datafiles{idir, iseg})
                % This segment is composed of multiple data files
                trigchan = [];
                neuraldata = [];
                for ifile = 1:length(datafiles{idir, iseg})
                    load([dirnames{idir} filesep datafiles{idir, iseg}{ifile}], 'board_dig_in_data', 'amplifier_data');
                    data{idir}.trigchan{iseg} = [data{idir}.trigchan{iseg}; board_dig_in_data(2, :)'];
                    data{idir}.neuraldata{iseg} = [data{idir}.neuraldata{iseg}; amplifier_data(chans{idir}, :)'];
                end
            else
                % This segment is composed of a single data file
                load([dirnames{idir} filesep datafiles{idir, iseg}], 'board_dig_in_data', 'amplifier_data');
                data{idir}.trigchan{iseg} = board_dig_in_data(2, :)';
                data{idir}.neuraldata{iseg} = amplifier_data(chans{idir}, :)';
            end
        end
        
        temp = load(matfiles{idir}, 'rseq', 'WID');
        rseq{idir} = temp.rseq;  % each set of stimulation parameters is repeated
        WID{idir} = temp.WID;    % unique sets of stimulation parameters
    end
    clear('board_dig_in_data', 'amplifier_data', 'temp');
end

%% determine location of stimulation artifacts
% get connected components for the digital trigger channel
[ccBefore, ccVec, ccAfter] = deal(cell(1, length(dirnames)));
for idir = 1:length(dirnames)
    ccVec{idir} = bwconncomp(data{idir}.trigchan{2}, 4);  % vecuronium
    ccBefore{idir} = bwconncomp(data{idir}.trigchan{1}, 4);  % before
    ccAfter{idir} = bwconncomp(data{idir}.trigchan{3}, 4);  % after
end

% first sample in each pulse
[sBefore, sVec, sAfter] = deal(cell(1, length(dirnames)));
for idir = 1:length(dirnames)
    sBefore{idir} = cellfun(@(x) x(1), ccBefore{idir}.PixelIdxList);
    sVec{idir} = cellfun(@(x) x(1), ccVec{idir}.PixelIdxList);
    sAfter{idir} = cellfun(@(x) x(1), ccAfter{idir}.PixelIdxList);
end

% band-pass filtered neural data
[b, a] = butter(3, [10 10e3]/(sr/2));
[fndBefore, fndVec, fndAfter] = deal(cell(1, length(dirnames)));
for idir = 1:length(dirnames)
    fndBefore{idir} = filtfilt(b, a, double(data{idir}.neuraldata{1}));  % before
    fndVec{idir} = filtfilt(b, a, double(data{idir}.neuraldata{2}));  % vecuronium
    fndAfter{idir} = filtfilt(b, a, double(data{idir}.neuraldata{3}));  % after
end

%% Run the impedance matching algorithm
n1 = 20;
n2 = 20;

if ~exist('adjustedBefore', 'var')
    adjustedBefore = cell(1, length(dirnames));
    for idir = 1:length(dirnames)
        adjustedBefore{idir} = impedance_adjustment( ...
            fndBefore{idir}(:, 1), fndBefore{idir}(:, 2), sr, ...
            [], 20, 20, 1, 9, sBefore{idir}', sBefore{idir}');
    end
end

if ~exist('adjustedVec', 'var')
    adjustedVec = cell(1, length(dirnames));
    for idir = 1:length(dirnames)
        adjustedVec{idir} = impedance_adjustment( ...
            fndVec{idir}(:, 1), fndVec{idir}(:, 2), sr, ...
            [], 20, 20, 1, 9, sVec{idir}', sVec{idir}');
    end
end

if ~exist('adjustedAfter', 'var')
    adjustedAfter = cell(1, length(dirnames));
    for idir = 1:length(dirnames)
        adjustedAfter{idir} = impedance_adjustment( ...
            fndAfter{idir}(:, 1), fndAfter{idir}(:, 2), sr, ...
            [], 20, 20, 1, 9, sAfter{idir}', sAfter{idir}');
    end
end

% bipolar simple subtraction
% [bi1, bi2, bi3] = deal(cell(1, length(dirnames)));
% for idir = 1:length(dirnames)
%     bi1{idir} = diff(fndBefore{idir}, 1, 2);
%     bi2{idir} = diff(fndVec{idir}, 1, 2);
%     bi3{idir} = diff(fndAfter{idir}, 1, 2);
% end

% unipolar
[uniBefore, uniVec, uniAfter] = deal(cell(1, length(dirnames)));
[uni2Before, uni2Vec, uni2After] = deal(cell(1, length(dirnames)));
for idir = 1:length(dirnames)
    uniBefore{idir} = fndBefore{idir}(:, 1);
    uni2Before{idir} = fndBefore{idir}(:, 2);
    uniVec{idir} = fndVec{idir}(:, 1);
    uni2Vec{idir} = fndVec{idir}(:, 2);
    uniAfter{idir} = fndAfter{idir}(:, 1);
    uni2After{idir} = fndAfter{idir}(:, 2);
end

%% Figure 4a (Stim artifact removal) 
% TODO this won't run unless I include all of the stimulation pulses
% go no farther
return

% get channel 1, channel 2, and adjusted channel each segment of the ENG with a stimulation artifact
nparams = length(WID{1});
ndirs = length(dirnames);
nsegs = 3;
[nLeft, nRight] = deal(10, 600);            % Window parameters
nWindow = nLeft + nRight + 1;
nRepeats = length(rseq{1}) / length(WID{1});  % each set of stimulation parameters repeats
% byRep's organize the data windows by stimulation repeats so that they can be averaged along that dimension
[byRepUni1, byRepUni2, byRepAdj] = deal( ...
    zeros(nWindow, nRepeats, ndirs, nsegs, nparams));
for idir = 1:length(dirnames)
    for iseg = 1:nsegs
        for iparams = 1:nparams
            inds = (rseq{idir} == iparams);
            if iseg == 1  % The variable name for each seg is different
                byRepUni1(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uniBefore{idir}(x - nLeft:x + nRight), [], 1), ...
                    sBefore{idir}(inds), 'UniformOutput', false));
                byRepUni2(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uni2Before{idir}(x - nLeft:x + nRight), [], 1), ...
                    sBefore{idir}(inds), 'UniformOutput', false));
                byRepAdj(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(adjustedBefore{idir}(x - nLeft:x + nRight), [], 1), ...
                    sBefore{idir}(inds), 'UniformOutput', false));
            elseif iseg == 2
                byRepUni1(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uniVec{idir}(x - nLeft:x + nRight), [], 1), ...
                    sVec{idir}(inds), 'UniformOutput', false));
                byRepUni2(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uni2Vec{idir}(x - nLeft:x + nRight), [], 1), ...
                    sVec{idir}(inds), 'UniformOutput', false));
                byRepAdj(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(adjustedVec{idir}(x - nLeft:x + nRight), [], 1), ...
                    sVec{idir}(inds), 'UniformOutput', false));
            else
                byRepUni1(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uniAfter{idir}(x - nLeft:x + nRight), [], 1), ...
                    sAfter{idir}(inds), 'UniformOutput', false));
                byRepUni2(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uni2After{idir}(x - nLeft:x + nRight), [], 1), ...
                    sAfter{idir}(inds), 'UniformOutput', false));
                byRepAdj(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(adjustedAfter{idir}(x - nLeft:x + nRight), [], 1), ...
                    sAfter{idir}(inds), 'UniformOutput', false));
            end
        end
    end
end

% take the median across stimulation parameter repeats
[idir, iseg, iparams] = deal(fig4dirInd, 1, 47);  % choose a stimulation artifact
h = figure;
t = 1000*((0:size(byRepUni1, 1) - 1) - 10)/sr;
figure(h);plot(t, [median(byRepUni1(:, :, idir, iseg, iparams), 2), ...
    median(byRepUni2(:, :, idir, iseg, iparams), 2)], 'LineWidth', 1);
xlabel('Time (ms)');ylabel('Voltage (\muV)');legend({'Channel 1', 'Channel 2'});
xlim([t(1) t(1) + 1000*0.02]);
title('Stimulation Artifact');

%% Figure 4b (Stim artifact removal)
% get simple subtraction and impedance adjusted subtraction for each segment of the ENG with a stimulation artifact
[simpSub, impAdjSub, mymeanuni] = deal(zeros(nWindow, nRepeats, ndirs, nsegs, nparams));
for idir = 1:length(dirnames)
    for iseg = 1:nsegs
        for iparams = 1:nparams
            inds = (rseq{idir} == iparams);
            if iseg == 1  % The variable name for each seg is different
                simpSub(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uni2Before{idir}(x - nLeft:x + nRight) - ...
                    uniBefore{idir}(x - nLeft:x + nRight), [], 1), ...
                    sBefore{idir}(inds), 'UniformOutput', false));
                impAdjSub(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uni2Before{idir}(x - nLeft:x + nRight) - ...
                    adjustedBefore{idir}(x - nLeft:x + nRight), [], 1), ...
                    sBefore{idir}(inds), 'UniformOutput', false));
            elseif iseg == 2
                simpSub(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uni2Vec{idir}(x - nLeft:x + nRight) - ...
                    uniVec{idir}(x - nLeft:x + nRight), [], 1), ...
                    sVec{idir}(inds), 'UniformOutput', false));
                impAdjSub(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uni2Vec{idir}(x - nLeft:x + nRight) - ...
                    adjustedVec{idir}(x - nLeft:x + nRight), [], 1), ...
                    sVec{idir}(inds), 'UniformOutput', false));
            else
                simpSub(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uni2After{idir}(x - nLeft:x + nRight) - ...
                    uniAfter{idir}(x - nLeft:x + nRight), [], 1), ...
                    sAfter{idir}(inds), 'UniformOutput', false));
                impAdjSub(:, :, idir, iseg, iparams) = cell2mat( ...
                    arrayfun(@(x) reshape(uni2After{idir}(x - nLeft:x + nRight) - ...
                    adjustedAfter{idir}(x - nLeft:x + nRight), [], 1), ...
                    sAfter{idir}(inds), 'UniformOutput', false));
            end
        end
    end
end

% take the median across stimulation parameter repeats
[idir, iseg, iparams] = deal(fig4dirInd, 1, 47);  % choose the same stimulation artifact
h = figure;
t = 1000*((0:size(simpSub, 1)-1) - 10)/sr;
figure(h);plot(t, [median(simpSub(:, :, idir, iseg, iparams), 2), ...
    median(impAdjSub(:, :, idir, iseg, iparams), 2)], 'LineWidth', 1);
xlabel('Time (ms)');ylabel('Voltage (\muV)');legend({'Simple Subtraction', 'Impedance Adjusted Subtraction'});
xlim([t(1) t(1) + 1000*0.02]);
title(sprintf('%d Amp=%f PW=%d', iparams, WID{idir}(iparams).Amp, WID{idir}(iparams).PW));

%% Figure 6 (EMG uni)

% filter by large amplitudes and long pulse widths
h = figure;
for idir = 1:length(dirnames)  % select a single recording
    inds_ = find([WID{idir}.Amp] > 2 & [WID{idir}.PW] > 200);  % select a subset of the stimulation parameters that elicit a response
    inds = ismember(rseq{idir}, inds_);  % each stimulation parameter is repeated 5 times
    % Channel 1 responses to stimulation artifacts in each segment (Before VEC, during VEC, after VEC)
    responsesUni1Before = cell2mat(arrayfun(@(x) reshape(uniBefore{idir}(x - nLeft:x + nRight), [], 1), ...
        sBefore{idir}(inds), 'UniformOutput', false));
    responsesUni1Vec = cell2mat(arrayfun(@(x) reshape(uniVec{idir}(x - nLeft:x + nRight), [], 1), ...
        sVec{idir}(inds), 'UniformOutput', false));
    responsesUni1After = cell2mat(arrayfun(@(x) reshape(uniAfter{idir}(x - nLeft:x + nRight), [], 1), ...
        sAfter{idir}(inds), 'UniformOutput', false));
    % Channel 2 responses to stimulation artifacts in each segment (Before VEC, during VEC, after VEC)
    responsesUni2Before = cell2mat(arrayfun(@(x) reshape(uni2Before{idir}(x - nLeft:x + nRight), [], 1), ...
        sBefore{idir}(inds), 'UniformOutput', false));
    responsesUni2Vec = cell2mat(arrayfun(@(x) reshape(uni2Vec{idir}(x - nLeft:x + nRight), [], 1), ...
        sVec{idir}(inds), 'UniformOutput', false));
    responsesUni2After = cell2mat(arrayfun(@(x) reshape(uni2After{idir}(x - nLeft:x + nRight), [], 1), ...
        sAfter{idir}(inds), 'UniformOutput', false));
    % Adjusted channel 2 responses to stimulation artifacts in each segment (Before VEC, during VEC, after VEC)
    responsesAdjBefore = cell2mat(arrayfun(@(x) reshape(adjustedBefore{idir}(x - nLeft:x + nRight), [], 1), ...
        sBefore{idir}(inds), 'UniformOutput', false));
    responsesAdjVec = cell2mat(arrayfun(@(x) reshape(adjustedVec{idir}(x - nLeft:x + nRight), [], 1), ...
        sVec{idir}(inds), 'UniformOutput', false));
    responsesAdjAfter = cell2mat(arrayfun(@(x) reshape(adjustedAfter{idir}(x - nLeft:x + nRight), [], 1), ...
        sAfter{idir}(inds), 'UniformOutput', false));
    
    % Figure 6A
    t = 1000*((0:size(responsesUni1Before, 1) - 1) - nLeft)/sr;
    % median response of
    % unipolar channels before VEC are solid lines
    % unipolar channels during VEC are dotted lines
    if idir == fig6dirInd
        figure(h);ax = gca;pl = plot(t, median(responsesUni1Before, 2), '-', t, median(responsesUni2Before, 2), '-', ...
            t, median(responsesUni1Vec, 2), ':', t, median(responsesUni2Vec, 2), ':', 'LineWidth', 1);
        xlim([t(1) t(end)]);grid('on');
        xlabel('Time (ms)');ylabel('Voltage (\muV)');
        title('Median Evoked Responses');
    end
    
    % determine where evoked response is different before and during VEC but the same between channel 1 and channel 2
    % [p1, tbl, stats] = anova1([responses1(:, ii) responses2(:, ii) responses1b(:, ii)]);[c, m] = multcompare(stats);
    p1 = arrayfun(@(ii) anova1([responsesUni1Before(ii, :)' responsesUni1Vec(ii, :)'], [], 'off'), 1:size(responsesUni1Before, 1));
    p2 = arrayfun(@(ii) anova1([responsesUni1Before(ii, :)' responsesUni2Before(ii, :)'], [], 'off'), 1:size(responsesUni1Before, 1));
    inds = p1 < 1e-5 & p2 > 1e-5;
    % I need the loop for inds2
    inds2{idir} = find(inds);  % convert to regular indices for the next figure
    ccv = bwconncomp(inds, 4);
    
    if idir == fig6dirInd
        yl = ylim(ax);
        hold(ax, 'on');
        for icc = 1:length(ccv.PixelIdxList)
            if length(ccv.PixelIdxList{icc}) < 5
                % don't draw really narrow regions
                continue
            end
            t0 = t(ccv.PixelIdxList{icc}(1)) - dt/2;
            tf = t(ccv.PixelIdxList{icc}(end)) + dt/2;
            patch(ax, 'XData', [t0 tf tf t0], 'YData', yl([1 1 2 2]), 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.15);
        end
        hold(ax, 'off');
        legend(pl, {'Channel 1 before VEC', 'Channel 2 before VEC', 'Channel 1 during VEC', 'Channel 2 during VEC'});
    end
end

%% Figure 5b (EMG removal)
t = 1000*((0:size(responsesUni1Before, 1) - 1) - nLeft)/sr;
figure;ax = gca;
pl = plot(t, median(responsesUni2Before - responsesAdjBefore, 2), '-', ...
    t, median(responsesUni2Vec - responsesAdjVec, 2), '-', 'LineWidth', 1);
xlim([t(1) t(end)]);grid('on');
xlabel('Time (ms)');ylabel('Voltage (\muV)');
title('Impedance Adjusted Median Evoked Responses');
yl = ylim(ax);

% add EMG region patches
hold(ax, 'on');
for icc = 1:length(ccv.PixelIdxList)
    t0 = t(ccv.PixelIdxList{icc}(1)) - dt/2;
    tf = t(ccv.PixelIdxList{icc}(end)) + dt/2;
    patch(ax, 'XData', [t0 tf tf t0], 'YData', yl([1 1 2 2]), ...
        'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.15);
end
hold(ax, 'off');
legend(pl, {'Before VEC' 'During VEC'});

%% Figures 5 and 7 (INR distributions)
% TODO what is 10, don't reproduce figures
emg = true;  % false for figure 5, true for figure 7
if emg
    % EMG regions for each experiment
    w = cellfun(@(x) x(x>10) - 10, inds2, 'UniformOutput', false);  % EMG using ANOVA bounds
end

[Suni, Sbi, Sadj] = deal([]);
[Nuni, Nbi, Nadj] = deal([]);
for idir = 1:length(dirnames) 
    % stimulation artifacts that evoke EMG responses
    inds_ = find([WID{idir}.Amp] > 2 & [WID{idir}.PW] > 200);
    inds = find(ismember(rseq{idir}, inds_));
    
    uniBefore_ = mean(fndBefore{idir}, 2)/1e6;
    biBefore_ = diff(fndBefore{idir}, 1, 2)/1e6;
    adjBefore_ = (fndBefore{idir}(:, 2) - adjustedBefore{idir})/1e6;
    if ~emg
        uniVec_ = mean(fndVec{idir}, 2)/1e6;
        biVec_ = diff(fndVec{idir}, 1, 2)/1e6;
        adjVec_ = (fndVec{idir}(:, 2) - adjustedVec{idir})/1e6;
    end
    uniAfter_ = mean(fndAfter{idir}, 2)/1e6;
    biAfter_ = diff(fndAfter{idir}, 1, 2)/1e6;
    adjAfter_ = (fndAfter{idir}(:, 2) - adjustedAfter{idir})/1e6;
    
    if ~emg
        pws = round(sr * arrayfun(@(x) WID{idir}(x).PW * 1e-6, rseq{idir}(inds)))';
        SuniBefore = arrayfun(@(x, pw) median(uniBefore_(x + (0:pw)).^2), sBefore{idir}(inds), pws);
        SbiBefore = arrayfun(@(x, pw) median(biBefore_(x + (0:pw)).^2), sBefore{idir}(inds), pws);
        SadjBefore = arrayfun(@(x, pw) median(adjBefore_(x + (0:pw)).^2), sBefore{idir}(inds), pws);
        SuniVec = arrayfun(@(x, pw) median(uniVec_(x + (0:pw)).^2), sVec{idir}(inds), pws);
        SbiVec = arrayfun(@(x, pw) median(biVec_(x + (0:pw)).^2), sVec{idir}(inds), pws);
        SadjVec = arrayfun(@(x, pw) median(adjVec_(x + (0:pw)).^2), sVec{idir}(inds), pws);
        SuniAfter = arrayfun(@(x, pw) median(uniAfter_(x + (0:pw)).^2), sAfter{idir}(inds), pws);
        SbiAfter = arrayfun(@(x, pw) median(biAfter_(x + (0:pw)).^2), sAfter{idir}(inds), pws);
        SadjAfter = arrayfun(@(x, pw) median(adjAfter_(x + (0:pw)).^2), sAfter{idir}(inds), pws);
    else
        SuniBefore = arrayfun(@(x) median(uniBefore_(x + w{idir}).^2), sBefore{idir}(inds));
        SbiBefore = arrayfun(@(x) median(biBefore_(x + w{idir}).^2), sBefore{idir}(inds));
        SadjBefore = arrayfun(@(x) median(adjBefore_(x + w{idir}).^2), sBefore{idir}(inds));        
        SuniAfter = arrayfun(@(x) median(uniAfter_(x + w{idir}).^2), sAfter{idir}(inds));
        SbiAfter = arrayfun(@(x) median(biAfter_(x + w{idir}).^2), sAfter{idir}(inds));
        SadjAfter = arrayfun(@(x) median(adjAfter_(x + w{idir}).^2), sAfter{idir}(inds));
    end
    
    if ~emg
        Suni = [Suni, SuniBefore, SuniVec, SuniAfter];
        Sbi = [Sbi, SbiBefore, SbiVec, SbiAfter];
        Sadj = [Sadj, SadjBefore, SadjVec, SadjAfter];
    else
        Suni = [Suni, SuniBefore, SuniAfter];
        Sbi = [Sbi, SbiBefore, SbiAfter];
        Sadj = [Sadj, SadjBefore, SadjAfter];
    end
    
    NuniBefore = arrayfun(@(xm1, x, xp1) median(uniBefore_([xm1+2000:x-500 x+2000:xp1-500]).^2), ...
        [0, sBefore{idir}(inds(1:end-1))], sBefore{idir}(inds), [sBefore{idir}(inds(2:end)), length(uniBefore_)]);
    NbiBefore = arrayfun(@(xm1, x, xp1) median(biBefore_([xm1+2000:x-500 x+2000:xp1-500]).^2), ...
        [0, sBefore{idir}(inds(1:end-1))], sBefore{idir}(inds), [sBefore{idir}(inds(2:end)), length(biBefore_)]);
    NadjBefore = arrayfun(@(xm1, x, xp1) median(adjBefore_([xm1+2000:x-500 x+2000:xp1-500]).^2), ...
        [0, sBefore{idir}(inds(1:end-1))], sBefore{idir}(inds), [sBefore{idir}(inds(2:end)), length(adjBefore_)]);
    if ~emg
        NuniVec = arrayfun(@(xm1, x, xp1) median(uniVec_([xm1+2000:x-500 x+2000:xp1-500]).^2), ...
            [0, sVec{idir}(inds(1:end-1))], sVec{idir}(inds), [sVec{idir}(inds(2:end)), length(uniVec_)]);
        NbiVec = arrayfun(@(xm1, x, xp1) median(biVec_([xm1+2000:x-500 x+2000:xp1-500]).^2), ...
            [0, sVec{idir}(inds(1:end-1))], sVec{idir}(inds), [sVec{idir}(inds(2:end)), length(biVec_)]);
        NadjVec = arrayfun(@(xm1, x, xp1) median(adjVec_([xm1+2000:x-500 x+2000:xp1-500]).^2), ...
            [0, sVec{idir}(inds(1:end-1))], sVec{idir}(inds), [sVec{idir}(inds(2:end)), length(adjVec_)]);
    end
    NuniAfter = arrayfun(@(xm1, x, xp1) median(uniAfter_([xm1+2000:x-500 x+2000:xp1-500]).^2), ...
        [0, sAfter{idir}(inds(1:end-1))], sAfter{idir}(inds), [sAfter{idir}(inds(2:end)), length(uniAfter_)]);    
    NbiAfter = arrayfun(@(xm1, x, xp1) median(biAfter_([xm1+2000:x-500 x+2000:xp1-500]).^2), ...
        [0, sAfter{idir}(inds(1:end-1))], sAfter{idir}(inds), [sAfter{idir}(inds(2:end)), length(biAfter_)]);
    NadjAfter = arrayfun(@(xm1, x, xp1) median(adjAfter_([xm1+2000:x-500 x+2000:xp1-500]).^2), ...
        [0, sAfter{idir}(inds(1:end-1))], sAfter{idir}(inds), [sAfter{idir}(inds(2:end)), length(adjAfter_)]);
    
    % don't include Vec for EMG
    if ~emg
        Nuni = [Nuni, NuniBefore, NuniVec, NuniAfter];
        Nbi = [Nbi, NbiBefore, NbiVec, NbiAfter];
        Nadj = [Nadj, NadjBefore, NadjVec, NadjAfter];
    else
        Nuni = [Nuni, NuniBefore, NuniAfter];
        Nbi = [Nbi, NbiBefore, NbiAfter];
        Nadj = [Nadj, NadjBefore, NadjAfter];
    end
end

% combine all 3 segments (before, during, after vecuronium)
[f1, xi1] = ksdensity(10*log10(Suni./Nuni));  % mean unipolar
[f2, xi2] = ksdensity(10*log10(Sbi./Nbi));  % simple subtraction
[f3, xi3] = ksdensity(10*log10(Sadj./Nadj));  % impedance adjusted subtraction
figure;ax = gca;hold on;
plot(xi1, f1, '-', 'LineWidth', 1);
plot(xi2, f2, '-', 'LineWidth', 1);
plot(xi3, f3, '-', 'LineWidth', 1);
hold off;

set(ax, 'XScale', 'linear');
xlabel('Average Power (dBm)');ylabel('Likelihood');
if ~emg
    title('Stimulation Artifact INR Distributions');
else
    title('EMG INR Distributions');
end
legend({'Mean Unipolar' 'Simple Subtraction' 'Impedance Adjusted Subtraction'});

p1 = signrank(10*log10(Suni./Nuni), 10*log10(Sbi./Nbi), 'Tail', 'right');  % sa: p = 1.9e-9,   emg: p = 1
p2 = signrank(10*log10(Sbi./Nbi), 10*log10(Sadj./Nadj), 'Tail', 'right');  % sa: p = 1.9e-187, emg: p = 1.5e-5

%%
% median EMG INR dB for Figure 6A
% w = cellfun(@(x) x(x>10) - 10, inds2, 'UniformOutput', false);  % EMG using ANOVA bounds
temp = [(responsesUni1Before + responsesUni2Before)/2 (responsesUni1Vec + responsesUni2Vec)/2];
S = median((temp(w{fig6dirInd}, :)/1e6).^2, 1);
medEmgInrAvgUni = median(10*log10(S ./ Nuni));

%% median EMG INR dB for Figure 6B
temp = [(responsesUni2Before - responsesAdjBefore) (responsesUni2Vec - responsesAdjVec)];
S = median((temp(w{fig6dirInd}, :)/1e6).^2, 1);
medEmgInrImpAdj = median(10*log10(S ./ Nuni));
