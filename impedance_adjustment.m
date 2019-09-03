function [adjusted, locs1_] = impedance_adjustment(x1, x2, sr, mpp, n1, n2, stride, Nmed, locs1, locs2)
% Input:
% x1       - unipolar channel 1
% x2       - unipolar channel 2
% sr       - sampling rate in Hz
% mpp      - MinPeakProminence for findpeaks
% n1       - samples before the peak for the window
% n2       - samples after the peak for the window
% stride   - number of peaks to process on each iteration of the loop
% locs1    - predetermined locations of peaks for x1 if available
% locs2    - predetermined locations of peaks for x2 if available
%
% Output:
% adjusted - adjusted x1 that should more closely match x2
% locs1_   - locations of the associated peaks between the two electrodes

assert(exist('istft', 'file') == 2, sprintf(['istft required\n' ...
    'https://www.mathworks.com/matlabcentral/fileexchange/' ...
    '45577-inverse-short-time-fourier-transform-istft-with-matlab']));

% default n1
if ~exist('n1', 'var') || isempty(n1)
    n1 = round(0.005*sr);
end

% default n2
if ~exist('n2', 'var') || isempty(n2)
    n2 = round(0.010*sr);
end

% ensure that the window length is odd
if mod(n1 + n2, 2) == 1
    n1 = n1 + 1;
end

% the stride determines how many consecutive peaks to apply the same correction
if ~exist('stride', 'var') || isempty(stride)
    stride = 1;
end

% number of windows to average (median) to make robust to outliers
if ~exist('Nmed', 'var') || isempty(Nmed)
    Nmed = 9;
end

%% find peaks
if ~exist('locs1', 'var') || isempty(locs1)
    fprintf(1, 'finding peaks in x1...');
    [~, locs1] = findpeaks(x1, 'MinPeakDistance', .05 * sr, 'MinPeakProminence', mpp);  % , 'MinPeakWidth', .001 * sr, 'MaxPeakWidth', .005 * sr);
    fprintf(1, 'done\n');
end

if ~exist('locs2', 'var') || isempty(locs2)
    fprintf(1, 'finding peaks in x2...');
    [~, locs2] = findpeaks(x2, 'MinPeakDistance', .05 * sr, 'MinPeakProminence', mpp);  % , 'MinPeakWidth', .001 * sr, 'MaxPeakWidth', .005 * sr);
    fprintf('done\n');
end

%% associate peaks
fprintf(1, 'associating peaks...');
[locs1_, ~] = associate_points(locs1, locs2, 0.1);
fprintf('done\n');

%% extract the detected cardiac waveforms
% remove peaks near the edges
locs1_ = locs1_(locs1_ > n1 & locs1_ <= length(x1) - n2);
% extract the data windows at each peak from each unipolar channel
y = cell2mat(arrayfun(@(ii) x1(locs1_(ii) + (-n1:n2)), 1:length(locs1_), 'UniformOutput', false));
y2 = cell2mat(arrayfun(@(ii) x2(locs1_(ii) + (-n1:n2)), 1:length(locs1_), 'UniformOutput', false));

%% compute the FFTs of the cardiac waveforms
f1 = fft(y, [], 1);
f2 = fft(y2, [], 1);
fftratio = f1./f2;
fftratio(1, :) = 1;  % Don't adjust the DC offset to avoid time aliasing

% movmedian should make fftratio robust to outliers
fftratio2 = movmedian(fftratio, Nmed, 2);

%% interpolate the ratio of the FFTs to each sample on an interval
% initialize the output
adjusted = zeros(size(x1));

% instantiate a message updater object to display the status during processing
hmsg = MessageUpdater();

% get the initial time for the status
t0 = clock();

% compute the window length
wlen = n1+n2+1;

% According to Table 1 in http://www.temjournal.com/content/81/TEMJournalFebruary2019_56_64.pdf
% I can determine the hop size that meets the weighted overlap add constraint and provide a speed improvement
hop = floor(wlen/14);

% the number of FFT samples is set to the window length
nfft = 1*wlen;

% use the synergy-complementary windows that were used in the example code 
% that came with Hristo Zhivomirov's istft function
% note istft was introduced in the Matlab Signal Processing Toolbox in R2019a
anal_win = blackmanharris(wlen, 'periodic');
synth_win = hamming(wlen, 'periodic');

for iseg = 1:stride:floor(size(y, 2))
    hmsg.update_message(sprintf('peak %d of %d in %f sec\n', iseg, size(y, 2), etime(clock(), t0)));
    bnd1 = max(1, iseg - 1);
    bnd2 = min(iseg + stride + 1, length(locs1_));
    % TODO can be made more efficient and more conducive to real-time operation by using hop for overlap instead of a whole peak
    inds = locs1_(bnd1):locs1_(bnd2)-hop+1;   % indices for the input data (intervals overlap to avoid edge affects)
    inds2 = locs1_(bnd1+1):locs1_(bnd2)-hop+1; % indices for the result
    
    % apply the analysis filter
    s = stft(x1(inds), anal_win, hop, nfft, sr);

    myinds = locs1_(bnd1):hop:locs1_(bnd2) - n1 - n2 - hop + 1;
    assert(size(s, 2) == length(myinds));
    bnds = interp1(myinds, 1:length(myinds), locs1_(bnd1:bnd2 - 1), 'nearest');
    % apply the correction to each interval
    for ii = 1:length(bnds)-1
        sinds = bnds(ii):bnds(ii+1)-1;
        ib = bnd1 + ii - 1;
        s(:, sinds) = s(:, sinds) ./ fftratio2(1:(wlen-1)/2+1, ib);
    end
    % the last interval goes to the end of the data
    sinds = bnds(end):size(s, 2);
    ib = ib + 1;
    s(:, sinds) = s(:, sinds) ./ fftratio2(1:(wlen-1)/2+1, ib);

    % apply the synthesis filter
    temp4 = istft(s, anal_win, synth_win, hop, nfft, sr);

    % recover the overlapped time-domain signal    
    if iseg == 1        
        adjusted(inds(1:length(temp4))) = temp4;
    else
        % recompute one interval to reduce edge effects but don't overwrite it        
        tempx = temp4(locs1_(bnd1+1) - locs1_(bnd1) + 1:end);
        adjusted(inds2(1:length(tempx))) = tempx;
    end
end

% delete the message updater object
delete(hmsg);
