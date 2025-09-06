function [] = plot_spectral(V)
Fs = 20000; % Sampling frequency (Hz)
N = size(V, 3); % Number of samples per recording (40000)
nChannels = size(V, 2); % Number of recordings (140)
halfPoints = N / 2; % Number of points in each half (20000)
% Parameters for spectrogram computation
windowSize = 1024; % Size of the FFT window
overlap = windowSize / 2; % Overlap between windows
nfft = 2048; % Number of FFT points
freqLimit = 1000; % Limit frequency range to 1000 Hz

psdFirstHalf = [];
psdSecondHalf = [];
avgSpecFirstHalf = [];
avgSpecSecondHalf = [];
for i = 1:nChannels
    % Extract the current recording
    signal = squeeze(V(1, i, :));
    firstHalf = signal(1:halfPoints);
    secondHalf = signal(halfPoints+1:end);

    % Compute PSD for each half using Welch's method
    [Pxx1, F] = pwelch(firstHalf, [], [], [], Fs);
    [Pxx2, ~] = pwelch(secondHalf, [], [], [], Fs);
    validFreqIdx = F <= freqLimit;
    Pxx1 = Pxx1(validFreqIdx);
    Pxx2 = Pxx2(validFreqIdx);
    psdFirstHalf(:, i) = Pxx1;
    psdSecondHalf(:, i) = Pxx2;

    % Compute spectrograms
    [S1, Fspec, T1] = spectrogram(firstHalf, windowSize, overlap, nfft, Fs);
    [S2, ~, T2] = spectrogram(secondHalf, windowSize, overlap, nfft, Fs);
    validFreqIdxSpec = Fspec <= freqLimit;
    S1 = abs(S1(validFreqIdxSpec, :));
    S2 = abs(S2(validFreqIdxSpec, :));
    if isempty(avgSpecFirstHalf)
        avgSpecFirstHalf = S1;
        avgSpecSecondHalf = S2;
    else
        avgSpecFirstHalf = avgSpecFirstHalf + S1;
        avgSpecSecondHalf = avgSpecSecondHalf + S2;
    end
end

% Average PSDs across all recordings
avgPsdFirstHalf = mean(psdFirstHalf, 2);
avgPsdSecondHalf = mean(psdSecondHalf, 2);
% Average the spectrograms across all recordings
avgSpecFirstHalf = avgSpecFirstHalf / nChannels;
avgSpecSecondHalf = avgSpecSecondHalf / nChannels;

figure;
subplot(1, 3, 1);
plot(squeeze(mean(V, 2)))

subplot(1, 3, 2);
plot(F(validFreqIdx), 10*log10(avgPsdFirstHalf), 'b', 'LineWidth', 1.5); hold on;
plot(F(validFreqIdx), 10*log10(avgPsdSecondHalf), 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Average Power Spectral Density (PSD)');
legend({'First Half', 'Second Half'});
grid on;

subplot(1, 3, 3);
T2 = T2 + T1(end);
imagesc([T1 T2], Fspec(validFreqIdxSpec), [10*log10(avgSpecFirstHalf), 10*log10(avgSpecSecondHalf)]);
axis xy;
colormap jet;
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Average Spectrogram of First and Second Halves');
hold on;
% Add dotted line to separate halves
xline(T1(end), 'w--', 'LineWidth', 2);
% Adjust visualization
% caxis([-100, 0]); % Adjust as necessary based on your data
hold off;
end

