clear all
close all

% Goal: Find the ASCII Message
% LFSR Reference: http://in.ncu.edu.tw/ncume_ee/digilogi/prbs.htm

load('Rcvd_Jaro.mat', 'Rcvd')
scatterplot(Rcvd)
title('Constellation Plot of Encoded Message')

%% Filter recieved signal with RRC Filter
beta = 0.75; % RRC rolloff factor
B_RCOS = [0.0038 0.0052 -0.0044 -0.0121 -0.0023 0.0143 0.0044 -0.0385...
    -0.0563 0.0363 0.2554 0.4968 0.6025 0.4968 0.2554 0.0363 -0.0563 -0.0385...
    0.0044 0.0143 -0.0023 -0.0121 -0.0044 0.0052 0.0038];
filtered = filter(B_RCOS,1,Rcvd);

%% Generate m sequence
mSeq = mSeqGen();

%% Oversample m sequence
upMseq = upsample(mSeq,4);

% figure
% plot((filter(fliplr(upMseq),1,repmat(mSeq,5,1))));
% This plot looks like a noisy line with positive slope

% figure
% plot(filter(fliplr(upMseq),1,real(filtered)));
% Bad correlation

figure
plot((xcorr(upMseq,real(filtered))));
title('Correlation: Oversampled M Sequence with Filtered Rcvd')
% first impulse at x = 576, y = -67; next at x = 1596, y = -101

%% Downsample filtered signal
downsampled = downsample(filtered,4);

%% Find and Correct Offset
figure
plot((xcorr(mSeq,real(downsampled(1:255)))));
title('Correlation: M Sequence with Filtered Rcvd Length of 1 Pilot (255)')
% There appears to be a offset of 144 
shiftmSeq = circshift(mSeq,144);

% figure
% plot(xcorr(shiftmSeq,real(downsampled(1:255))));
% First impulse at 255, so I assume the offset is corrected (?)

%% Apply m Sequence
nFrame = length(downsampled)/255;
longshiftmSeq = repmat(shiftmSeq,nFrame,1);
post_m = longshiftmSeq'.*downsampled;
scatterplot(post_m)
title('Constellation Plot of Signal with M Sequence Applied')

%% Phase Shift
% The constellation plot of the signal after the m-sequence is applied is a
% circle. The points are brought down to the x-axis.
sign = ((real(post_m)>0) - 0.5).*2;
freqAdjusted = sign.*(abs(post_m));
scatterplot(freqAdjusted)
title('Frequency-Adjusted Signal')

%% Shift to zeroes and ones
% The signal is now spread across the real line. The following moves points
% closest to 0 to 0, points closest to -1 to -1, and points closest to 1 to
% 1.
for i = 1:length(freqAdjusted)
    x = freqAdjusted(i);
    if (x > -.5) && (x < 0.5)
        freqAdjusted(i)=0;
    elseif x >0.5
        freqAdjusted(i)=1;
    elseif x<-0.5
        freqAdjusted(i)=-1;
    end
end
binaryAdjusted = freqAdjusted;
%scatterplot(binaryAdjusted) % Dots at -1, 0, 1

%% Reshape into an 8x510 matrix
% The new signal must be reshaped for the hadamard matrix to be applied. 
signalMat8 = reshape(binaryAdjusted,[8,510]);

%% Generate 8-ary matrix 
% The hadamard function generates the hadamard matrix
hMat=hadamard(8);

%% Multiple 8-ary Hadamard Matrix by 8x510 Signal Matrix
out = hMat*signalMat8;
data = out(6,:); % channel 5
%scatterplot(data) % data scattered at integers along on x-axis

%% Convert to -1 and 1 
for i = 1:length(data)
    x=data(i);
    if x <= 0 
        data(i) = 1;
    elseif x > 0
        data(i) = -1;
%     else
%         data(i) = 0;
    end
end
preBPSK = data;
%scatterplot(preBPSK) % Dots at -1, 0, 1

%% BPSK Demod
for i = 1:length(preBPSK)
    x = preBPSK(i);
    if x == -1
        preBPSK(i) = 1;
    elseif x == 1
        preBPSK(i) = 0;
    end
end
demodData = preBPSK;
%scatterplot(demodData) % Dots at 0, 1

%% Obtain Characters
out = reshape(transpose(demodData(1:504)),[8 63]); % divisible by 8
out = transpose(out);
for i = 1:63
    c(i) = char(bi2de(out(i,:),'left-msb'));
end

SecretMessage = c % Prints SecretMessage to command window
