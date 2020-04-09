clc
clear all
close all
load('Rcvd_Jaro.mat', 'Rcvd')

% Filter recieved signal with RRC Filter
B_RCOS = [0.0038 0.0052 -0.0044 -0.0121 -0.0023 0.0143 0.0044 -0.0385 -0.0563 0.0363 0.2554 ...
    0.4968 0.6025 0.4968 0.2554 0.0363 -0.0563 -0.0385 0.0044 0.0143 -0.0023 -0.0121 -0.0044 ...
    0.0052 0.0038];
unfilterRcvd = filter(B_RCOS,1,Rcvd);

% Generate m sequence
G = [1 1 1 0 0 0 0 1]; % G(X)= X^8 + X^7 + X^6 + X + 1
M = mSeqGen(G);

% Oversample the m sequence 
spreadM = (reshape([1-2*M;zeros(3,length(M))],1,[]));

% Correlate oversampled m against oversampled unfiltered Rcvd
plot(abs(filter(fliplr(spreadM),1,unfilterRcvd)));
title('Oversampled M Sequence Correlated with Unfiltered Rcvd')
% observe first max at x = 1044

% Begin downsampling by 4 from there
downsampleRcvd = unfilterRcvd(1044:4:end);

% Repeat m to match length of downsampled unfiltered Rcvd
repM = repmat(1-2*M,1,ceil(length(downsampleRcvd)/255));

% Descramble by applying m sequence
descrambled = downsampleRcvd.*repM(1:length(downsampleRcvd));

% Phase Shift (rotate to zero degrees)
theta = angle(descrambled);
rotated = exp(-1j*theta).*(descrambled);
scatterplot(rotated)
title('Rotated Message')

% Reshape to apply hadamard matrix
frames = (floor(length(rotated)/255)); % Finds number of complete frames
a = rotated(1:frames*255); % Take only complete frames
b = reshape(a,255,[]); % 1 Column = 1 Frame = 255 Chips
c = reshape(b(1:192,:),[],8); % Data is stored only in the first 192 chips
msg = reshape(c,8,[]); % 1 Column = 8 Chips

% Decode with 8-ary hadamard matrix
h = hadamard(8);
decoded = msg.'*h;

% BPSK Demod
demod = pskdemod(decoded,2);

% Take the output from Walsh Channel 5 
channel_5 =(demod(:,6)).';

% Convert binary to decimal to char
c_out = zeros(1,length(channel_5)/8); 
for i = 1:length(channel_5)/8
    c_out(i) = (bi2de(channel_5((i-1)*8+1:i*8),'right-msb'));
end

%Print Secret Message
SecretMessage = char(c_out) 

%'I've got a theory, it could be bunnies...þ'
% Buffy the Vampite Slayer

%% Functions %%
function [M] = mSeqGen(G)
    lfsr = zeros(1,length(G));
    lfsr(end) = 1;

    M = zeros(1,255);
    M(end) = lfsr(end);

    for m = 1:length(M)-1
        last = lfsr(end);
        for i = length(lfsr):-1:2
            lfsr(i) = mod(lfsr(i-1)+G(i)*last,2);
        end
        lfsr(1)=last;
        M(length(M)-m) = lfsr(end);
    end
end
