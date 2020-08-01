clear
clc
%% init
j = sqrt(-1);
L = 128;
NTrain = 10000;
snrs=-10:2:20

Phase_low = 0;
Phase_high = pi/8;
NClass = 4;
train_data = zeros(NClass * NTrain, L);
train_label = zeros(NClass * NTrain, 2);
%% get train data
% bpsk
disp('bpsk begin');
signal_data = zeros(NTrain, L);
mod = comm.BPSKModulator();
x = constellation(mod);
xN = length(x);
for row = 1:NTrain
    EsNo = snrs(randperm(length(snrs),1)); % EsNo ~ U(EsNoLow, EsNoHigh)
    P = 10^(EsNo/10);
    train_label(row, 2) = EsNo;
    for col = 1:L
        s = x(unidrnd(xN));
        phase_jitter = 2*Phase_high*rand() - Phase_high; % phase jitter ~ U(-PhaseHigh, PhaseHigh)
        signal_data(row, col) = sqrt(P)*exp(j*phase_jitter)*s + sqrt(1/2)*(randn+j*randn);
    end
end
for row = 1:NTrain
    train_data(row, :) = signal_data(row, :);
    train_label(row, 1) = 0;
end
% pam4
disp('pam4 begin');
signal_data = zeros(NTrain, L);
noise_data = zeros(NTrain, L);
mod = comm.PAMModulator('ModulationOrder',4,'NormalizationMethod',...
        'Average power','AveragePower',1);
x = constellation(mod);
xN = length(x);
for row = 1:NTrain
    EsNo = snrs(randperm(length(snrs),1));
    P = 10^(EsNo/10);
    train_label(row+NTrain, 2) = EsNo;
    for col = 1:L
        s = x(unidrnd(xN));
        phase_jitter = 2*Phase_high*rand() - Phase_high;
        signal_data(row, col) = sqrt(P)*exp(j*phase_jitter)*s + sqrt(1/2)*(randn+j*randn);
    end
end
for row = 1:NTrain
    train_data(row + NTrain, :) = signal_data(row, :);
    train_label(row + NTrain, 1) = 1;
end
% psk8
disp('8psk begin');
signal_data = zeros(NTrain, L);
mod = comm.PSKModulator('ModulationOrder',8,'PhaseOffset',0);
x = constellation(mod);
xN = length(x);
for row = 1:NTrain
    EsNo = snrs(randperm(length(snrs),1));
    P = 10^(EsNo/10);
    train_label(row+ 2*NTrain, 2) = EsNo;
    for col = 1:L
        s = x(unidrnd(xN));
        phase_jitter = 2*Phase_high*rand() - Phase_high;
        signal_data(row, col) = sqrt(P)*exp(j*phase_jitter)*s + sqrt(1/2)*(randn+j*randn);
    end
end
for row = 1:NTrain
    train_data(row + 2*NTrain, :) = signal_data(row, :);
    train_label(row + 2*NTrain, 1) = 2;
end
% qam4
disp('qam4 begin');
signal_data = zeros(NTrain, L);
mod = comm.RectangularQAMModulator('ModulationOrder',4,...
        'NormalizationMethod','Average power','AveragePower',1);
x = constellation(mod);
xN = length(x);
for row = 1:NTrain
    EsNo = snrs(randperm(length(snrs),1));
    P = 10^(EsNo/10);
    train_label(row+ 3*NTrain, 2) = EsNo;
    for col = 1:L
        s = x(unidrnd(xN));
        phase_jitter = 2*Phase_high*rand() - Phase_high;
        signal_data(row, col) = sqrt(P)*exp(j*phase_jitter)*s + sqrt(1/2)*(randn+j*randn);
    end
end
for row = 1:NTrain
    train_data(row + 3*NTrain, :) = signal_data(row, :);
    train_label(row + 3*NTrain, 1) = 3;
end

save('train_data_awgn_10k.mat', 'train_data', '-mat');
save('train_label_awgn_10k.mat', 'train_label', '-mat');
