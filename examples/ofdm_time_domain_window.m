clc;
clear;
close all;

% 定义时域样本点数量和采样间隔
N = 1024; % 样本数量
T = 1/1000; % 采样间隔，对应1000Hz的采样频率

width = 10; % 窗口边沿宽度

%% 创建时域矩形信号 (boxcar function)
t = (-N/2:N/2-1)*T; % 时域坐标（考虑到对称性和0频分量位置）
rect_width = 100*T; % 矩形信号的宽度，例如100个采样间隔
rect_win = double(abs(t) < rect_width/2); % 生成矩形波形
rect_win = circshift(rect_win, 50);
% plot(t, rect_win);

% 第二个符号
rect_win2 = circshift(rect_win, 100);
% plot(t, rect_win2);

% 第三个符号
rect_win3 = circshift(rect_win2, 100);

%% 创建时域升余弦信号
rcos_win_interval = (0:(1/width):(1-1/width))*pi;
rcos_win1 = 0.5 - 0.5*cos(rcos_win_interval);
% plot(rcos_win1);
rcos_win2 = 0.5 + 0.5*cos(rcos_win_interval);
% plot(rcos_win2);
rcos_win = [zeros(1,462-width/2), rcos_win1, ones(1,100-width), rcos_win2, zeros(1,462-width/2)];
rcos_win = circshift(rcos_win, 50);
% plot(t, rcos_win);
% xlim([-0.1, 0.15])
% ylim([-2, 2]);

% 第二个符号
rcos_win2 = circshift(rcos_win, 100);
% plot(t, rcos_win2);
% 第三个符号
rcos_win3 = circshift(rcos_win2, 100);
% plot(t, rcos_win3);

%% 对信号加矩形窗
f = 10;
rect_sin_signal11 = sin(t*f*2*pi).*rect_win;
rect_sin_signal12 = sin(t*2*f*2*pi+pi/4).*rect_win;
rect_sin_signal13 = sin(t*3*f*2*pi+3*pi/4).*rect_win;
figure(1);
plot(t, rect_sin_signal11, 'g-', t, rect_sin_signal12, 'b-', t, rect_sin_signal13, 'r-');
hold on;
rect_sin_signal21 = sin(t*f*2*pi+pi/2).*rect_win2;
rect_sin_signal22 = sin(t*2*f*2*pi-pi/4).*rect_win2;
rect_sin_signal23 = sin(t*3*f*2*pi+pi).*rect_win2;
plot(t, rect_sin_signal21, 'g-', t, rect_sin_signal22, 'b-', t, rect_sin_signal23, 'r-');
hold on;
rect_sin_signal31 = sin(t*f*2*pi-pi/2).*rect_win3;
rect_sin_signal32 = sin(t*2*f*2*pi+pi/4).*rect_win3;
rect_sin_signal33 = -sin(t*3*f*2*pi-pi).*rect_win3;
plot(t, rect_sin_signal31, 'g-', t, rect_sin_signal32, 'b-', t, rect_sin_signal33, 'r-');

xlim([-0.01, 0.31])
ylim([-1.2, 1.2]);

%% 画加了矩形窗后的信号时域图
rect_sin_signal1 = rect_sin_signal11 + rect_sin_signal12 + rect_sin_signal13;
rect_sin_signal2 = rect_sin_signal21 + rect_sin_signal22 + rect_sin_signal23;
rect_sin_signal3 = rect_sin_signal31 + rect_sin_signal32 + rect_sin_signal33;

figure(2);
plot(t, rect_sin_signal1);
hold on;
plot(t, rect_sin_signal2);
hold on;
plot(t, rect_sin_signal3);
xlim([-0.01, 0.31])
ylim([-3, 3]);


%% 对信号加升余弦窗
f = 10;
rcos_sin_signal11 = sin(t*f*2*pi).*rcos_win;
rcos_sin_signal12 = sin(t*2*f*2*pi+pi/4).*rcos_win;
rcos_sin_signal13 = sin(t*3*f*2*pi+3*pi/4).*rcos_win;
figure(3);
plot(t, rcos_sin_signal11, 'g-', t, rcos_sin_signal12, 'b-', t, rcos_sin_signal13, 'r-');
hold on;
rcos_sin_signal21 = sin(t*f*2*pi+pi/2).*rcos_win2;
rcos_sin_signal22 = sin(t*2*f*2*pi-pi/4).*rcos_win2;
rcos_sin_signal23 = sin(t*3*f*2*pi+pi).*rcos_win2;
plot(t, rcos_sin_signal21, 'g-', t, rcos_sin_signal22, 'b-', t, rcos_sin_signal23, 'r-');
hold on;
rcos_sin_signal31 = sin(t*f*2*pi-pi/2).*rcos_win3;
rcos_sin_signal32 = sin(t*2*f*2*pi+pi/4).*rcos_win3;
rcos_sin_signal33 = -sin(t*3*f*2*pi-pi).*rcos_win3;
plot(t, rcos_sin_signal31, 'g-', t, rcos_sin_signal32, 'b-', t, rcos_sin_signal33, 'r-');

hold on;
plot(t, rcos_win);
hold on;
plot(t, rcos_win2);
hold on;
plot(t, rcos_win3);

xlim([-0.01, 0.31])
ylim([-1.2, 1.2]);

%% 画加了升余弦窗后的信号时域图
rcos_sin_signal1 = rcos_sin_signal11 + rcos_sin_signal12 + rcos_sin_signal13;
rcos_sin_signal2 = rcos_sin_signal21 + rcos_sin_signal22 + rcos_sin_signal23;
rcos_sin_signal3 = rcos_sin_signal31 + rcos_sin_signal32 + rcos_sin_signal33;

figure(4);
plot(t, rcos_sin_signal1);
hold on;
plot(t, rcos_sin_signal2);
hold on;
plot(t, rcos_sin_signal3);
xlim([-0.01, 0.31])
ylim([-3, 3]);

%% 使用FFT将时域信号转换到频域
%% 矩形窗
rect_sin_signal = rect_sin_signal1 + rect_sin_signal2 + rect_sin_signal3;

figure(5);
subplot(4,1,1);
plot(t, rect_sin_signal);
xlim([-0.01, 0.31]);
ylim([-3, 3]);

% 第一个符号
start_point = 513;
rect_sin_sig = rect_sin_signal1(start_point:start_point+100-1);
X_rect_sin = fft(rect_sin_sig); % 进行FFT
X_rect_sin_shifted = fftshift(X_rect_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

tmp = max(abs(X_rect_sin_shifted));
X = X_rect_sin_shifted / tmp;  % 取模并正规化至0至1

Xf = (-100/2:100/2-1)*(1/(T*100)); % 频率坐标
subplot(4,1,2);
stem(Xf, abs(X), '.');
disp('矩形窗第一个符号')
disp(X(51:60));

% 第二个符号
start_point = 513 + 100;
rect_sin_sig = rect_sin_signal2(start_point:start_point+100-1);
X_rect_sin = fft(rect_sin_sig); % 进行FFT
X_rect_sin_shifted = fftshift(X_rect_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

tmp = max(abs(X_rect_sin_shifted));
X = X_rect_sin_shifted / tmp;  % 取模并正规化至0至1

Xf = (-100/2:100/2-1)*(1/(T*100)); % 频率坐标
subplot(4,1,3);
stem(Xf, abs(X), '.');
disp('矩形窗第二个符号')
disp(X(51:60));

% 第三个符号
start_point = 513 + 100 + 100;
rect_sin_sig = rect_sin_signal3(start_point:start_point+100-1);
X_rect_sin = fft(rect_sin_sig); % 进行FFT
X_rect_sin_shifted = fftshift(X_rect_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

tmp = max(abs(X_rect_sin_shifted));
X = X_rect_sin_shifted / tmp;  % 取模并正规化至0至1

Xf = (-100/2:100/2-1)*(1/(T*100)); % 频率坐标
subplot(4,1,4);
stem(Xf, abs(X), '.');
disp('矩形窗第三个符号')
disp(X(51:60));

%% 升余弦窗
rcos_sin_signal = rcos_sin_signal1 + rcos_sin_signal2 + rcos_sin_signal3;

figure(6);
subplot(4,1,1);
plot(t, rcos_sin_signal);
xlim([-0.01, 0.31])
ylim([-3, 3]);

% 第一个符号
start_point = 513;
rcos_sin_sig = rcos_sin_signal1(start_point:start_point+100-1);
X_rcos_sin = fft(rcos_sin_sig); % 进行FFT
X_rcos_sin_shifted = fftshift(X_rcos_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

tmp = max(abs(X_rcos_sin_shifted));
X = X_rcos_sin_shifted / tmp;  % 取模并正规化至0至1

Xf = (-100/2:100/2-1)*(1/(T*100)); % 频率坐标
subplot(4,1,2);
stem(Xf, abs(X), '.');
disp('升余弦窗第一个符号')
disp(X(51:60));

% 第二个符号
start_point = 513 + 100;
rcos_sin_sig = rcos_sin_signal2(start_point:start_point+100-1);
X_rcos_sin = fft(rcos_sin_sig); % 进行FFT
X_rcos_sin_shifted = fftshift(X_rcos_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

tmp = max(abs(X_rcos_sin_shifted));
X = X_rcos_sin_shifted / tmp;  % 取模并正规化至0至1

Xf = (-100/2:100/2-1)*(1/(T*100)); % 频率坐标
subplot(4,1,3);
stem(Xf, abs(X), '.');
disp('升余弦窗第二个符号')
disp(X(51:60));

% 第三个符号
start_point = 513 + 100 + 100;
rcos_sin_sig = rcos_sin_signal3(start_point:start_point+100-1);
X_rcos_sin = fft(rcos_sin_sig); % 进行FFT
X_rcos_sin_shifted = fftshift(X_rcos_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

tmp = max(abs(X_rcos_sin_shifted));
X = X_rcos_sin_shifted / tmp;  % 取模并正规化至0至1

Xf = (-100/2:100/2-1)*(1/(T*100)); % 频率坐标
subplot(4,1,4);
stem(Xf, abs(X), '.');
disp('升余弦窗第三个符号')
disp(X(51:60));


%% 功率分析
%% 矩形窗
Y_rect_sin = fft(rect_sin_signal); % 进行FFT
Y_rect_sin_shifted = fftshift(Y_rect_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

% 创建频域坐标系
Yf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标

% 绘制时域信号
figure(7);
subplot(2,1,1);
plot(t, rect_sin_signal);
title('Time Domain Signal - Rectangle OFDM ');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-0.1, 0.4]);

% 绘制频域信号
subplot(2,1,2);
plot(Yf, 10*log10(abs(Y_rect_sin_shifted).^2/max(abs(Y_rect_sin_shifted))^2)); % 取模并正规化至0至1
title('Power/Frequency(dB/Hz) - Rectangle OFDM');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
ylim([-140,0]);

%% 升余弦窗
Y_rcos_sin = fft(rcos_sin_signal); % 进行FFT
Y_rcos_sin_shifted = fftshift(Y_rcos_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

% 创建频域坐标系
Yf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
figure(8);
subplot(2,1,1);
plot(t, rcos_sin_signal);
title('Time Domain Signal - Rcos OFDM ');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-0.1, 0.4]);

% 绘制频域信号
subplot(2,1,2);
plot(Yf, 10*log10(abs(Y_rcos_sin_shifted).^2/max(abs(Y_rcos_sin_shifted))^2)); % 取模并正规化至0至1
title('Power/Frequency(dB/Hz) - Rcos OFDM ');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
ylim([-140,0]);

% 用 pwelch 画功率谱对比
figure(9);
subplot(2,1,1);
pwelch(rect_sin_signal);
title('Power Spectral Density for rect OFDM');
subplot(2,1,2);
pwelch(rcos_sin_signal);
title('Power Spectral Density for rcos OFDM');

% 提示：由于离散傅立叶变换带来的影响，使用FFT产生的频率信号是对实际连续sinc函数的近似。