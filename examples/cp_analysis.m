clc;
clear;
close all;

% 定义时域样本点数量和采样间隔
N = 1024; % 样本数量
T = 1/1000; % 采样间隔，对应1000Hz的采样频率


%% 创建时域矩形信号 (boxcar function)
t = (-N/2:N/2-1)*T; % 时域坐标（考虑到对称性和0频分量位置）
rect_width = 100*T; % 矩形信号的宽度，例如100个采样间隔
rect_win = double(abs(t) < rect_width/2); % 生成矩形波形
% use 20 dot for cp.
rect_win = circshift(rect_win, 70);
% plot(t, rect_win);
rect_win2 = circshift(rect_win, 100);

% 对 sin 函数加矩形窗
f = 20;
rect_sin_signal1 = sin((t-0.02)*f*2*pi).*rect_win;
rect_sin_signal2 = sin((t-0.02)*1.5*f*2*pi).*rect_win;

%% 不加 CP 时
figure(1);
subplot(3,1,1);
plot(t, rect_sin_signal1);
hold on;
plot(t, rect_sin_signal2);

% 使用FFT将时域信号转换到频域
X_rect_sin1 = fft(rect_sin_signal1); % 进行FFT
X_rect_sin_shifted1 = fftshift(X_rect_sin1); % 频域信号以0Hz为中心（将直流分量移到中央）
X_rect_sin2 = fft(rect_sin_signal2); % 进行FFT
X_rect_sin_shifted2 = fftshift(X_rect_sin2); % 频域信号以0Hz为中心（将直流分量移到中央）

% 绘制频域信号
subplot(3,1,2);
Xf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
plot(Xf, abs(X_rect_sin_shifted1) / max(abs(X_rect_sin_shifted1))); % 取模并正规化至0至1
hold on;
plot(Xf, abs(X_rect_sin_shifted2) / max(abs(X_rect_sin_shifted2))); % 取模并正规化至0至1

title('Frequency - No CP');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
xlim([-100, 100]);

subplot(3,1,3);
rect_sin_signal = rect_sin_signal1 + rect_sin_signal2;
X_rect_sin = fft(rect_sin_signal); % 进行FFT
X_rect_sin_shifted = fftshift(X_rect_sin); % 频域信号以0Hz为中心（将直流分量移到中央）
Xf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
plot(Xf, abs(X_rect_sin_shifted) / max(abs(X_rect_sin_shifted))); % 取模并正规化至0至1
xlim([-100, 100]);

%% 加 CP 时
add_cp_signal1 = rect_sin_signal1;
add_cp_signal1(514:533) = rect_sin_signal1(613:632);
add_cp_signal2 = rect_sin_signal2;
add_cp_signal2(514:533) = rect_sin_signal2(613:632);

figure(2);
subplot(3,1,1);
plot(t, add_cp_signal1);
hold on;
plot(t, add_cp_signal2);

% 使用FFT将时域信号转换到频域
X_cp1 = fft(add_cp_signal1); % 进行FFT
X_cp_shifted1 = fftshift(X_cp1); % 频域信号以0Hz为中心（将直流分量移到中央）
X_cp2 = fft(add_cp_signal2); % 进行FFT
X_cp_shifted2 = fftshift(X_cp2); % 频域信号以0Hz为中心（将直流分量移到中央）

% 绘制频域信号
subplot(3,1,2);
Xf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
plot(Xf, abs(X_cp_shifted1) / max(abs(X_cp_shifted1))); % 取模并正规化至0至1
hold on;
Xf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
plot(Xf, abs(X_cp_shifted2) / max(abs(X_cp_shifted2))); % 取模并正规化至0至1
title('Frequency - with CP');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
xlim([-100, 100]);

subplot(3,1,3);
add_cp_signal = add_cp_signal1 + add_cp_signal2;
X_cp = fft(add_cp_signal); % 进行FFT
X_cp_shifted = fftshift(X_cp); % 频域信号以0Hz为中心（将直流分量移到中央）
Xf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
plot(Xf, abs(X_cp_shifted) / max(abs(X_cp_shifted))); % 取模并正规化至0至1
xlim([-100, 100]);
