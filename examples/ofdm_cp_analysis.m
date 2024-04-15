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
rect_sin_signal11 = sin((t-0.02)*f*2*pi).*rect_win;
rect_sin_signal12 = sin((t-0.02)*1.5*f*2*pi).*rect_win;
rect_sin_signal13 = sin((t-0.02)*2*f*2*pi).*rect_win;

rect_sin_signal21 = -sin((t-0.02)*f*2*pi).*rect_win2;
rect_sin_signal22 = sin((t-0.02)*1.5*f*2*pi).*rect_win2;
rect_sin_signal23 = -sin((t-0.02)*2*f*2*pi).*rect_win2;

%% 不加 CP 时
figure(1);
subplot(4,1,1);
plot(t, rect_sin_signal11);
hold on;
plot(t, rect_sin_signal12);
hold on;
plot(t, rect_sin_signal13);
hold on;
plot(t, rect_sin_signal21);
hold on;
plot(t, rect_sin_signal22);
hold on;
plot(t, rect_sin_signal23);
xlim([0, 0.3]);
ylim([-1.2, 1.2]);

rect_sin_signal1 = rect_sin_signal11 + rect_sin_signal12 + rect_sin_signal13;
rect_sin_signal2 = rect_sin_signal21 + rect_sin_signal22 + rect_sin_signal23;
rect_sin_signal = rect_sin_signal1 + rect_sin_signal2;

subplot(4,1,2);
plot(t, rect_sin_signal);
xlim([0, 0.3]);

% 使用FFT将时域信号转换到频域
X_rect_sin = fft(rect_sin_signal); % 进行FFT
X_rect_sin_shifted = fftshift(X_rect_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

% 绘制频域信号
subplot(4,1,3);
Xf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
plot(Xf, abs(X_rect_sin_shifted) / max(abs(X_rect_sin_shifted))); % 取模并正规化至0至1

title('Frequency - No CP');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
xlim([-100, 100]);

% 绘制功率谱密度
subplot(4,1,4);
Xf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
plot(Xf, (abs(X_rect_sin_shifted).^2 / max(abs(X_rect_sin_shifted))^2)); % 取模并正规化至0至1
xlim([-100, 100]);

%% 加 CP 时
add_cp_signal11 = rect_sin_signal11;
add_cp_signal11(513:532) = add_cp_signal11(613:632);
add_cp_signal12 = rect_sin_signal12;
add_cp_signal12(513:532) = add_cp_signal12(613:632);
add_cp_signal13 = rect_sin_signal13;
add_cp_signal13(513:532) = add_cp_signal13(613:632);

add_cp_signal21 = rect_sin_signal21;
add_cp_signal21 = circshift(add_cp_signal21, 20);
add_cp_signal21(633:652) = add_cp_signal21(733:752);
add_cp_signal22 = rect_sin_signal22;
add_cp_signal22 = circshift(add_cp_signal22, 20);
add_cp_signal22(633:652) = add_cp_signal22(733:752);
add_cp_signal23 = rect_sin_signal23;
add_cp_signal23 = circshift(add_cp_signal23, 20);
add_cp_signal23(633:652) = add_cp_signal23(733:752);

figure(2);
subplot(4,1,1);
plot(t, add_cp_signal11);
hold on;
plot(t, add_cp_signal12);
hold on;
plot(t, add_cp_signal13);
hold on;
plot(t, add_cp_signal21);
hold on;
plot(t, add_cp_signal22);
hold on;
plot(t, add_cp_signal23);
xlim([-0.05, 0.3]);
ylim([-1.2, 1.2]);

add_cp_signal1 = add_cp_signal11 + add_cp_signal12 + add_cp_signal13;
add_cp_signal2 = add_cp_signal21 + add_cp_signal22 + add_cp_signal23;

add_cp_signal = add_cp_signal1 + add_cp_signal2;

subplot(4,1,2);
plot(t, add_cp_signal);
xlim([0, 0.3]);

% 使用FFT将时域信号转换到频域
X_cp = fft(add_cp_signal); % 进行FFT
X_cp_shifted = fftshift(X_cp); % 频域信号以0Hz为中心（将直流分量移到中央）

% 绘制频域信号
subplot(4,1,3);
Xf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
plot(Xf, abs(X_cp_shifted) / max(abs(X_cp_shifted))); % 取模并正规化至0至1
title('Frequency - with CP');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
xlim([-100, 100]);

subplot(4,1,4);
Xf = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标
plot(Xf, abs(X_cp_shifted).^2 / max(abs(X_cp_shifted))^2); % 取模并正规化至0至1
xlim([-100, 100]);


