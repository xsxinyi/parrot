% 定义时域样本点数量和采样间隔
N = 1024; % 样本数量
T = 1/1000; % 采样间隔，对应1000Hz的采样频率

width = 10; % 窗口边沿宽度

% 创建时域矩形信号 (boxcar function)
t = (-N/2:N/2-1)*T; % 时域坐标（考虑到对称性和0频分量位置）
rect_width = 100*T; % 矩形信号的宽度，例如100个采样间隔
rect_win = double(abs(t) < rect_width/2); % 生成矩形波形

% 创建时域升余弦信号
rcos_win_interval = (0:(1/width):(1-1/width))*pi;
rcos_win1 = 0.5 - 0.5*cos(rcos_win_interval);
% plot(rcos_win1);
rcos_win2 = 0.5 + 0.5*cos(rcos_win_interval);
% plot(rcos_win2);
rcos_win = [zeros(1,462-width/2), rcos_win1, ones(1,100-width), rcos_win2, zeros(1,462-width/2)];

% 对 sin 函数分别加矩形窗和升余弦窗
f = 30;
rect_sin_signal = sin(t*f*2*pi).*rect_win;
rcos_sin_signal = sin(t*f*2*pi).*rcos_win;

% 使用FFT将时域信号转换到频域
X_rect_sin = fft(rect_sin_signal); % 进行FFT
X_rect_sin_shifted = fftshift(X_rect_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

X_rcos_sin = fft(rcos_sin_signal); % 进行FFT
X_rcos_sin_shifted = fftshift(X_rcos_sin); % 频域信号以0Hz为中心（将直流分量移到中央）

% 创建频域坐标系
f = (-N/2:N/2-1)*(1/(T*N)); % 频率坐标

% 绘制时域信号
figure(1);
subplot(2,1,1);
plot(t, rect_sin_signal);
title('Time Domain Signal - Rectangle * sin');
xlabel('Time (s)');
ylabel('Amplitude');

% 绘制频域信号
subplot(2,1,2);
plot(f, 10*log10(abs(X_rect_sin_shifted).^2/max(abs(X_rect_sin_shifted))^2)); % 取模并正规化至0至1
title('Power/Frequency(dB/Hz) - Rectangle * sin');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
ylim([-200,0]);

figure(2);
subplot(2,1,1);
plot(t, rcos_sin_signal);
title('Time Domain Signal - Rcos * sin');
xlabel('Time (s)');
ylabel('Amplitude');

% 绘制频域信号
subplot(2,1,2);
plot(f, 10*log10(abs(X_rcos_sin_shifted).^2/max(abs(X_rcos_sin_shifted))^2)); % 取模并正规化至0至1
title('Power/Frequency(dB/Hz) - Rcos * sin');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
ylim([-200,0]);

% 用 pwelch 画功率谱对比
figure(3);
subplot(2,1,1);
pwelch(rect_sin_signal);
title('Power Spectral Density for rect * sin');
subplot(2,1,2);
pwelch(rcos_sin_signal);
title('Power Spectral Density for rcos * sin');

% 提示：由于离散傅立叶变换带来的影响，使用FFT产生的频率信号是对实际连续sinc函数的近似。