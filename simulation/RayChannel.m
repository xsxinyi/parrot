%%%%%%%%%%%%%%%%%%%%%  瑞利 (Rayleigh) 信道仿真 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  RayChannel.m %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%程序说明
% 瑞利分布是一个均值为0，方差为σ²的平稳窄带高斯过程，其包络的一维分布是瑞利分布。
% 瑞利分布是最常见的用于描述平坦衰落信号接收包络或独立多径分量接受包络统计时变特性的一种分布类型。两个正交高斯噪声信号之和的包络服从瑞利分布。
% 瑞利衰落能有效描述存在能够大量散射无线电信号的障碍物的无线传播环境。若传播环境中存在足够多的散射，则冲激信号到达接收机后表现为大量统计独立的随机变量的叠加，
% 根据中心极限定理，则这一无线信道的冲激响应将是一个高斯过程。如果这一散射信道中不存在主要的信号分量，通常这一条件是指不存在直射信号（LoS），
% 则这一过程的均值为0，且相位服从0 到2π的均匀分布。即，信道响应的能量或包络服从瑞利分布。若信道中存在一主要分量，例如直射信号（LoS），则信道响应的包络服从莱斯分布，
% 对应的信道模型为莱斯衰落信道。(以上来自百度百科)

clc;
clear;
close all;
tic;

%% 瑞利信道
N = 100000;
level = 1000;

% randn 返回一个从标准正态分布中得到的随机标量
H = (randn(1, N) + 1j * randn(1, N)) / sqrt(2);
% H = ones(1,N);

figure(1); % 新建图形窗口
% 创建分布图
subplot(3,2,1);
% histogram(X) 基于 X 创建直方图。histogram 函数使用自动分 bin 算法，然后返回均匀宽度的 bin，这些 bin 可涵盖 X 中的元素范围并显示分布的基本形状。
% histogram 将 bin 显示为矩形条，这样每个矩形的高度就表示 bin 中的元素数量。
% histogram(X,nbins) 指定 bin 的数量。
H_hist = histogram(abs(H), level); % 绘制分布图
title('幅度分布图'); % 图形标题
xlabel('幅度大小'); % x轴标签
ylabel('分布'); % y轴标签
hold on;
ray_theory_x = 0:0.01:4;
% 瑞利分布的标准差会和实部或虚部的高斯分布的标准差一样，这里是 sqrt(2)
ray_theory_p = ray_theory_x.*exp(-ray_theory_x.^2/2/0.5) / 0.5;
ray_theory_p = (N*H_hist.BinWidth).*ray_theory_p;
plot(ray_theory_x, ray_theory_p, 'r');

subplot(3,2,2);
histogram(angle(H), level); % 绘制分布图
title('角度分布图'); % 图形标题
xlabel('角度大小'); % x轴标签
ylabel('分布'); % y轴标签

%% 信号经过瑞利信道
nsym = 4;
fs = 1;
sps = N/nsym;
x = [zeros(1, sps) ones(1, sps) zeros(1, sps) zeros(1, sps)];
x = 2*x - 1;
t = 0:nsym*sps-1;
t = t / fs / sps;
subplot(3,2,3);
plot(t, x);
ylim([-1.2, 1.2]);
title('输入信号时域图形'); % 图形标题
xlabel('时间（秒）'); % x轴标签
ylabel('幅度'); % y轴标签

X = fft(x);

% 给 X 作图
f_Nyquist = sps * fs;
% 1Hz 是几个采样点。
dot_per_f = length(X) / f_Nyquist;
% 方波的频率是 sinc 函数，只画到 20Hz 就够了。
f_plot = 20*dot_per_f;
X1 = X(1:f_plot);
subplot(3,2,4);
plot(abs(X1));
title('输入信号频域图形'); % 图形标题
xlabel('频率（Hz）'); % x轴标签
ylabel('幅度'); % y轴标签

Y = X.*H;

% 给 Y 作图
Y1 = Y(1:f_plot);
subplot(3,2,5);
plot(abs(Y1));
title('输出信号频域图形'); % 图形标题
xlabel('频率（Hz）'); % x轴标签
ylabel('幅度'); % y轴标签

y = ifft(Y);
subplot(3,2,6);
plot(t, y);
title('输出信号时域图形'); % 图形标题
xlabel('时间（秒）'); % x轴标签
ylabel('幅度'); % y轴标签

toc;


% figure(2);
% h = ifft(H);
% plot(abs(h));
