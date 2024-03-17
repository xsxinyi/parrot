clear all;
close all;

format long;
tic;

%****************************  程序主体 *************
%%%%%%%%%%%%%%%%%  参数设定   %%%%%%%%%%%%%%
fc = 2000; %载波频率

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   信源    %%%%%%%%%%%%%%%%%%%
msg_source = [1 0 1 0 1 0 1 0 1 0 1 0 ];

figure(1);
subplot(3,1,1);
n = 0:length(msg_source)-1;
stem(n, msg_source, 'filled');
axis([0,12,-1,2]);  %画针状图
title('(a)信源时域波形');xlabel('n');ylabel('x(n)');
subplot(3,1,2);
Xk = fft(msg_source);
wn = 0:length(Xk)-1;
stem(2*wn./length(Xk), abs(Xk), 'filled');
axis([0,2,-1,10]);  %画针状图
title('(b)12点FFT得到的信源频域波形');xlabel('\omega/\pi');ylabel('X(k)');
subplot(3,1,3);
[Xk2,w]=freqz(msg_source, 1, 'whole'); % 'whole' 将区间设置为0-2pi, 否则是 0-pi。
plot(w/pi,abs(Xk2)); %创建 Y 中数据对 X 中对应值的二维线图
title('(c)freqz得到的信源频域波形');xlabel('\omega/\pi');ylabel('X(k)');
















