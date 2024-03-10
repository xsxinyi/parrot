%% 用窗函数法设计线性相位高通FIRDF，要求通带截止频率wp = pi/2，阻带截止频率ws = pi/4，通带最大衰减ap=1dB，阻带最小衰减as=40dB。
%% as = 40dB，可选用汉宁窗。过渡带宽度 DB = wp - ws = pi/4，汉宁窗的精确过渡带宽度 DB_ = 6.2pi/N。所以 6.2pi/N <= pi/4, N >= 24.8。

clear global;
close all;
wp = pi/2;
ws = pi/4;
DB = wp-ws;

N0 = ceil(6.2*pi/DB); %向上取整。
N = N0 + mod(N0+1, 2); %对于第一类相伴特性的FIRDF，高通滤波器N只能为奇数。
wc = (wp+ws)/2/pi;
hn = fir1(N-1, wc, 'high', hanning(N));

n = 0:N-1;

subplot(2,1,1);
stem(n, hn, 'filled');
axis([0,30,-1,1]);  %画针状图
title('(a)');xlabel('n');ylabel('hn(n)');


[H,w1]=freqz(hn,1);
% plot(w1/pi,20*log10(abs(H))); 
subplot(2,1,2);
plot(w1/pi,20*log10(abs(H))); %创建 Y 中数据对 X 中对应值的二维线图
axis([0,1,-80,10]); 
xlabel('归一化频率/p') ;
ylabel('幅度/dB') ;

