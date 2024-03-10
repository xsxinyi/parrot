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

subplot(3,1,1);
stem(n, hn, 'filled');
axis([0,30,-1,1]);  %画针状图
title('(a)h(n)波形');xlabel('n');ylabel('h(n)');


[H,w1]=freqz(hn,1, 'whole'); % 'whole' 将区间设置为0-2pi, 否则是 0-pi。
% plot(w1/pi,20*log10(abs(H))); 
subplot(3,1,2);
plot(w1/pi,20*log10(abs(H))); %创建 Y 中数据对 X 中对应值的二维线图
axis([0,2,-80,10]); 
title('(b)损耗函数曲线');
xlabel('归一化频率/p') ;
ylabel('幅度/dB') ;

subplot(3,1,3);
plot(w1/pi,angle(H)); %创建 Y 中数据对 X 中对应值的二维线图
axis([0,2,-5,5]); 
title('(c)相频响应曲线');
xlabel('\omega/\pi') ;
ylabel('\phi(\omega)') ;

