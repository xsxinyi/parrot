%%%%%%%%%程序说明
%%% b = rcosdesign(beta,span,sps,shape)
%%% beta：滚降因子
%%% span: 表示截断的符号范围，对滤波器取了几个Ts的长度
%%% sps:每个符号的采样数
%%% shape:可选择'normal'或者'sqrt'
%%% b:1*（sps*span）的行向量，升余弦或余弦滤波器的系数

%****************************  程序主体 ****************%
clear global;
hn = rcosdesign(0.8,8,16,'sqrt');
fvtool(hn,'Analysis','impulse'); %将脉冲响应可视化

[H, w] = freqz(hn,1);

% figure(1);
% plot(hn);
% 
% figure(2);
% plot(w/pi, 20*log10(abs(H)));
% 
% figure(3);
% plot(w/pi, angle(H));

%****************************  仿真结论 ***************%
%%%%%成型滤波器rcosdesign的使用
%%%新增freqz使用