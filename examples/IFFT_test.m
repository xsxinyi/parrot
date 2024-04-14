clc;
clear;
close all;

X = [0 0 2];
xn = ifft(X, 1024);
plot(real(xn));
