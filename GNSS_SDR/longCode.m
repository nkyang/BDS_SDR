function fftCode = longCode( settings )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
longCode = makeCaTable(settings,1,0);
longCodeDelay = longCode.*circshift(longCode,[0,22]);
longCodeDelay = double(longCodeDelay);
fftCode  = conj(fft(longCodeDelay,[],2));
fftCode = fftCode(settings.acqSatelliteList,:);