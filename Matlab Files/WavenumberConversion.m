%% Wavelength to Wavenumber Converter
%%% I'm going crazy trying to convert cm^-1 to wavelength, and viceversa.

close all
clear all

%% Convert Wavelength to Wavenumber
lambdaRange = [680e-9];          %[m]
waveNum = 1./lambdaRange;               %[m^-1]
waveNumcm = waveNum./100                %[cm^-1]

%% Convert Linewidth in Wavelength to Resolution in Wavenumber
% deltaLambda = [0.01e-9];                %[nm] linewidth
% LambdaCenter = [1060e-9];               %[m] center wavelength
% waveNumCenter = 1./LambdaCenter;        %[m^-1]
% waveNumCenter = waveNumCenter./100      %[cm^-1]
% deltaWaveNum = deltaLambda./LambdaCenter.^2

%% Convert Wavenumber to Wavelength
waveNuminversecm = [10000 20000]; %[cm^-1]
lambda = 1./(100.*waveNuminversecm); %[m]
lambdanm = lambda.*1e9 %[nm]