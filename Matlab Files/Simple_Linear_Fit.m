%% Overview
%%% This script will generate a linear fit to a set of imported x and y
%%% data.  Typically, I will use it to fit photodiode calibrations and the
%%% like.

close all
clear all

fontsize=16;

%% Import data
filepath = char('/Users/work/Documents/Spin Squeezing/Equipment/Thorlabs PDA100A Photodiode/20110927_VoltagetoPower_Calibration.dat');
[voltagePD opticalpower]=textread(filepath,'%f%f','commentstyle','matlab'); %[Volts Watts]; read in data file, voltagePD in volts and optical power in watts

%% Fit the data
InitialGuess=[1 0.01]; %[Volts/Watts Watts] initial guess for waist value and position of waist
errorvector = opticalpower./opticalpower; %weight all points equally in absence of better information
[P,r,J]=nlinfitweight(voltagePD(1:5),opticalpower(1:5),@fitalinedangit,InitialGuess,errorvector(1:5));

fitslope = P(1) %[Watts/Volt] Slope of line that fits watts/voltage calibration
fitycross = P(2) %[Watts] Power read by photodiode even when there's no light

%% Generate vector that represents the fit to the points
stepsize = (max(voltagePD)-min(voltagePD))/100;
fitx = min(voltagePD):stepsize:max(voltagePD);
fitline = fitslope.*fitx + fitycross; %[Watts/Volt]

%% Plot the results
figure(1)
plot(voltagePD(:),opticalpower(:),'sr','MarkerFaceColor','r','MarkerSize',10)
hold on
plot(fitx(:),fitline(:),'-k','LineWidth',1.2)
xlim([0 0.3]);
set(gca,'FontSize',fontsize,'FontWeight','bold');
xlabel('Photodiode Voltage [Voltage]','FontSize',fontsize,'FontWeight','bold');
ylabel('Power Meter Reading [Watts]','FontSize',fontsize,'FontWeight','bold');
text(mean(voltagePD),mean(opticalpower)/2,strcat('P=',num2str(fitslope),'V+',num2str(fitycross)),'FontSize',fontsize,'FontWeight','bold');
title('2011.09.26 PDA100A Photodiode Calibration','FontSize',fontsize,'FontWeight','bold');