%%%This program takes the results of knife edge measurements at 10% and 90%
%%%levels of the maximum power and turns them into 1/e^2 beam sizes.  This
%%%technique is a shortcut to full bema profiling.
close all
clear all

fontsize=16;

measurementposition = [9.5 17 30 40.5 119.3]; %[cm] position relative to beam's zero point
numfiles = 5; %length(measurementposition); %there should be one data file per measurement position

%Initialize vectors before the for loop
waistsize = zeros(1,length(measurementposition));
peakpower = zeros(1,length(measurementposition));
for i=1:numfiles
    %% Load Data
    if i==1
        filename = char('/Users/work/Documents/Spin Squeezing/Equipment/Fiber Lasers/Calibrations/Beam Size/PM Measurement/100percent_9p5cm.dat');
    elseif i==2
        filename = char('/Users/work/Documents/Spin Squeezing/Equipment/Fiber Lasers/Calibrations/Beam Size/PM Measurement/100percent_17cm.dat');
    elseif i==3
        filename = char('/Users/work/Documents/Spin Squeezing/Equipment/Fiber Lasers/Calibrations/Beam Size/PM Measurement/100percent_30cm.dat');
    elseif i==4
        filename = char('/Users/work/Documents/Spin Squeezing/Equipment/Fiber Lasers/Calibrations/Beam Size/PM Measurement/100percent_40p5cm.dat');
    elseif i==5
        filename = char('/Users/work/Documents/Spin Squeezing/Equipment/Fiber Lasers/Calibrations/Beam Size/AM Measurement/7Whighpowercoll47.dat');
    end
    [xdata ydata]=textread(filename,'%f%f','commentstyle','matlab'); %read in data file, z in cm and beamsize in mm
    
    %% Fitting
    %Guesses for data
    Aguess = (max(ydata)-min(ydata))/2;
    oguess = mean(xdata)/(max(ydata)-min(ydata));
    cguess = mean(xdata);
    InitialGuesses = [Aguess cguess oguess];
    
    [P,r,J] = nlinfit(xdata,ydata,@ErfFitFunction,InitialGuesses); %fit error function to knife edge data
    
    A=P(1); %amplitude (power)
    c=P(2); %position of center of beam
    sigma=P(3); %sigma is one standard deviation of gaussian beam (1/e^2 radius = 2*sigma)
    
    %Create vectors for plotting fits
    x=min(xdata):.1:max(xdata);
    erffit = (A./2).*(1-erf((x-c)./(sqrt(2).*sigma))); %vector of error function fit to knife edge measurements
    gaussianfromfit = (A/2).*exp(-(x-c).^2./(2.*sigma.^2)); %vector of gaussian shape of cloud as function of position
 
    %% Plotting
    figure(i)
    scatter(xdata,ydata,'or','MarkerFaceColor','r') %plot data
    hold on
    plot(x,erffit,':b');
    plot(x,gaussianfromfit,'-b');
    
    set(gca,'FontSize',fontsize,'FontWeight','bold');
    xlabel('Knife Edge Position [cm]','FontSize',fontsize,'FontWeight','bold');
    ylabel('Measured Power [W]','FontSize',fontsize,'FontWeight','bold');
    text(mean(xdata),2,strcat('1/e^2 radius = ',num2str(abs(2*sigma)),' mm'),'FontSize',fontsize,'FontWeight','bold');
    title(strcat('Beam Size at z = ',num2str(measurementposition(i)),' cm'),'FontSize',fontsize,'FontWeight','bold');

    waistsize(i) = abs(2*sigma);
    peakpower(i) = abs(A);
end

outputdata = [measurementposition(:) waistsize(:)]; %[cm mm] data is output
dlmwrite('/Users/work/Documents/Spin Squeezing/Equipment/Fiber Lasers/Calibrations/Beam Size/Beam_Profile_IPGLaserOutput_20110422.txt',outputdata,'\t');
%Done