%%%This program takes the results of knife edge measurements, which should
%%%have the form of an error function and turns them into a 1/e^2 beam
%%%size.

close all
clear all

fontsize=16;
shortcut1090 = 1; %1 if you want to use positions of only 90% and 10% maximum power to determine beam size (shortcut for beam profiling)

if shortcut1090==1
    %% Import data
    filepath = char('/Users/work/Documents/Nanofibers/Laser Profiling/Blue-detuned Laser Profiling/Telescope 1 Output Side/');
    inputfilename = char('ProbeLaser_NanofiberOutput_Telescope1_20130912_redcoord.dat');
    outputfilename = char('ProbeLaser_NanofiberOutput_Telescope1_20130912_redcoord.txt');
    filename = strcat(filepath,inputfilename);
    [measurementposition,data90,data10]=textread(filename,'%f%f%f','commentstyle','matlab'); %read in data file, z in cm and beamsize in um
    
    %% Calculate waist size based on Siegman formula (see p. 94 Pascal's Lab Notebook 1)
    waistsize = abs(data90-data10)./1.28; %[um] 1/e^2 beam radius as a function of measurement position along length of beam
    
    %% Plot waist size versus distance from lens
    figure(1)
    plot(measurementposition(:),waistsize(:),'sr','MarkerFaceColor','r');
    hold on
    set(gca,'FontSize',fontsize,'FontWeight','bold');
    xlabel('Longitudinal Position [cm]','FontSize',fontsize,'FontWeight','bold');
    ylabel('Waist Size [mm]','FontSize',fontsize,'FontWeight','bold');
    %xlim([0 1.13]);
    ylim([0 max(waistsize(:))+0.1.*max(waistsize(:))]);
    %legend('Vertical Size','Horizontal Size','Location','best');
    %title('Beam Profile of 2988 Laser','FontSize',fontsize,'FontWeight','bold');
    
else
    measurementposition = [4.32 6.7]; %[cm] position relative to beam's zero point
    numfiles = length(measurementposition); %there should be one data file per measurement position
    %Initialize vectors before the for loop
    waistsize = zeros(1,length(measurementposition)); %[mm] initialize waist size vector
    peakpower = zeros(1,length(measurementposition)); %[W] initialize laser beam power vector
    for i=1:numfiles
        %% Load Data
        if i==1
            filename = char('/Users/work/Documents/Nanofibers/MOT Lasers/Beam Characterization/2988/BeamSize_Toptica_2988_DirectOutput_V.dat');
            oguess=-1;
        elseif i==2
            filename = char('/Users/work/Documents/Nanofibers/MOT Lasers/Beam Characterization/2988/BeamSize_Toptica_2988_DirectOutput_H.dat');
            oguess=1;
        elseif i==3
            filename = char('/Users/work/Documents/Spin Squeezing/Equipment/Fiber Lasers/Calibrations/Beam Size/PM Measurement/100percent_30cm.dat');
        elseif i==4
            filename = char('/Users/work/Documents/Spin Squeezing/Equipment/Fiber Lasers/Calibrations/Beam Size/PM Measurement/100percent_40p5cm.dat');
        elseif i==5
            filename = char('/Users/work/Documents/Spin Squeezing/Equipment/Fiber Lasers/Calibrations/Beam Size/AM Measurement/7Whighpowercoll47.dat');
        end
        [xdata ydata]=textread(filename,'%f%f','commentstyle','matlab'); %read in data file, z in cm and beamsize in mm
        xdata=xdata;
        %% Fitting
        %Guesses for data
        Aguess = (max(ydata)-min(ydata))/2;
        %oguess = mean(xdata)/(max(ydata)-min(ydata));
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
        xlabel('Knife Edge Position [mm]','FontSize',fontsize,'FontWeight','bold');
        ylabel('Measured Power [W]','FontSize',fontsize,'FontWeight','bold');
        text(mean(xdata)+0.2,mean(ydata),strcat('1/e^2 radius = ',num2str(abs(2*sigma)),' mm'),'FontSize',fontsize,'FontWeight','bold');
        %title(strcat('Beam Size at z = ',num2str(measurementposition(i)),' cm'),'FontSize',fontsize,'FontWeight','bold');
        
        waistsize(i) = abs(2*sigma); %[mm] 1/e^2 radius of laser beam in mm
        peakpower(i) = abs(A); %[W] power of laser beam in watts
    end
end %end of shortcut1090 option

outputdata = [measurementposition(:) waistsize(:)]; %[cm mm] data is output
outputfilepath = strcat(filepath,outputfilename);
dlmwrite(outputfilepath,outputdata,'\t');
%Done