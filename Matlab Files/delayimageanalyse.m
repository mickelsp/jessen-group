%%%Duplicates delay_image_analyse from LabView: Created 2007.06.14
%%%Also adds fitting of lineshapes, which LabView version doesn't do.

%%%Import number, size, etc... from imagefit.  Then fit data to get
%%%temperature OR lifetime OR linewidth (Lorentzian lineshape).  See
%%%fittingtype option
%%%Outputs a ____TEMP.txt file OR a _____LIFE.txt file OR a
%%%____PEAKPARAMETERS.txt file.  See fittingtype option.

%%%Plot rms vs changed variable, number vs changed variable, peak OD vs
%%%changed variable.

%%%TO DO LIST:
%%% 1)use nlinfitweight for lifetime fitting and linewidth fitting.

close all
clear all
%clc

%% Initialization and Options
%%%Constants
lambda=461*10^(-9); %intercombination line wavelength
%lambda=689*10^(-9); %intercombination line wavelength
kBoltz=1.38*10^(-23);
mass=84*1.672*10^(-27); %strontium mass
templow=1*10^(-6);
h=6.626*10^(-34);
hbar=h/(2*pi);
w0=100*10^(-6);     %dipole trap beam waist
lambdaODT=1064*10^(-9); %dipole trap wavelength
factor=2*(h*150*10^3/kBoltz)*(50*10^(-6))^2; %To compute ODT trap depth, for 1064nm light; see MMA notebook on ODT; factor of 2 for crossed config
width1microkelvin=10^(-6)*sqrt(2*kBoltz*templow/mass)/lambda; %not used?
numberfactor=10^6; %scales number of atoms for fitting
timefactor=10^3; %scales time from ms to s
sizefactor=10^3; %scales size from m to mm5

%%%Options and initialization
fitting=0; %run fittin g routine?
readaveragedinput=0; %reads input from AveragingMFile.m (ODT Laser Cooling Studies); result of binning values together.
fittingtype=1; %1 = Temperature fit; 2 = Lifetime fit; 3 = Linewidth fit
sumfiles=0; %1 means we sum number of atoms for all files in run
scancenter=7.58;
if fittingtype==1
    xaxislabel='Time (ms)';
elseif fittingtype==2
    xaxislabel='Time (ms)';
elseif fittingtype==3
    xaxislabel='Binding Energy (MHz)';
    xaxislabel='DAC Voltage';
end
%%%Fit options for lifetime fits (fitting==2)
bodyloss = 1;   %bodyloss=1 for 1-body fits, 2 for 2-body fits
plot1body = 0;  %plot1body=0 plots the 1-body fits with data, 0 does not
density=0;  %density=0 to use number loss equation, 1 to use density loss equation
subsetsize=[18;18;5]; %specify what sizes subsets are; should be in DESCENDING order
plotsubset=0;  %1 means the SMALLEST subset is used as the lifetime fit function
optionsnlinfit=statset('TolFun',1e-9,'MaxIter',1000,'TolX',1e-9);
%%%Fit options for linewidth fits (fitting==3)
singlepeakfit=1; %1 to do the single lorentzian fitting; 0 to do 2-lorentzian fitting
convertfreq=0; %1 to perform voltage to frequency conversion (voltage to wavemeter reading)
wavemeteractive=0;  %1 wavemeter readings on; 2 wavemeter readings off, manually read beginning/end freqs; 0 load volt2freq data from file and fit
plotvoltfreq=0; %1 to plot voltage to frequency conversion fit; 0 fit, but not plot conversion
fitsubset=0; %1 takes the first numberoflines for fitting instead of all the data.
numberoflines=20; %determines how many points are used for fitting.
composite=1; %overlay data on single plot?
save2places=0;  %1 means we save output files to both primary and secondary directories for backup purposes
splitdelay=1; %In case we split the delay time between two columns, we can fix the delay time read in from the batch file;
psd=0;  %calculate psd, assumes you provide the appropriate effective volumes
pos=0; %plot cloud position (x and y).
%%%Plotting options
markersize=6;
fontsize=14; %font size for everything but text on plots
textsize=14; %font size for text on plots
tempconstant=1.4; %changes the vertical height where the temperature is displayed

%% Directory and Files
%%%File path information
%directory=char('\\magnesium.rice.edu\neutraldata\ODT\2009.12.17\');
%%Windows directory
%directory=char('\\Magnesium\neutraldata\Sr87\2010.01.19\'); %Windows directory
directory=char('/Users/Work/Documents/Analysis/Unforced Evap of 84Sr/');
basehead=char('Files_unforcedevap84Sr');

volt2freq=char(strcat(deblank('\\magnesium.rice.edu\neutraldata\Calibrations\'),'volt2freqZOS150-at107MHz-20080626.txt'));
basenamelist=char(strcat(deblank(directory),deblank(basehead),deblank('.txt')));

[directoryvector basenamevector timevector powerV splittime delaytimeBF OFRon timeoffset rfpower]=...
    textread(basenamelist, '%s%s%f%f%f%f%f%f%f','commentstyle','matlab'); %read in data file for atoms
numberofbasenames=size(basenamevector,1);
powerW=splittime;

% for i=1:numberofbasenames
%     if OFRon(i)==1
%         char(strcat('legendentry',num2str(i)))='OFR On';
%     else
%         char(strcat('legendentry',num2str(i)))='OFR Off';
%     end %of
% end %of for loop for legend creation

%%%Initialize vectors before looping
filecounter=0;   %this will keep track of all the files analyzed in all the batches

%% Loop Through The Basenames
for basenamenumber=1:numberofbasenames
    %% Load and Prep Data
    filecounter=filecounter+1
    
    %directory=basedirectory(basenamenumber);
    %directory=char('/Volumes/MAGNESIUM25/ODT/2007.12.18C/'); %Mac directory syntax
    clear basehead;
    basehead=char(basenamevector(basenamenumber))
    basenamenumber
    
    %OUTPUT FILES
    if fittingtype==1
        outputdatafile=char(strcat(deblank(directoryvector(filecounter)),deblank(basenamevector(filecounter)),deblank('TEMP.txt')));
        if save2places==1
            outputdatafile2=char(strcat(deblank(directory2),deblank(basenamevector(filecounter)),deblank('TEMP.txt')));
        end
    elseif fittingtype==2
        outputdatafile=char(strcat(deblank(directoryvector(filecounter)),deblank(basenamevector(filecounter)),deblank('LIFE.txt')));
        if save2places==1
            outputdatafile2=char(strcat(deblank(directory2),deblank(basenamevector(filecounter)),deblank('LIFE.txt')));
        end
    elseif fittingtype==3
        outputdatafile = char(strcat(deblank(directoryvector(filecounter)),deblank(basenamevector(filecounter)),'peakparameters.fit'));
        if save2places==1
            outputdatafile2 = char(strcat(deblank(directory2),deblank(basenamevector(filecounter)),'peakparameters.fit'));
        end
    end
    
    clear e1 xaxis e3 e4 e5 e6 e7 e8 e9 delaytime e11 e12
    clear f1 f2 f3 f4 f5 f6 positionx positiony f9 peakOD sigmax sigmay f13 f14 f15 f16 f17 f18 sigmaxerror sigmayerror f21 f22 f23 f24 f25 numatomserror numatoms chi3D sizeparameter fbarModel Tfermi
    %INPUT FILES
    if readaveragedinput==1 %reads input from AveragingMFile.m (ODT Laser Cooling Studies); result of binning values together.
        datafile=char(strcat(deblank(directoryvector(filecounter)),deblank(basehead),deblank('.txt')))
        [xaxis numatoms sigmax sigmay]=textread(datafile,'%f %f %f %f','commentstyle','matlab');
        numatoms = numatoms'./numberfactor; %number of atoms in millions of atoms
        xaxis=xaxis';
    else
        
        xaxisfile=char(strcat(deblank(directoryvector(filecounter)),deblank(basehead),deblank('.batch')))
        imagefile2D=char(strcat(deblank(directoryvector(filecounter)),deblank(basehead),deblank('2Dfitparams.txt')))
        %We use changed values from batch file and number of atoms from 2D fit parameters file (output of imagefit).
        [e1 xaxis e3 e4 ODTvoltage e6 e7 e8 e9 filedelaytime e11 e12]=textread(xaxisfile,...
            '%s%f%f%s%f%f%f%f%f%f%f%s','commentstyle','matlab');
        [f1 f2 f3 f4 f5 f6 positionx positiony f9 peakOD sigmax sigmay f13 f14 f15 f16 f17 f18 sigmaxerror sigmayerror f21 f22 f23 f24 f25 numatomserror numatoms chi3D sizeparameter]=textread(imagefile2D,...
            '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','commentstyle','matlab');
        if splitdelay==1
            %splittime=7.5; %for Nov. 2009 data, we split the delay time between two columns (12ms and 3 ms);
            filedelaytime=filedelaytime+splittime(filecounter);
        end
            
        imagesize=size(xaxis,1); %number of data points
        
        %%%Data preparation
        numatoms = numatoms(:)./numberfactor; %number of atoms in millions of atoms
        numsystematic=0.0; %systematic error in number of atoms
        %Add systematic and statistical error in quadrature to get total error in number
        numatomserror = (numsystematic.^2+numatomserror.^2).^(0.5)./numberfactor;
        xaxis=xaxis(:); delaytime=filedelaytime(imagesize)*ones(length(filedelaytime),1); %delaytime=delaytime(:);
        sigmax = abs(sigmax.*sizeparameter); %size is spit out of imagefit in # of pixels.
        sigmay = abs(sigmay.*sizeparameter); %sizeparameter converts to meters.
        sigmasystematic = 10.4*10^(-6); %one PixelFly pixel is our best resolution
        %Add systematic and statistical error in quadrature to get total error in sizes
        sigmaxerror = sqrt((sigmaxerror.*sizeparameter).^2+sigmasystematic.^2);
        sigmayerror = sqrt((sigmayerror.*sizeparameter).^2+sigmasystematic.^2);
        sigmainitial=sigmax(1); %initial size of cloud; using x instead of y is arbitrary; all for lifetime fitting
        samplesizeinitial=ones(1,length(sigmax)).*sigmainitial;   %vector containing initial size of cloud
        xsizetemp=(sigmax.^2.*mass)./(kBoltz.*(delaytime(imagesize)./timefactor).^2)./templow; %calculate temperature based on x size
        ysizetemp=(sigmay.^2.*mass)./(kBoltz.*(delaytime(imagesize)./timefactor).^2)./templow; %calculate temperature based on y size
        
        if psd==1
            %calculate phase space density for 2817 on 2009/07/22
            V1eff=[5.74*10^-16, 5.74*10^-16, 5.74*10^-16, 5.74*10^-16, 6.07*10^-16,6.77*10^-16, 7.43*10^-16,...
                8.05*10^-16, 8.66*10^-16, 9.3*10^-16,9.95*10^-16, 1.07*10^-15, 1.14*10^-15, 1.22*10^-15, 1.31*10^-15,...
                1.41*10^-15, 1.53*10^-15, 1.69*10^-15, 1.93*10^-15, 2.31*10^-15,2.92*10^-15, 3.88*10^-15, 5.38*10^-15,...
                7.59*10^-15, 1.08*10^-14];
            peakdensity=numatoms(:).*numberfactor./V1eff(:);
            debroglielambdax=(h.^2./(2.*pi.*mass.*kBoltz.*xsizetemp.*templow)).^(3/2);
            debroglielambday=(h.^2./(2.*pi.*mass.*kBoltz.*ysizetemp.*templow)).^(3/2);
            psdx=peakdensity.*debroglielambdax;
            psdy=peakdensity.*debroglielambday;
        end
        
        delaytime(1)
        %Get rid of any "bad" fits; negative fits
        counter=1;
        clear xaxismod nummod rmsymod peakODmod
        for i=1:length(xaxis)
            if numatoms(i)>-100000000 && numatoms(i)<500000000
                xaxismod(counter)=xaxis(i);
                nummod(counter)=numatoms(i);
                peakODmod(counter)=peakOD(i);
                rmsymod(counter)=sigmay(i);
                counter=counter+1;
            end
        end
        xaxis=xaxismod; numatoms=nummod; sigmay=rmsymod; peakOD=peakODmod;
    end %END if readaveragedinput==1 LOOP
    
    %% Fitting Section
    if convertfreq==1
        %Prepare data for fitting
        if wavemeteractive==1  %use if wavemeter was on
            wavemetererror=wavemeter./wavemeter; %need error to make nlinfit happy. this will weight all pts equally b/c we have no error estimate for now.
            InitialGuess = [(max(wavemeter)-min(wavemeter))./(min(xaxis)-max(xaxis)),max(wavemeter)]; %negative sloped guess b/c increased voltage reduces wavenumber
            % Calculate fit parameters
            [Pc,rc,Jc]=nlinfitweight(xaxis(:),wavemeter(:),@fitalinedangit,InitialGuess,wavemetererror(:));
            cic=nlparci(Pc,rc,Jc);
            slopeerror=(cic(1,2)-cic(1,1))./4;
            % Evaluate the fit
            pts=100;
            for t=1:pts
                indevariable(t) = min(xaxis)+(max(xaxis)-min(xaxis)).*(t./pts);
                fit(t) = Pc(2)+Pc(1).*indevariable(t); %absolute value b/c sometimes complex value is put in. not sure why...
            end
            if plotvoltfreq==1
                figure(100)
                plot(xaxis,wavemeter,'kx')
                hold on
                plot(indevariable,fit,'-b')
            end
            xaxis=Pc(2)+xaxis.*Pc(1);
            xaxis=xaxis.*30-9.9606*10^4; %convert from wavenumber to GHz
        elseif wavemeteractive==2 %use if wavemeter reading was manually read at beginning and end
            Pc(1) = (startfreq(filecounter)-stopfreq(filecounter))./(min(xaxis)-max(xaxis)); %negative sloped guess b/c increased voltage reduces wavenumber
            Pc(2) = startfreq(filecounter)-Pc(1).*xaxis(1);
            xaxis=Pc(2)+xaxis.*Pc(1);
            xaxis=xaxis.*30-9.9606*10^4; %convert from wavenumber to GHz
        else %if we load voltage to frequency conversion from file
            [volt freq]=textread(volt2freq,'%f%f','commentstyle','matlab');
            %volt=volt;
            freq=81.9-(2.*abs(freq));
            %Add 75.1 for 2008.02.07 PAS line data
            %subtract 75.1 for 2008.02.06 feature data (we're 110
            %MHZ blue) and make it negative.
            %add 75.1 for 2008.01.15 feature data (we're 110 MHz red)
            %add 106.5 MHz for 2007.12.18 data to get binding energy.
            %add 0.0 MHz (nothing) to 2008.01.05 data to get AOM#2 freq
            %add 98.5 to 2008.01.09 data to get binding energy
            conversionerror=freq./freq; %need error to make nlinfit happy. this will weight all pts equally b/c we have no error estimate for now.
            InitialGuess = [(max(freq)-min(freq))./(min(volt)-max(volt)),max(freq)]; %negative sloped guess b/c increased voltage reduces wavenumber
            % Calculate fit parameters
            [Pc,rc,Jc]=nlinfitweight(volt(:),freq(:),@fitalinedangit,InitialGuess,conversionerror(:));
            cic=nlparci(Pc,rc,Jc);
            slopeerror=(cic(1,2)-cic(1,1))./4;
            % Evaluate the fit
            pts=100;
            for t=1:pts
                freqfit(t) = min(volt)+(max(volt)-min(volt)).*(t./pts);
                fit(t) = Pc(2)+Pc(1).*freqfit(t); %absolute value b/c sometimes complex value is put in. not sure why...
            end
            if plotvoltfreq==1
                figure(100)
                plot(volt,freq,'kx')
                hold on
                plot(freqfit,fit,'-b')
                xlabel('Voltage','FontSize',fontsize,'FontWeight','bold');
                ylabel('Frequency [MHz]','FontSize',fontsize,'FontWeight','bold');
                title(char(volt2freq));
                set(gca,'FontSize',fontsize,'FontWeight','bold');
                %xlim([7 7.7]);
                hold off
            end
            xaxisfreq=Pc(2)+xaxis.*Pc(1);
        end %end of wavemeter active choice
    end %of voltage-to-frequency conversion
    
    if fitting==1
        %% Temperature Fitting
        if fittingtype==1 %fit sizes to extract temperatures
            %Prepare data for fitting
            xaxisfit = xaxis./timefactor; %convert hold time from ms to s
            xx = xaxisfit.*xaxisfit; sizex = sigmax.*sigmax; sizey = sigmay.*sigmay;
            sizexerror=sigmaxerror.*sigmaxerror; sizeyerror=sigmayerror.*sigmayerror;
            InitialGuessX = [(max(sizex)-min(sizex))./(max(xx)-min(xx)),0];
            InitialGuessY = [(max(sizey)-min(sizey))./(max(xx)-min(xx)),0];
            % Calculate fit parameters
            [px,r1,J1]=nlinfitweight(xx(:),sizex(:),@fitalinedangit, InitialGuessX,sizexerror(:));
            [py,r2,J2]=nlinfitweight(xx(:),sizey(:),@fitalinedangit, InitialGuessY,sizeyerror(:));
            cix=nlparci(px,r1,J1); ciy=nlparci(py,r2,J2);
            xslopeerror=(cix(1,2)-cix(1,1))./4; yslopeerror=(ciy(1,2)-ciy(1,1))./4;
            % Evaluate the fit
            pts=100;
            for t=1:pts
                indevariable(t) = min(xaxis)+(max(xaxis)-min(xaxis)).*(t./pts);
                fitx(t) = abs(sqrt(px(2)+px(1).*(indevariable(t)./timefactor).^2).*(sizefactor./timefactor)); %absolute value b/c sometimes complex value is put in. not sure why...
                fity(t) = abs(sqrt(py(2)+py(1).*(indevariable(t)./timefactor).^2).*(sizefactor./timefactor));
            end
            xtemperature = px(1).*mass.*1000./kBoltz;
            ytemperature = py(1).*mass.*1000./kBoltz;
            xtemperror = xslopeerror.*mass.*1000./kBoltz;
            ytemperror = yslopeerror.*mass.*1000./kBoltz;
            
            %Prepare to write to output file for temperature fitting
            Uodt(filecounter)=factor*powerW(filecounter)/(w0^2); %Peak trap depth for give power and waist of ODT beam
            Temperature(1)=xtemperature*10^3; %Temp in uK
            Temperature(2)=ytemperature*10^3; %Temp in uK;
            Temperature(3)=xtemperror*10^3; %Temp error in uK;
            Temperature(4)=ytemperror*10^3; %Temp error in uK;
            Temperature(5)=Uodt(filecounter)*10^6; %Peak trap depth in uK
            %Write all Temperature parameters to outputdatafile
            dlmwrite(outputdatafile, Temperature,'\t');
            if save2places==1
                dlmwrite(outputdatafile2, Temperature,'\t');
            end
            
        elseif fittingtype==2 %fit numbers to extract lifetimes
            %% Lifetime Fitting
            %Prepare data for fitting (one body decay rate only)
            %Initialguess2=[numatoms(1,1),0.00059]; %number, decay rate in s-1
            Initialguess2=[numatoms(1,1),0.1]; %number, decay rate in s-1
            data2=transpose([xaxis;numatoms]);
            TIME=[0:2:max(xaxis)];
            
            %Create subsets of the entire dataset and determine one-body
            %fits to subsets of data.
            xdatalength=size(data2,1);
            if xdatalength > subsetsize(1)
                dataLAST40 = transpose([xaxis(xdatalength-subsetsize(1):xdatalength);numatoms(xdatalength-subsetsize(1):xdatalength)]);
                dataLAST30 = transpose([xaxis(xdatalength-subsetsize(2):xdatalength);numatoms(xdatalength-subsetsize(2):xdatalength)]);
                dataLAST20 = transpose([xaxis(xdatalength-subsetsize(3):xdatalength);numatoms(xdatalength-subsetsize(3):xdatalength)]);
                [coeff40,r,J] = nlinfit(dataLAST40(:,1).*1e-3,dataLAST40(:,2),@Lifetimefitfunction,Initialguess2);
                fit40=coeff40(1).*exp(-coeff40(2)*dataLAST40(:,1));
                coeff40(2)
                N1body40=coeff40(1).*exp(-coeff40(2)*TIME*1e-3);
                [coeff30,r,J] = nlinfit(dataLAST30(:,1).*1e-3,dataLAST30(:,2),@Lifetimefitfunction,Initialguess2);
                fit30=coeff30(1).*exp(-coeff30(2)*dataLAST30(:,1));
                coeff30(2)
                N1body30=coeff30(1).*exp(-coeff30(2)*TIME*1e-3);
                [coeff20,r,J] = nlinfit(dataLAST20(:,1).*1e-3,dataLAST20(:,2),@Lifetimefitfunction,Initialguess2);
                fit20=coeff20(1).*exp(-coeff20(2)*dataLAST20(:,1));
                coeff20(2)
                N1body20=coeff20(1).*exp(-coeff20(2)*TIME*1e-3);
            elseif xdatalength  > subsetsize(2)
                dataLAST30 = transpose([xaxis(xdatalength-subsetsize(2):xdatalength);numatoms(xdatalength-subsetsize(2):xdatalength)]);
                dataLAST20 = transpose([xaxis(xdatalength-subsetsize(3):xdatalength);numatoms(xdatalength-subsetsize(3):xdatalength)]);
                [coeff30,r,J] = nlinfit(dataLAST30(:,1).*1e-3,dataLAST30(:,2),@Lifetimefitfunction,Initialguess2);
                fit30=coeff30(1).*exp(-coeff30(2)*dataLAST30(:,1));
                N1body30=coeff30(1).*exp(-coeff30(2)*TIME*1e-3);
                [coeff20,r,J] = nlinfit(dataLAST20(:,1).*1e-3,dataLAST20(:,2),@Lifetimefitfunction,Initialguess2);
                fit20=coeff20(1).*exp(-coeff20(2)*dataLAST20(:,1));
                N1body20=coeff20(1).*exp(-coeff20(2)*TIME*1e-3);
            elseif xdatalength > subsetsize(3)
                dataLAST20 = transpose([xaxis(xdatalength-subsetsize(3):xdatalength);numatoms(xdatalength-subsetsize(3):xdatalength)]);
                [coeff20,r,J] = nlinfit(dataLAST20(:,1).*1e-3,dataLAST20(:,2),@Lifetimefitfunction,Initialguess2);
                fit20=coeff20(1).*exp(-coeff20(2)*dataLAST20(:,1));
                N1body20=coeff20(1).*exp(-coeff20(2)*TIME*1e-3);
            end
            
            %One-body fits to entire dataset, not just subsets
            %[coeff2,r,J] = nlinfit(data2(:,1).*1e-3,data2(:,2),@Lifetimefitfunction,Initialguess2);
            [coeff2,r,J] = nlinfit(data2(:,1).*1e-3,data2(:,2),@Lifetimefitfunction,Initialguess2);
            ci2=nlparci(coeff2,r,J);
            fit2=coeff2(1).*exp(-coeff2(2).*data2(:,1));
            lifetime1body(filecounter)=1./coeff2(2);
            onebodydecay(filecounter)=coeff2(2);
            gammafitdev=(ci2(2,2)-ci2(2,1))./4;
            onebodydecayerror(filecounter)=gammafitdev;
            %lifetime1bodyerror =
            %sqrt((onebodydecayerror./onebodydecay.^2).^2);
            lifetime1bodyerror = onebodydecayerror./onebodydecay.^2;
            TIME=[0:2:max(xaxis)];
            N1body=coeff2(1).*exp(-coeff2(2)*TIME*1e-3);
            
            %%%Prepare data for fitting (two body decay rate)
            N0=ones(1,length(numatoms))*numatoms(1);  %in multiples of 10^6
            gamma=ones(1,length(numatoms))*coeff2(2); %determine gamma from one body fit to final x # pts.  then, use as constant for 2-body fits
            beta=ones(1,length(numatoms))*1.9733e-5;  %%%WHERE DID THIS NUMBER in the BETA formula come from???
            
            if bodyloss==2 %FOR 2-BODY LOSSES IN SAMPLE
                if density==1 %From density loss equation
                    Initialguess2body=[beta(1)]; %two-body loss rate; ONLY FITS BETA:  initial number and gamma are constants
                    data2body=transpose([xaxis;numatoms;samplesize;N0;gamma]); %time, #, initial size, fitted initial #, fitted 1-body decay rate
                    [coeff2body,r,J]=nlinfit([data2body(:,1),data2body(:,3), data2body(:,4), data2body(:,5)], data2body(:,2),@Lifetimefitfunction2bodydensityloss, Initialguess2body,optionsnlinfit);
                    factor2body=(sqrt(2*pi)*sigma)^3; %2*sqrt(2); %sigma is initial size of sample, defined above
                    TIME=data2body(:,1);
                    N2body=(coeff2(1).*exp(-coeff20(2).*TIME))./(1+(coeff2(1).*coeff2body(1)./(coeff20(2).*factor2body)).*(1-exp(-coeff20(2).*TIME)));
                    lifetime1body(numberofbasenames)=1/coeff20(2); onebodydecay(numberofbasenames)=coeff20(2); twobodydecay(numberofbasenames)=coeff2body(1);
                elseif density==0 %From number loss equation; FITS BETA, GAMMA, and INITIAL NUMBER of atoms
                    Initialguess2body=[numatoms(1),gamma(1),beta(1)]; %coeff20(2),(numatoms(1,1)/10^6) %initial # divided by 10^6, 1-body loss rate in s-1, 2-body loss rate in s-1
                    data2body=transpose([xaxis;numatoms;beta]);
                    [coeff2body,r,J]=nlinfitweight(data2body(:,1),data2body(:,2),@Lifetimefitfunction2body, Initialguess2body,numatomserror(:));
                    ci2=nlparci(coeff2body,r,J);
                    initialnumber=coeff2body(1);
                    initialnumberdev=(ci2(1,2)-ci2(1,1))/4;
                    gammafit=coeff2body(2);
                    gammafitdev=(ci2(2,2)-ci2(2,1))/4;
                    betafit=coeff2body(3);
                    betafitdev=(ci2(3,2)-ci2(3,1))/4;
                    TIME=[0:50:max(xaxis)];
                    N2body=(initialnumber.*exp(-gammafit.*TIME))./(1+((initialnumber.*betafit./gammafit).*(1-exp(-gammafit.*TIME))));
                    lifetime1body(numberofbasenames)=1/(gammafit*timefactor); onebodydecay(numberofbasenames)=gammafit*timefactor; twobodydecay(numberofbasenames)=betafit*timefactor;
                    onebodydecayerror(numberofbasenames) = gammafitdev*timefactor; twobodydecayerror(numberofbasenames) = betafitdev*timefactor;
                    lifetime1bodyerror = sqrt((onebodydecayerror/onebodydecay^2)^2); %error in the lifetime determined by propagation of error and using error in gammafit.
                end %end density==0 (From number loss equation)
            end %end bodyloss==2 (fitting)
            
            %onebodydecay=coeff30(2);lifetime1body=1./coeff30(2);twobodydecay=0;onebodydecayerror=0;lifetime1bodyerror=0;twobodydecayerror=0;
            %Prepare to write to output file for lifetime fitting
            Lifetime(1)=onebodydecay(filecounter); %one body decay rate
            Lifetime(2)=lifetime1body(filecounter); %lifetime of atoms
            if bodyloss==2
                Lifetime(3)=twobodydecay; %two body decay rate
            end
            Lifetime(4)=onebodydecayerror(filecounter); %one body decay rate error
            Lifetime(5)=lifetime1bodyerror(filecounter); %lifetime of atoms error
            if bodyloss==2
                Lifetime(6)=twobodydecayerror(filecounter); %two body decay rate error
            end
            %Write all Lifetime parameters to outputdatafile
            dlmwrite(outputdatafile, Lifetime,'\t');
            if save2places==1
                dlmwrite(outputdatafile2, Lifetime,'\t');
            end
        elseif fittingtype==3
            %% Linewidth Fitting
            %%%Find voltage-to-frequency conversion
            if convertfreq==1
                %% Prepare data for fitting
                if wavemeteractive==1  %use if wavemeter was on
                    wavemetererror=wavemeter./wavemeter; %need error to make nlinfit happy. this will weight all pts equally b/c we have no error estimate for now.
                    InitialGuess = [(max(wavemeter)-min(wavemeter))./(min(xaxis)-max(xaxis)),max(wavemeter)]; %negative sloped guess b/c increased voltage reduces wavenumber
                    % Calculate fit parameters
                    [Pc,rc,Jc]=nlinfitweight(xaxis(:),wavemeter(:),@fitalinedangit,InitialGuess,wavemetererror(:));
                    cic=nlparci(Pc,rc,Jc);
                    slopeerror=(cic(1,2)-cic(1,1))./4;
                    % Evaluate the fit
                    pts=100;
                    for t=1:pts
                        indevariable(t) = min(xaxis)+(max(xaxis)-min(xaxis)).*(t./pts);
                        fit(t) = Pc(2)+Pc(1).*indevariable(t); %absolute value b/c sometimes complex value is put in. not sure why...
                    end
                    if plotvoltfreq==1
                        figure(100)
                        plot(xaxis,wavemeter,'kx')
                        hold on
                        plot(indevariable,fit,'-b')
                    end
                    xaxis=Pc(2)+xaxis.*Pc(1);
                    xaxis=xaxis.*30-9.9606*10^4; %convert from wavenumber to GHz
                elseif wavemeteractive==2 %use if wavemeter reading was manually read at beginning and end
                    Pc(1) = (startfreq(filecounter)-stopfreq(filecounter))./(min(xaxis)-max(xaxis)); %negative sloped guess b/c increased voltage reduces wavenumber
                    Pc(2) = startfreq(filecounter)-Pc(1).*xaxis(1);
                    xaxis=Pc(2)+xaxis.*Pc(1);
                    xaxis=xaxis.*30-9.9606*10^4; %convert from wavenumber to GHz
                else %if we load voltage to frequency conversion from file
                    [volt freq]=textread(volt2freq,'%f%f','commentstyle','matlab');
                    %volt=volt;
                    freq=81.9-(2.*abs(freq));
                    %freq=(abs(freq));
                    %81.9-2.*freq for 2008.05.12 PAS data
                    %Add 75.1 for 2008.02.07 PAS line data
                    %subtract 75.1 for 2008.02.06 feature data (we're 110
                    %MHZ blue) and make it negative.
                    %add 75.1 for 2008.01.15 feature data (we're 110 MHz red)
                    %add 106.5 MHz for 2007.12.18 data to get binding
                    %energy.
                    %add 0.0 MHz (nothing) to 2008.01.05 data to get AOM#2 freq
                    %add 98.5 to 2008.01.09 data to get binding energy
                    conversionerror=freq./freq; %need error to make nlinfit happy. this will weight all pts equally b/c we have no error estimate for now.
                    InitialGuess = [(max(freq)-min(freq))./(min(volt)-max(volt)),max(freq)]; %negative sloped guess b/c increased voltage reduces wavenumber
                    % Calculate fit parameters
                    [Pc,rc,Jc]=nlinfitweight(volt(:),freq(:),@fitalinedangit,InitialGuess,conversionerror(:));
                    cic=nlparci(Pc,rc,Jc);
                    slopeerror=(cic(1,2)-cic(1,1))./4;
                    % Evaluate the fit
                    pts=100;
                    for t=1:pts
                        freqfit(t) = min(volt)+(max(volt)-min(volt)).*(t./pts);
                        fit(t) = Pc(2)+Pc(1).*freqfit(t); %absolute value b/c sometimes complex value is put in. not sure why...
                    end
                    if plotvoltfreq==1
                        figure(100)
                        plot(volt,freq,'kx')
                        hold on
                        plot(freqfit,fit,'-b')
                        xlabel('Voltage','FontSize',fontsize,'FontWeight','bold');
                        ylabel('Frequency [MHz]','FontSize',fontsize,'FontWeight','bold');
                        title(char(volt2freq));
                        set(gca,'FontSize',fontsize,'FontWeight','bold');
                        %xlim([7 7.7]);
                        hold off
                    end
                    xaxisfreq=Pc(2)+xaxis.*Pc(1);
                end %end of wavemeter active choice
            end %of voltage-to-frequency conversion
            if fitsubset==1
                for numel=1:numberoflines
                    xaxistemp(numel)=xaxis(numel);
                    numatomstemp(numel)=numatoms(numel);
                end
                xaxisfit=xaxistemp(:);numatomsfit=numatomstemp(:);
            else
                xaxisfit=xaxis(:); %xaxisfit=xaxisfreq;
                numatomsfit=numatoms;
            end
            if singlepeakfit==1;
                %% Fit only with single Lorentzian
                %Make initial guesses
                center1 = scancenter;
                amplitude1 = -(max(numatoms)-min(numatoms))./2;
                width1 = 0.01; %rough estimate because it doesn't seem to vary that much usually without washing out
                offset = mean(numatoms(1:5));
                InitialGuesses=[amplitude1 center1 width1 offset];
                %xaxis(29:imagesize),numatoms(29:imagesize)
                [P r J] = nlinfit(xaxisfit(:),numatomsfit(:),@LorentzFitFunction1peak, InitialGuesses);
                
                amplitude1 = P(1);
                center1 = P(2);
                width1 = P(3);
                offset= P(4);
                
                ci(:,:,1)=nlparci(P(:),r,J); %error in the fit parameters 95% confidence interval
                sdev(:,1)=(ci(:,2,1)-ci(:,1,1))/4; %turn confidence interval into 1 sigma
                %sdev is matrix w/ num. of rows equal to total num. of fit param. and num.
                %of columns equal to num. of data files
                
                pts=100;
                for t=1:pts
                    indevariable(t)=min(xaxisfit)+t.*(max(xaxisfit)-min(xaxisfit))./pts;
                    fitcurve(t) = amplitude1*1./(1+(2*(indevariable(t)-center1)/width1).^2)+offset;
                end
            elseif singlepeakfit==2
                %% Fit with 2 peak Lorentzian
                %Make initial guesses
                displacement = 0.02; %amount that peaks are separated by
                locationleft = 0.63; %specify position of leftmost peak. set manually
                locationright = locationleft-displacement;
                amplitude = (max(numatoms)-min(numatoms));
                width = 0.03; %rough estimate of linewidth because it doesn't seem to vary that much usually without washing out
                offsetguess = mean(numatoms(1:5));
                InitialGuesses=[amplitude locationleft width amplitude locationright width offsetguess];
                
                [P r J] = nlinfit(xaxis,numatoms,@LorentzFitFunction2peaks, InitialGuesses);
                amplitude1 = P(1);
                center1 = P(2);
                width1 = P(3);
                amplitude2 = P(4);
                center2 = P(5);
                width2 = P(6);
                offset = P(7);
                
                ci(:,:,1)=nlparci(P(:),r,J); %error in the fit parameters 95% confidence interval
                sdev(:,1)=(ci(:,2,1)-ci(:,1,1))/4; %turn confidence interval into 1 sigma
                %sdev is matrix w/ num. of rows equal to total num. of fit param. and num.
                %of columns equal to num. of data files
                
                pts=100;
                for t=1:pts
                    indevariable(t)=min(xaxis)+t.*(max(xaxis)-min(xaxis))./pts;
                    fitcurve(t) = amplitude1*1./(1+(2*(indevariable(t)-center1)/width1).^2)+amplitude2*1./(1+(2*(indevariable(t)-center2)/width2).^2)+offset;
                end
            elseif singlepeakfit==3
                %% Fit with 3 peak Lorentzian
                %Make initial guesses
                displacement = .4;
                locationmiddle = 164.3; %center frequency of the whole thing. set manually.
                locationleft = locationmiddle-displacement;
                locationright = locationmiddle+displacement;
                amplitude = (max(numatoms)-min(numatoms));
                %amplitude=-1;
                width = 0.2; %rough estimate because it doesn't seem to vary that much usually without washing out
                offsetguess = mean(numatoms(1:5));
                InitialGuesses=[amplitude locationleft width amplitude locationmiddle width amplitude locationright width offsetguess];
                
                [P r J] = nlinfit(xaxis,numatoms,@LorentzFitFunction3peaks, InitialGuesses);
                amplitude1 = P(1);
                center1 = P(2);
                width1 = P(3);
                amplitude2 = P(4);
                center2 = P(5);
                width2 = P(6);
                amplitude3 = P(7);
                center3 = P(8);
                width3 = P(9);
                offset = P(10);
                
                ci(:,:,1)=nlparci(P(:),r,J); %error in the fit parameters 95% confidence interval
                sdev(:,1)=(ci(:,2,1)-ci(:,1,1))/4; %turn confidence interval into 1 sigma
                %sdev is matrix w/ num. of rows equal to total num. of fit param. and num.
                %of columns equal to num. of data files
                
                pts=100;
                for t=1:pts
                    indevariable(t)=min(xaxis)+t.*(max(xaxis)-min(xaxis))./pts;
                    fitcurve(t) = amplitude1*1./(1+(2*(indevariable(t)-center1)/width1).^2)+amplitude2*1./(1+(2*(indevariable(t)-center2)/width2).^2)+amplitude3*1./(1+(2*(indevariable(t)-center3)/width3).^2)+offset;
                end
            end %end of singlepeakfit==1
            %%%Write output file for linewidth fitting
            dlmwrite(outputdatafile,P,'\t');
            if save2places==1
                dlmwrite(outputdatafile2,P,'\t');
            end
        end %of fit sizes/numbers section
    end %of fitting section
    
    if convertfreq==0 || fittingtype ~=3
        xaxisfreq=xaxis; %temporary place for this to solve programming issue; resolve this!!!
    end
    %% Plotting section
    if readaveragedinput==1
        figure(filecounter)
        semilogy(xaxisfreq,numatoms,'.b')
        hold on
        
        if fitting==1 && fittingtype==2
            if bodyloss==2
                plot(TIME,N2body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
            elseif bodyloss==1
                plot(TIME,N1body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
            end
            hold on
            
            text(mean(xaxisfreq)/1.5,1*mean(numatoms),strcat(char('\tau ='),num2str(lifetime1body(filecounter),3),'\pm',num2str(lifetime1bodyerror(filecounter),2), ' s'),'FontSize',textsize,'FontWeight','bold');
        elseif fitting==1 && fittingtype==3
            plot(indevariable,fitcurve,'Color',[1 0 0],'LineStyle','-','LineWidth',1);
            hold on
            text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.8,ztext,strcat('Width=',num2str(width1),'+/-',num2str(sdev(3))),'FontSize',textsize);
            text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.6,ztext,strcat('Center=',num2str(center1),'+/-',num2str(sdev(2))),'FontSize',textsize);
            text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.4,ztext,strcat('FracAmpl=',num2str(abs(amplitude1/offset))),'FontSize',textsize);
        end
        xlim([min(xaxisfreq) max(xaxisfreq)]);
        ylim([0+(min(numatoms)-min(numatoms)/4) max(numatoms)+max(numatoms)/4]);
        xlabel('Hold Time [ms]','FontSize',fontsize,'FontWeight','bold');
        ylabel('# Atoms (10^6)','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        title(basehead,'FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        hold on
    else
        %xaxisfreq=xaxis;
        xtext=median(xaxisfreq)-(max(xaxisfreq)-min(xaxisfreq))/2;
        ytext = max(numatoms)-1.5.*(max(numatoms)-min(numatoms));
        ztext=0;
        figure(filecounter)
        %%%Plot 1 is number of atoms vs changed value
        subplot(2,2,1)
        %errorbar(xaxisfreq,numatoms,numatomserror,'xb')
        plot(xaxisfreq,numatoms,'.b')
        hold on
        GRID On
        if fitting==1 && fittingtype==2
            if bodyloss==2
                plot(TIME,N2body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
            elseif bodyloss==1
                if plotsubset==1
                    plot(TIME,N1body20,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
                else
                    plot(TIME,N1body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
                end
            end
            hold on
            if plotsubset==1
                text(mean(xaxisfreq)/1.5,1.5*mean(numatoms),strcat(char('\tau ='),num2str(1./coeff20(2)), ' s'),'FontSize',textsize,'FontWeight','bold');
            else
                text(mean(xaxisfreq)/1.5,1*mean(numatoms),strcat(char('\tau ='),num2str(lifetime1body(filecounter),3),'\pm',num2str(lifetime1bodyerror(filecounter),2), ' s'),'FontSize',textsize,'FontWeight','bold');
            end
        elseif fitting==1 && fittingtype==3
            plot(indevariable,fitcurve,'Color',[1 0 0],'LineStyle','-','LineWidth',1);
            hold on
            %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.8,ztext,strcat('Width=',num2str(width1),'+/-',num2str(sdev(3))),'FontSize',textsize);
            %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.6,ztext,strcat('Center=',num2str(center1),'+/-',num2str(sdev(2))),'FontSize',textsize);
            %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.4,ztext,strcat('FracAmpl=',num2str(abs(amplitude1/offset))),'FontSize',textsize);
        end
        %xlim([0+(min(xaxis)-abs(min(xaxis))/4) max(xaxis)+max(xaxis)/4]);
        %xlim([min(xaxis) max(xaxis)]);
%        xlim([min(xaxisfreq) max(xaxisfreq)]);
        %xlim([162 166]);
        ylim([0+(min(numatoms)-min(numatoms)/4) max(numatoms)+max(numatoms)/4]);
        %ylim([14 20]);
        xlabel('Frequency [Hz]','FontSize',fontsize,'FontWeight','bold');
        xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        %xlabel('some freq (not AOM#1) [MHz]','FontSize',fontsize,'FontWeight','bold');
        ylabel('# Atoms (10^6)','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        title(basehead,'FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        hold off
        %%%Plot 2 is Peak OD vs changed value
        subplot(2,2,2)
        plot(xaxis,peakOD,'db');
        hold on
        GRID On
%        xlim([min(xaxis) max(xaxis)]);
        %xlim([162 166]);
        ylim([0 max(peakOD)+max(peakOD)/4]);
        xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        ylabel('Peak OD','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        hold off
        %%%%Plot 3 is x size vs changed value
        subplot(2,2,3)
        errorbar(xaxis,sigmax.*sizefactor,sigmaxerror.*sizefactor,'xb');
        hold on
        errorbar(xaxis,sigmay.*sizefactor,sigmayerror.*sizefactor,'ob');
        hold on
        GRID On
        if fitting==1 && fittingtype==1
            plot(indevariable,fitx.*sizefactor,'-r');
            hold on
            plot(indevariable,fity.*sizefactor,'-r');
            hold on
            text(mean(xaxis)/4,tempconstant*max(sigmax.*sizefactor)-min(sigmax.*sizefactor),strcat(char('Tx='),num2str(xtemperature,3), ' mK'),'FontSize',textsize,'FontWeight','bold');
            text(mean(xaxis)/4,tempconstant*max(sigmay.*sizefactor)-min(sigmay.*sizefactor),strcat(char('Ty='),num2str(ytemperature,3), ' mK'),'FontSize',textsize,'FontWeight','bold');
        end
%        xlim([min(xaxis) max(xaxis)]);
        %xlim([0 15000]);
        %ylim([min(sigmax.*sizefactor)-min(sigmax.*sizefactor)/4 max(sigmax.*sizefactor)+max(sigmax.*sizefactor)/4]);
        %ylim([0 0.15]);
        xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        ylabel('\sigma_x & \sigma_y (mm)','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        title('Sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        hold off
        %Plot 4 is y size vs changed value
        subplot(2,2,4)
        ysizetemp=(sigmay(:).^2.*mass)./(kBoltz.*(delaytime(:)./timefactor).^2)./templow;
        xsizetemp=(sigmax(:).^2.*mass)./(kBoltz.*(delaytime(:)./timefactor).^2)./templow;
        numatoms'
        plot(xaxis,ysizetemp,'ob');
        hold on
        plot(xaxis,xsizetemp,'xb');
        hold on
        GRID On
%        xlim([min(xaxis) max(xaxis)]);
        %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
        %ylim([0 12]);
        xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        ylabel('Temp (\muK)','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        title('Temp. calc. from sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        
        
%         subplot(3,2,5)
         
%         plot(xaxisfreq,PhaseSpaceDensity,'-r')
%         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
%         ylabel('PSD','FontSize',fontsize,'FontWeight','bold');
%         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        %title('Temp. calc. from sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        hold off
        
        if psd==1
            figure(100+filecounter)
            subplot(2,2,1)
            plot(xaxisfreq,numatoms,'.b')
            hold on
            GRID On
            if fitting==1 && fittingtype==2
                if bodyloss==2
                    plot(TIME,N2body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
                elseif bodyloss==1
                    if plotsubset==1
                        plot(TIME,N1body20,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
                    else
                        plot(TIME,N1body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
                    end
                end
                hold on
                if plotsubset==1
                    text(mean(xaxisfreq)/1.5,2*mean(numatoms),strcat(char('\tau ='),num2str(1./coeff20(2)), ' s'),'FontSize',textsize,'FontWeight','bold');
                else
                    text(mean(xaxisfreq)/1.5,2*mean(numatoms),strcat(char('\tau ='),num2str(lifetime1body(filecounter),3),'\pm',num2str(lifetime1bodyerror(filecounter),2), ' s'),'FontSize',textsize,'FontWeight','bold');
                end
            elseif fitting==1 && fittingtype==3
                plot(indevariable,fitcurve,'Color',[1 0 0],'LineStyle','-','LineWidth',1);
                hold on
                %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.8,ztext,strcat('Width=',num2str(width1),'+/-',num2str(sdev(3))),'FontSize',textsize);
                %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.6,ztext,strcat('Center=',num2str(center1),'+/-',num2str(sdev(2))),'FontSize',textsize);
                %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.4,ztext,strcat('FracAmpl=',num2str(abs(amplitude1/offset))),'FontSize',textsize);
            end
            %xlim([0+(min(xaxis)-abs(min(xaxis))/4) max(xaxis)+max(xaxis)/4]);
            %xlim([min(xaxis) max(xaxis)]);
            xlim([min(xaxisfreq) max(xaxisfreq)]);
            %xlim([162 166]);
            ylim([0+(min(numatoms)-min(numatoms)/4) max(numatoms)+max(numatoms)/4]);
            %ylim([14 20]);
            xlabel('Frequency [Hz]','FontSize',fontsize,'FontWeight','bold');
            xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
            %xlabel('some freq (not AOM#1) [MHz]','FontSize',fontsize,'FontWeight','bold');
            ylabel('# Atoms (10^6)','FontSize',fontsize,'FontWeight','bold');
            set(gca,'FontSize',fontsize-4,'FontWeight','bold');
            title(basehead,'FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
            hold off
            subplot(2,2,2)
            plot(xaxisfreq,psdx,'xb')
            hold on
            plot(xaxisfreq,psdy,'ob')
            GRID On
            xlabel('Frequency [Hz]','FontSize',fontsize,'FontWeight','bold');
            ylabel('Phase Space Density','FontSize',fontsize,'FontWeight','bold');
        end
        
        if pos==1
            figure(200+filecounter)
            subplot(2,2,1)
            plot(xaxisfreq,numatoms,'.b')
            hold on
            GRID On
            if fitting==1 && fittingtype==2
                if bodyloss==2
                    plot(TIME,N2body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
                elseif bodyloss==1
                    if plotsubset==1
                        plot(TIME,N1body20,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
                    else
                        plot(TIME,N1body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
                    end
                end
                hold on
                if plotsubset==1
                    text(mean(xaxisfreq)/1.5,2*mean(numatoms),strcat(char('\tau ='),num2str(1./coeff20(2)), ' s'),'FontSize',textsize,'FontWeight','bold');
                else
                    text(mean(xaxisfreq)/1.5,2*mean(numatoms),strcat(char('\tau ='),num2str(lifetime1body(filecounter),3),'\pm',num2str(lifetime1bodyerror(filecounter),2), ' s'),'FontSize',textsize,'FontWeight','bold');
                end
            elseif fitting==1 && fittingtype==3
                plot(indevariable,fitcurve,'Color',[1 0 0],'LineStyle','-','LineWidth',1);
                hold on
                %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.8,ztext,strcat('Width=',num2str(width1),'+/-',num2str(sdev(3))),'FontSize',textsize);
                %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.6,ztext,strcat('Center=',num2str(center1),'+/-',num2str(sdev(2))),'FontSize',textsize);
                %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.4,ztext,strcat('FracAmpl=',num2str(abs(amplitude1/offset))),'FontSize',textsize);
            end
            %xlim([0+(min(xaxis)-abs(min(xaxis))/4) max(xaxis)+max(xaxis)/4]);
            %xlim([min(xaxis) max(xaxis)]);
            xlim([min(xaxisfreq) max(xaxisfreq)]);
            %xlim([162 166]);
            ylim([0+(min(numatoms)-min(numatoms)/4) max(numatoms)+max(numatoms)/4]);
            %ylim([14 20]);
            xlabel('Frequency [Hz]','FontSize',fontsize,'FontWeight','bold');
            xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
            %xlabel('some freq (not AOM#1) [MHz]','FontSize',fontsize,'FontWeight','bold');
            ylabel('# Atoms (10^6)','FontSize',fontsize,'FontWeight','bold');
            set(gca,'FontSize',fontsize-4,'FontWeight','bold');
            title(basehead,'FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
            hold off
            subplot(2,2,2)
            plot(xaxisfreq(:),positionx(:),'xb')
            hold on
            plot(xaxisfreq(:),positiony(:),'ob')
            GRID On
            xlim([0 900]);
            xlabel('Hold time [ms]','FontSize',fontsize,'FontWeight','bold');
            ylabel('Cloud Position [m]','FontSize',fontsize,'FontWeight','bold');
            subplot(2,2,3)
            plot(xaxisfreq(:),positionx(:),'xb')
            hold on
            plot(xaxisfreq(:),positiony(:),'ob')
            GRID On
            xlim([0 900]);ylim([0.0004 0.0006]);
            xlabel('Hold time [ms]','FontSize',fontsize,'FontWeight','bold');
            ylabel('Cloud Position [m]','FontSize',fontsize,'FontWeight','bold');
        end
        
        %         %%%%Plot 3 is x size vs changed value
        %         figure(557)
        %          %%%Plot 1 is number of atoms vs changed value
        %         %errorbar(xaxisfreq,numatoms,numatomserror,'xb')
        %         plot(xaxisfreq,numatoms,'.b')
        %         hold on
        %         GRID On
        %         if fitting==1 && fittingtype==2
        %             if bodyloss==2
        %                 plot(TIME,N2body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
        %             elseif bodyloss==1
        %                 if plotsubset==1
        %                     plot(TIME,N1body20,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
        %                 else
        %                 plot(TIME,N1body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
        %                 end
        %             end
        %             hold on
        %              if plotsubset==1
        %                  text(mean(xaxisfreq)/1.5,2*mean(numatoms),strcat(char('\tau ='),num2str(1./coeff20(2)), ' s'),'FontSize',textsize,'FontWeight','bold');
        %                           else
        %             text(mean(xaxisfreq)/1.5,1*mean(numatoms),strcat(char('\tau ='),num2str(lifetime1body(filecounter),3),'\pm',num2str(lifetime1bodyerror(filecounter),2), ' s'),'FontSize',textsize,'FontWeight','bold');
        %              end
        %         elseif fitting==1 && fittingtype==3
        %             plot(indevariable,fitcurve,'Color',[1 0 0],'LineStyle','-','LineWidth',1);
        %             hold on
        %             %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.8,ztext,strcat('Width=',num2str(width1),'+/-',num2str(sdev(3))),'FontSize',textsize);
        %             %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.6,ztext,strcat('Center=',num2str(center1),'+/-',num2str(sdev(2))),'FontSize',textsize);
        %             %             text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.4,ztext,strcat('FracAmpl=',num2str(abs(amplitude1/offset))),'FontSize',textsize);
        %         end
        %         %xlim([0+(min(xaxis)-abs(min(xaxis))/4) max(xaxis)+max(xaxis)/4]);
        %         %xlim([min(xaxis) max(xaxis)]);
        %         xlim([min(xaxisfreq) max(xaxisfreq)]);
        %         %xlim([162 166]);
        %         ylim([0+(min(numatoms)-min(numatoms)/4) max(numatoms)+max(numatoms)/4]);
        %         %ylim([14 20]);
        %         xlabel('Frequency [Hz]','FontSize',fontsize,'FontWeight','bold');
        %         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        %         %xlabel('some freq (not AOM#1) [MHz]','FontSize',fontsize,'FontWeight','bold');
        %         ylabel('# Atoms (10^6)','FontSize',fontsize,'FontWeight','bold');
        %         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        %         title(basehead,'FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        %         hold off
        
    end %END if readaveragedinput==1 LOOP
    if composite==1
        COLORS=[0 0 0;1 0 0;0 1 0; 0 0 1; .25 .25 .25; 0 .75 0; 0 .75 .75; .5 0 .5; .75 0 .75;0 0 1;0 0 0;1 0 .25;0 .75 0; .5 0 .5; 0 1 0; 0 .75 0; 0 .75 .75; .5 0 .5; .75 0 .75;0 0 1;0 0 0;1 0 .25;0 .75 0; .5 0 .5; 0 1 0; 0 .75 0; 0 .75 .75; .5 0 .5; .75 0 .75;0 0 1;0 0 0;1 0 .25;0 .75 0; .5 0 .5; 0 1 0; 0 .75 0; 0 .75 .75; .5 0 .5; .75 0 .75;0 0 1];
        MARKERS={'-o','-s','-^','-d','-v','-p','-<','-h','->','-x','-*','-+','-.','-o','-s','-^','-d','-v','-p','-<','-h','->','-x','-*','-+','-.','-o','-s','-^','-d','-v','-p','-<','-h','->','-x','-*','-+','-.','-o','-s','-^','-d','-v','-p','-<','-h','->','-x','-*','-+','-.'};
        
        figure(555)
%         subplot(2,2,2)
%         plot(xaxisfreq,xsizetemp,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         hold on
%         GRID On
%         %xlim([min(xaxis) max(xaxis)]);
%         %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
%         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
%         ylabel('x axis temp (\muK)','FontSize',fontsize,'FontWeight','bold');
%         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
%         %title('Temp. calc. from sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
%         %str(rfpower(5))
        xaxisfreq=xaxisfreq+timeoffset(basenamenumber);

        subplot(3,3,1)
        %semilogx(xaxisfreq,numatoms,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        plot(xaxisfreq,numatoms,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        %plot(xaxisfreq,numatoms./max(numatoms),MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        %text(mean(xaxisfreq)/1.35,mean(numatoms),num2str(timevector(basenamenumber)),'FontSize',textsize,'Color',COLORS(basenamenumber,:));
%         if basenamenumber==1
%             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR On','FontSize',textsize,'Color',COLORS(basenamenumber,:));
%         elseif basenamenumber==2
%             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR Off','FontSize',textsize,'Color',COLORS(basenamenumber,:));
%         end
        hold on
        GRID On
        %xlim([min(xaxis) max(xaxis)]);
        %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
        xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        %ylabel('Num atoms [10^6]','FontSize',fontsize,'FontWeight','bold');
        ylabel('Num of atoms [10^6]','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        if OFRon(filecounter)==1
            stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr dual species-',num2str(timevector(basenamenumber)));
        else
            stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr alone-',num2str(timevector(basenamenumber)));
        end
        if filecounter==numberofbasenames
            l1=legend(stringmatrix,'Location','Best');
        end
        %ylim([0 1]);
        %xmax=50;
        %xlim([0 xmax]);
        %stringmatrix=get(l1,'String');
        %set(l1,'String',stringmatrix);
        %title('Temp. calc. from sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        
        subplot(3,3,2)
%         semilogx(xaxisfreq,ysizetemp,'o','Color',COLORS(basenamenumber,:),'MarkerSize',7);
%         plot(xaxisfreq,ysizetemp,'o','Color',COLORS(basenamenumber,:),'MarkerSize',7);
%         hold on
%         semilogx(xaxisfreq,xsizetemp,'x','Color',COLORS(basenamenumber,:),'MarkerSize',7);
%         plot(xaxisfreq,xsizetemp,'x','Color',COLORS(basenamenumber,:),'MarkerSize',7);
        plot(xaxisfreq,(xsizetemp+ysizetemp)./2,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7);
        hold on
        %xsizetemp
        %plot(xaxisfreq,sigmax.*sizefactor,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        %hold on
        %plot(xaxisfreq,sigmax.*sizefactor,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        GRID On
        %xlim([min(xaxis) max(xaxis)]);
        %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
        %ylim([0 1]);
        %xlim([8000 9000]);
        xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        ylabel('Temperatures [\muK]','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        %title('Temp. calc. from sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        %str(rfpower(5))

        subplot(3,3,3)
        ODTCalFactor=1.122; % Power of the input beam = ODTCalFactor * DAC Voltage, where 1.122 is the calibration taken on Jan. 23, 2010
        ODTpower=ODTvoltage(:).*ODTCalFactor;
%         fbarModel=-8947.61078+25026.69095.*ODTpower-26177.9322.*ODTpower.
%         ^2+12205.76901.*ODTpower.^3-2136.*ODTpower.^4; % for short range
        fbarModel=-189.6739+419.55801.*ODTpower-278.879.*ODTpower.^2+106.32647.*ODTpower.^3-23.68352.*...
            ODTpower.^4+3.06336.*ODTpower.^5-0.21312.*ODTpower.^6+0.00616.*ODTpower.^7;
      %  Tfermi=((xsizetemp+ysizetemp)./2).*10^(-6).*kBoltz./((6.*numatoms(:).*(10^6)./10).^(1/3).*hbar.*(2.*pi.*fbarModel));
        %semilogx(xaxisfreq,numatoms,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        %plot(xaxisfreq,Tfermi,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7);
        %plot(xaxisfreq,numatoms./max(numatoms),MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        %text(mean(xaxisfreq)/1.35,mean(numatoms),num2str(timevector(basenamenumber)),'FontSize',textsize,'Color',COLORS(basenamenumber,:));
%         if basenamenumber==1
%             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR On','FontSize',textsize,'Color',COLORS(basenamenumber,:));
%         elseif basenamenumber==2
%             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR Off','FontSize',textsize,'Color',COLORS(basenamenumber,:));
%         end
        hold on
        GRID On
        %xlim([min(xaxis) max(xaxis)]);
        %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
        xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        %ylabel('Num atoms
        %[10^6]','FontSize',fontsize,'FontWeight','bold');
        ylabel('T / T_f','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
%         if OFRon(filecounter)==1
%             stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr dual species-',num2str(timevector(basenamenumber)));
%         else
%             stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr alone-',num2str(timevector(basenamenumber)));
%         end
%         if filecounter==numberofbasenames
%             l1=legend(stringmatrix,'Location','Best');
%         end
        
        subplot(3,3,4)
        ODTCalFactor=1.122; % Power of the input beam = ODTCalFactor * DAC Voltage, where 1.122 is the calibration taken on Jan. 23, 2010
%         ODTpower=ODTvoltage(:).*ODTCalFactor;
%         % PhaseSpaceDensity=SrPSD(ODTpower(:),xsizetemp(:).*10.^(-6),numatoms(:).*10.^6, 86);
%         for n=1:length(ODTpower)  %It did not work for "PhaseSpaceDensity=SrPSD(ODTpower(:),10.^(-6),numatoms(:).*10.^6,86);" which may be ascribed to the "if" code in the V1eff function, so I employed this "stupid" method to evaluate PSD
%             PhaseSpaceDensity(n)=SrPSD20091222(ODTpower(n),xsizetemp(n).*10.^(-6),numatoms(n).*10.^6,87);
%         end
%         semilogy(xaxisfreq,PhaseSpaceDensity,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        hold on
        GRID On
        %xlim([min(xaxis) max(xaxis)]);
        %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
        xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        ylabel('PSD','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        %xlim([0 xmax]);
        %ylim([10^-2 10^-1]);
        xlim('auto');
        ylim('auto');
        %title('Temp. calc. from sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        clear PhaseSpaceDensity
        
        subplot(3,3,5)
        % PhaseSpaceDensity=(numatoms(:).*ODTpower.^(3/2))./xsizetemp(:).^3;
        plot(xaxisfreq,sigmax,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        hold on
        plot(xaxisfreq,sigmay,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        hold on
        GRID On
        %xlim([min(xaxis) max(xaxis)]);
        %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
        xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
        ylabel('Size [m]','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
        %xlim([0 xmax]);
        %ylim([0 1.5*10^-4]);
        %xlim('auto');
        ylim('auto');
        %title('Temp. calc. from
        %sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
        %clear PowerInput ODTpower 
        
        subplot(3,3,6)
        ODTCalFactor=1.122; % Power of the input beam = ODTCalFactor * DAC Voltage, where 1.122 is the calibration taken on Jan. 23, 2010
        ODTpower=ODTvoltage(:).*ODTCalFactor;
%         fbarModel=-8947.61078+25026.69095.*ODTpower-26177.9322.*ODTpower.
%         ^2+12205.76901.*ODTpower.^3-2136.*ODTpower.^4; % for short range
        fbarModel=-189.6739+419.55801.*ODTpower-278.879.*ODTpower.^2+106.32647.*ODTpower.^3-23.68352.*...
            ODTpower.^4+3.06336.*ODTpower.^5-0.21312.*ODTpower.^6+0.00616.*ODTpower.^7;
        Tfermi=((xsizetemp+ysizetemp)./2).*10^(-6).*kBoltz./((6.*numatoms(:).*(10^6)./10).^(1/3).*hbar.*(2.*pi.*fbarModel));
        %semilogx(xaxisfreq,numatoms,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        plot(numatoms,Tfermi,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7);
        %plot(xaxisfreq,numatoms./max(numatoms),MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
        %text(mean(xaxisfreq)/1.35,mean(numatoms),num2str(timevector(basenamenumber)),'FontSize',textsize,'Color',COLORS(basenamenumber,:));
%         if basenamenumber==1
%             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR On','FontSize',textsize,'Color',COLORS(basenamenumber,:));
%         elseif basenamenumber==2
%             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR Off','FontSize',textsize,'Color',COLORS(basenamenumber,:));
%         end
        hold on
        GRID On
        %xlim([min(xaxis) max(xaxis)]);
        %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
        xlabel('Num of atoms [10^6]','FontSize',fontsize,'FontWeight','bold');
        %ylabel('Num atoms [10^6]','FontSize',fontsize,'FontWeight','bold');
        ylabel('T / T_f','FontSize',fontsize,'FontWeight','bold');
        set(gca,'FontSize',fontsize-4,'FontWeight','bold');
%         if OFRon(filecounter)==1
%             stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr dual species-',num2str(timevector(basenamenumber)));
%         else
%             stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr alone-',num2str(timevector(basenamenumber)));
%         end
%         if filecounter==numberofbasenames
%             l1=legend(stringmatrix,'Location','Best');
%         end
        
         subplot(3,3,7)
         ODTCalFactor=1.122; % Power of the input beam = ODTCalFactor * DAC Voltage
         ODTpower=ODTvoltage(:).*ODTCalFactor;

         % Calibration Before Dec 22, 2009
         %          trapdepth=-7.11787+7.20234.*ODTpower-0.52412.*ODTpower.^2-0.31157.*ODTpower.^3+0.13548.*ODTpower.^4-...
%              0.02166.*ODTpower.^5+0.0016.*ODTpower.^6-0.000045194.*ODTpower.^7; 

% Calibraition On Dec. 22, 2009
         trapdepth=-19.09089+34.08091.*ODTpower-23.08458.*ODTpower.^2+9.80717.*ODTpower.^3-2.55853.*ODTpower.^4+...
             0.427.*ODTpower.^5-0.04577.*ODTpower.^6+0.00305.*ODTpower.^7-0.000114988.*ODTpower.^8+0.00000187438.*ODTpower.^9;
         
         %trapdepthNoGravity=4.00763.*ODTpower;

         eta=trapdepth./xsizetemp;
         plot(xaxisfreq,trapdepth,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
         %plot(xaxisfreq,trapdepthNoGravity,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
         hold on
         GRID On
         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
         ylabel('Trap Depth [\muK]','FontSize',fontsize,'FontWeight','bold');
         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
         %xlim([0 xmax]);
         %ylim([0 3*10^-4]);
         xlim('auto');
         clear trapdepth
        
                 
         subplot(3,3,8)
         plot(xaxisfreq,eta,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
         hold on
         GRID On
         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
         ylabel('\eta','FontSize',fontsize,'FontWeight','bold');
         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
         %xlim([0 xmax]);
         %ylim([0 3*10^-4]);
         xlim('auto');
         
         subplot(3,3,9)
         plot(xaxisfreq,ODTvoltage,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
         hold on
         GRID On
         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
         ylabel('DAC [Voltage]','FontSize',fontsize,'FontWeight','bold');
         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
         %xlim([0 xmax]);
         %ylim([0 3*10^-4]);
         xlim('auto');
         
%           figure(888)
% 
%         xaxisfreq=xaxisfreq+timeoffset(basenamenumber);
% 
%         subplot(2,2,1)
%         %semilogx(xaxisfreq,numatoms,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         plot(xaxisfreq,numatoms,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         %plot(xaxisfreq,numatoms./max(numatoms),MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         %text(mean(xaxisfreq)/1.35,mean(numatoms),num2str(timevector(basenamenumber)),'FontSize',textsize,'Color',COLORS(basenamenumber,:));
% %         if basenamenumber==1
% %             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR On','FontSize',textsize,'Color',COLORS(basenamenumber,:));
% %         elseif basenamenumber==2
% %             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR Off','FontSize',textsize,'Color',COLORS(basenamenumber,:));
% %         end
%         hold on
%         GRID On
%         %xlim([min(xaxis) max(xaxis)]);
%         %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
%         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
%         ylabel('Num atoms [10^6]','FontSize',fontsize,'FontWeight','bold');
%         %ylabel('Normalized Num atoms [10^6]','FontSize',fontsize,'FontWeight','bold');
%         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
%         if OFRon(filecounter)==1
%             stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr dual species-',num2str(timevector(basenamenumber)));
%         else
%             stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr alone-',num2str(timevector(basenamenumber)));
%         end
%         if filecounter==numberofbasenames
%             l1=legend(stringmatrix,'Location','Best');
%         end
%         %ylim([0 1]);
%         %xmax=50;
%         %xlim([0 xmax]);
%         %stringmatrix=get(l1,'String');
%         %set(l1,'String',stringmatrix);
%         %title('Temp. calc. from sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
%         
%         subplot(2,2,2)
% %         semilogx(xaxisfreq,ysizetemp,'o','Color',COLORS(basenamenumber,:),'MarkerSize',7);
%         plot(xaxisfreq,ysizetemp,'o','Color',COLORS(basenamenumber,:),'MarkerSize',7);
%         hold on
% %         semilogx(xaxisfreq,xsizetemp,'x','Color',COLORS(basenamenumber,:),'MarkerSize',7);
%         plot(xaxisfreq,xsizetemp,'x','Color',COLORS(basenamenumber,:),'MarkerSize',7);
%         %xsizetemp
%         %plot(xaxisfreq,sigmax.*sizefactor,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         %hold on
%         %plot(xaxisfreq,sigmax.*sizefactor,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         GRID On
%         %xlim([min(xaxis) max(xaxis)]);
%         %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
%         %ylim([0 1]);
%         %xlim([8000 9000]);
%         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
%         ylabel('Temperatures [\muK]','FontSize',fontsize,'FontWeight','bold');
%         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
%         %title('Temp. calc. from sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
%         %str(rfpower(5))
% 
%         subplot(2,2,3)
%         %semilogx(xaxisfreq,numatoms,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         %plot(xaxisfreq,numatoms,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         plot(xaxisfreq,numatoms./max(numatoms),MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         %text(mean(xaxisfreq)/1.35,mean(numatoms),num2str(timevector(basenamenumber)),'FontSize',textsize,'Color',COLORS(basenamenumber,:));
% %         if basenamenumber==1
% %             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR On','FontSize',textsize,'Color',COLORS(basenamenumber,:));
% %         elseif basenamenumber==2
% %             text(mean(xaxisfreq)/1.35,mean(numatoms),'OFR Off','FontSize',textsize,'Color',COLORS(basenamenumber,:));
% %         end
%         hold on
%         GRID On
%         %xlim([min(xaxis) max(xaxis)]);
%         %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
%         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
%         %ylabel('Num atoms [10^6]','FontSize',fontsize,'FontWeight','bold');
%         ylabel('Normalized Num atoms [10^6]','FontSize',fontsize,'FontWeight','bold');
%         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
%         if OFRon(filecounter)==1
%             stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr dual species-',num2str(timevector(basenamenumber)));
%         else
%             stringmatrix{filecounter}=strcat(num2str(rfpower(basenamenumber)),'Sr alone-',num2str(timevector(basenamenumber)));
%         end
%         if filecounter==numberofbasenames
%             l1=legend(stringmatrix,'Location','Best');
%         end
%         
%         subplot(2,2,4)
%         % PhaseSpaceDensity=(numatoms(:).*ODTpower.^(3/2))./xsizetemp(:).^3;
%         plot(xaxisfreq,sigmax,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         hold on
%         plot(xaxisfreq,sigmay,MARKERS{basenamenumber},'Color',COLORS(basenamenumber,:),'MarkerSize',7)
%         hold on
%         GRID On
%         %xlim([min(xaxis) max(xaxis)]);
%         %ylim([max(xsizetemp)-min(xsizetemp)/4 max(xsizetemp)+max(xsizetemp)/4]);
%         xlabel(xaxislabel,'FontSize',fontsize,'FontWeight','bold');
%         ylabel('Size [m]','FontSize',fontsize,'FontWeight','bold');
%         set(gca,'FontSize',fontsize-4,'FontWeight','bold');
%         %xlim([0 xmax]);
%         %ylim([0 1.5*10^-4]);
%         %xlim('auto');
%         ylim('auto');
%         %title('Temp. calc. from
%         %sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
%         %clear PowerInput ODTpower      
                  
    end %end composite==1
    
    %     figure(1035+filecounter)
    %     %subplot(2,1,1)
    %         plot(xaxisfreq,numatoms,'*b','MarkerSize',10)
    %         hold on
    %          plot(TIME,N1body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',2);
    %          hold on
    %           text(mean(xaxisfreq)/1.5,1*mean(numatoms),strcat(char('\tau ='),num2str(lifetime1body(filecounter),3), ' s'),'FontSize',textsize+2);
    %         xlim([min(xaxis) max(xaxis)]);
    %         ylabel('Number of atoms','FontSize',fontsize+2);
    %         set(gca,'FontSize',fontsize);
    %         hold on
    %         subplot(2,1,2)
    %         errorbar(xaxis,sigmay.*sizefactor,sigmayerror.*sizefactor,'ob');
    %         hold on
    %         if fitting==1 && fittingtype==1
    %             plot(indevariable,fity.*sizefactor,'-k');
    %             hold on
    %             text(mean(xaxis)/4,tempconstant*max(sigmay.*sizefactor)-min(sigmay.*sizefactor),strcat(char('Ty='),num2str(ytemperature*10^3,3), ' \muK'),'FontSize',textsize+2);
    %         end
    %         xlim([min(xaxis) max(xaxis)]);
    %xlabel(xaxislabel,'FontSize',fontsize+2);
    %ylabel('\sigma_y (mm)','FontSize',fontsize+2);
    %set(gca,'FontSize',fontsize);
    %title('Sizes','FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
    %hold off
    
    
end %End loop through batchlist

clear numatomssum sigmaxsum sigmaysum filecounter
filecounter=0

%if sumfiles==1 && powerW==scancenter  %%%Loop through the basenames AGAIN in cases where we want to sum separate files.
if sumfiles==1
    if powerW==scancenter
        for basenamenumber=1:numberofbasenames
            filecounter=filecounter+1
            
            clear basehead xaxistemp peakODtemp sigmaxtemp sigmaytemp numatomserrortemp
            basehead=char(basenamevector(basenamenumber));
            
            clear e1 xaxistemp e3 e4 e5 e6 e7 e8 e9 e10 e11 e12
            clear f1 f2 f3 f4 f5 f6 positionx positiony f9 peakODtemp sigmaxtemp sigmaytemp f13 f14 f15 f16 f17 f18 sigmaxerrortemp sigmayerrortemp f21 f22 f23 f24 f25 numatomserrortemp numatomstemp chi3Dtemp sizeparametertemp
            %INPUT FILES
            xaxisfile=char(strcat(deblank(directoryvector(filecounter)),deblank(basehead),deblank('.batch')));
            imagefile2D=char(strcat(deblank(directoryvector(filecounter)),deblank(basehead),deblank('2Dfitparams.txt')));
            %We use changed values from batch file and number of atoms from 2D fit parameters file (output of imagefit).
            [e1 xaxistemp e3 e4 e5 e6 e7 e8 e9 e10 e11 e12]=textread(xaxisfile,...
                '%s%f%f%s%f%f%f%f%f%f%f%s','commentstyle','matlab');
            [f1 f2 f3 f4 f5 f6 positionx positiony f9 peakODtemp sigmaxtemp sigmaytemp f13 f14 f15 f16 f17 f18 sigmaxerrortemp sigmayerrortemp f21 f22 f23 f24 f25 numatomserrortemp numatomstemp chi3Dtemp sizeparametertemp]=textread(imagefile2D,...
                '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','commentstyle','matlab');
            
            if basenamenumber==1
                numatomssum=numatomstemp;
                sigmaxsum=sigmaxtemp;
                sigmaysum=sigmaytemp;
            else
                numatomssum=numatomssum+numatomstemp;
                sigmaxsum=sigmaxsum+sigmaxtemp;
                sigmaysum=sigmaysum+sigmaytemp;
            end
        end
        
    end %end of loop through basenames
    
    numatomssum=numatomssum./numberofbasenames./numberfactor; %average number of all data files
    sigmaxsum=sigmaxsum.*sizeparameter./numberofbasenames; %average size of all data files
    sigmaysum=sigmaysum.*sizeparameter./numberofbasenames; %average size of all data files
    if singlepeakfit==1; %Fit only with single Lorentzian
        %Make initial guesses
        center1 = scancenter;
        %center1 = -24;
        amplitude1 = -(max(numatomssum)-min(numatomssum))./2;
        width1 = 0.1; %rough estimate because it doesn't seem to vary that much usually without washing out
        offset = mean(numatomssum(1:5));
        InitialGuesses=[amplitude1 center1 width1 offset];
        %xaxis(29:imagesize),numatoms(29:imagesize)
        [P r J] = nlinfit(xaxisfit(:),numatomssum(:),@LorentzFitFunction1peak, InitialGuesses);
        
        amplitude1 = P(1);
        center1 = P(2);
        width1 = P(3);
        offset= P(4);
        
        ci(:,:,1)=nlparci(P(:),r,J); %error in the fit parameters 95% confidence interval
        sdev(:,1)=(ci(:,2,1)-ci(:,1,1))/4; %turn confidence interval into 1 sigma
        %sdev is matrix w/ num. of rows equal to total num. of fit param. and num.
        %of columns equal to num. of data files
        
        pts=100;
        for t=1:pts
            indevariable(t)=min(xaxisfit)+t.*(max(xaxisfit)-min(xaxisfit))./pts;
            fitcurve(t) = amplitude1*1./(1+(2*(indevariable(t)-center1)/width1).^2)+offset;
        end
    end
    figure(998) %plot the number of atoms from the summed images
    plot(xaxisfreq,numatomssum,'.b')
    hold on
    if fitting==1 && fittingtype==2
        if bodyloss==2
            plot(TIME,N2body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
        elseif bodyloss==1
            plot(TIME,N1body,'Color',[1 0 0],'LineStyle','-', 'LineWidth',1);
        end
        hold on
        text(mean(xaxisfreq)/1.5,2*mean(numatomssum),strcat(char('\tau ='),num2str(lifetime1body(filecounter),3),'\pm',num2str(lifetime1bodyerror(filecounter),2), ' s'),'FontSize',textsize,'FontWeight','bold');
    elseif fitting==1 && fittingtype==3
        plot(indevariable,fitcurve,'Color',[1 0 0],'LineStyle','-','LineWidth',1);
        hold on
        %                 text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.8,ztext,strcat('Width=',num2str(width1),'+/-',num2str(sdev(3))),'FontSize',textsize);
        %                 text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.6,ztext,strcat('Center=',num2str(center1),'+/-',num2str(sdev(2))),'FontSize',textsize);
        %                 text(xtext,ytext+(max(numatoms)-min(numatoms)).*0.4,ztext,strcat('FracAmpl=',num2str(abs(amplitude1/offset))),'FontSize',textsize);
    end
    %xlim([0+(min(xaxis)-abs(min(xaxis))/4) max(xaxis)+max(xaxis)/4]);
    %xlim([min(xaxis) max(xaxis)]);
    %xlim([min(xaxisfreq) max(xaxisfreq)]);
    %xlim([130 155]);
    ylim([0+(min(numatomssum)-min(numatomssum)/4) max(numatomssum)+max(numatomssum)/4]);
    %ylim([0 1.5]);
    xlabel('Binding Energy [MHz]','FontSize',fontsize,'FontWeight','bold');
    ylabel('# Atoms (10^6)','FontSize',fontsize,'FontWeight','bold');
    set(gca,'FontSize',fontsize-4,'FontWeight','bold');
    title(strcat(num2str(powerW(filecounter)),' MHz Summed Images'),'FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
    hold on
    
    figure(999) %plot the temperatures from the summed images
    ysizetempsum=(sigmaysum.^2.*mass)./(kBoltz.*(delaytime./timefactor).^2)./templow;xsizetempsum=(sigmaxsum.^2.*mass)./(kBoltz.*(delaytime./timefactor).^2)./templow;
    plot(xaxisfreq,xsizetempsum,'xb')
    hold on
    plot(xaxisfreq,ysizetempsum,'ob')
    
    %ylim([0+(min(numatomssum)-min(numatomssum)/4) max(numatomssum)+max(numatomssum)/4]);
    ylim([7 18]);
    xlabel('Binding Energy [MHz]','FontSize',fontsize,'FontWeight','bold');
    ylabel('Temp [\muK]','FontSize',fontsize,'FontWeight','bold');
    set(gca,'FontSize',fontsize-4,'FontWeight','bold');
    title(strcat(num2str(powerW(filecounter)),' MHz Summed Images'),'FontSize',fontsize,'FontWeight','bold','FontAngle','italic');
    hold on
    outputdata=[xaxisfreq(:) numatomssum(:)];
    outputfilename = strcat('aggregate',num2str(scancenter));
    outputfilepath = char(strcat(deblank(directoryvector(filecounter)),deblank(outputfilename),'.txt'));
    dlmwrite(outputfilepath,outputdata,'\t');
end

status = 'done'
basehead=char('Files_20091222');