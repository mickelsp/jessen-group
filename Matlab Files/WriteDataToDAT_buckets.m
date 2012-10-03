%%%This reads the batch file and the 2Dfitparams.txt file and plots 
%%%number and temperature versus time.  Results are written to file
%%%in correct format for fitting with MMA evaporation model.
%% Initialization
close all
clear all

%%%Constants
lambda=689*10^(-9); %intercombination line wavelength
kBoltz=1.38*10^(-23);
protonmass=1.672*10^(-27); %proton mass
templow=1*10^(-6);
h=6.626*10^(-34);
w0=100*10^(-6);     %dipole trap beam waist
lambdaODT=1064*10^(-9); %dipole trap wavelength
factor=2*(h*150*10^3/kBoltz)*(50*10^(-6))^2; %To compute ODT trap depth, for 1064nm light; see MMA notebook on ODT; factor of 2 for crossed config
numberfactor=10^6; %scales number of atoms for fitting
timefactor=10^3; %scales time from ms to s
sizefactor=10^3; %scales size from m to mm5

%% File and Directory information
datestamp='20091010';
batchdir=char('/Users/work/Documents/Analysis/Unforced Evap of 84Sr/');
batchname=char('Files_unforcedevap84Sr.txt');
batchpath=strcat(batchdir,batchname);
%datadirbase='/Users/work/Documents/Analysis/Unforced Evap of Mixed Species/'; %location to save data
datadirbase=batchdir;

[directoryvector basenamevector timevector feshfreq feshbachpower droptime ODTdepth dummy01 dummy02]=...
    textread(batchpath, '%s%s%f%f%f%f%f%f%f','commentstyle','matlab'); %read in data file for atoms
numfiles=length(timevector);

%% Loop through all files
for filecounter=1:numfiles
    basehead=strcat(directoryvector(filecounter),basenamevector(filecounter));
    filename=strcat(basenamevector(filecounter));
    if ODTdepth(filecounter)==0
        color='Alone';
    else
        color='Mixed';
    end %of Fesh detuning label
    delaytime=droptime(filecounter);
    timestamp=num2str(timevector(filecounter));
    isotope=num2str(dummy02(filecounter));
    mass=dummy02(filecounter).*protonmass;
    trapdepth=feshfreq(filecounter);
    idtag=strcat(isotope,'Sr',color,num2str(trapdepth),'Vtrap',datestamp,'_',timestamp);

    %%
    clear e1 xaxis e3 e4 e5 e6 e7 e8 e9 e10 e11 e12
    clear f1 f2 f3 f4 f5 f6 f7 f8 f9 peakOD sigmax sigmay f13 f14 f15 f16 f17 f18 sigmaxerror sigmayerror f21 f22 f23 f24 f25 numatomserror numatoms chi3D sizeparameter
    %%%INPUT FILES
    xaxisfile=char(strcat(deblank(basehead),deblank('.batch')));
    imagefile2D=char(strcat(deblank(basehead),deblank('2Dfitparams.txt')));
    %xaxisfile=char(strcat(deblank(filename),deblank('.batch')));
    %imagefile2D=char(strcat(deblank(filename),deblank('2Dfitparams.txt')));
    %%%We use changed values from batch file and number of atoms from 2D fit parameters file (output of imagefit).
    [e1 xaxis e3 e4 e5 e6 e7 e8 e9 e10 e11 e12]=textread(xaxisfile,...
        '%s%f%f%s%f%f%f%f%f%f%f%s','commentstyle','matlab');
    [f1 f2 f3 f4 f5 f6 f7 f8 f9 peakOD sigmax sigmay f13 f14 f15 f16 f17 f18 sigmaxerror sigmayerror f21 f22 f23 f24 f25 numatomserror numatoms chi3D sizeparameter]=textread(imagefile2D,...
        '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','commentstyle','matlab');
    imagesize=size(xaxis,1); %number of data points
    
    if dummy02(filecounter)==87
        numatoms=numatoms.*2; %factor of two because atoms distributed between different spin states
    end

    %%
    %%%Data preparation
    numatoms = numatoms(:)./numberfactor; %number of atoms in millions of atoms
    numsystematic=0.0; %systematic error in number of atoms
    %Add systematic and statistical error in quadrature to get total error in number
    numatomserror = (numsystematic.^2+numatomserror.^2).^(0.5)./numberfactor;
    xaxis=xaxis(:);
    sigmax = abs(sigmax.*sizeparameter); %size is spit out of imagefit in # of pixels.
    sigmay = abs(sigmay.*sizeparameter); %sizeparameter converts to meters.
    sigmasystematic = 10.4*10^(-6); %one PixelFly pixel is our best resolution
    %Add systematic and statistical error in quadrature to get total error in
    %sizes
    sigmaxerror = sqrt((sigmaxerror.*sizeparameter).^2+sigmasystematic.^2);
    sigmayerror = sqrt((sigmayerror.*sizeparameter).^2+sigmasystematic.^2);
    sigmainitial=sigmax(1); %initial size of cloud; using x instead of y is arbitrary; all for lifetime fitting
    samplesizeinitial=ones(1,length(sigmax)).*sigmainitial;   %vector containing initial size of cloud

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
    %%%Calculate temperature
    ysizetemp=((sigmay).^2.*mass)./(kBoltz.*(delaytime./timefactor).^2)./templow;
    xsizetemp=((sigmax).^2.*mass)./(kBoltz.*(delaytime./timefactor).^2)./templow;
    avgtemp=(ysizetemp(:)+xsizetemp(:))./2;
    
    [bucketsx numatomsbucket numbererrorbucket]=buckets(xaxis,numatoms,210); %84Sr dataset binned
    [bucketsx xsizetempbucket xsizetemperrorbucket]=buckets(xaxis,xsizetemp,210); %84Sr dataset binned
    [bucketsx ysizetempbucket ysizetemperrorbucket]=buckets(xaxis,ysizetemp,210); %84Sr dataset binned
    avgtempbucket=(xsizetempbucket+ysizetempbucket)./2;

    %%
    %%%Plot as a double check
    figure(filecounter)
    %plot(xaxis,numatoms,'ok',xaxis,xsizetemp,'xr',xaxis,ysizetemp,'sb',xaxis,avgtemp,'.c');
    hold on
    plot(bucketsx,numatomsbucket,'ok',bucketsx,avgtempbucket,'xr');
    %ylim([0 10]);
    legend('N','T_x','T_y','Location','Best');
    set(gca,'FontSize',20,'FontWeight','bold');
    title(filename,'FontSize',10);

    %%
    %%%Write data to .dat file
    xaxis = bucketsx';
    numatoms = numatomsbucket';
    ysizetemp = ysizetempbucket';
    xsizetemp = xsizetempbucket';
    avgtemp = avgtempbucket';

    %%%Build path name; make directory if it doesn't already exist.
    datadir=char(directoryvector(filecounter));
    
    if isdir(datadir)==0
        mkdir(datadir);
    end

    atomsname=strcat(datadir,'/numAtoms_',idtag,'.dat');
    xsizename=strcat(datadir,'/xSizeTemp_',idtag,'.dat');
    ysizename=strcat(datadir,'/ySizeTemp_',idtag,'.dat');
    avgsizename=strcat(datadir,'/avgSizeTemp_',idtag,'.dat');

    dlmwrite(atomsname,[xaxis/10^3,numatoms*10^(6)], 'delimiter', ' ', 'newline', 'pc');
    dlmwrite(xsizename,[xaxis/10^3,xsizetemp*10^(-6)], 'delimiter', ' ', 'newline', 'pc');
    dlmwrite(ysizename,[xaxis/10^3,ysizetemp*10^(-6)], 'delimiter', ' ', 'newline', 'pc');
    dlmwrite(avgsizename,[xaxis/10^3,avgtemp*10^(-6)], 'delimiter', ' ', 'newline', 'pc');
end %of loop through all files

status = 'done'