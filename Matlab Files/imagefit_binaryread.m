%This program is designed to read several datafiles and corresponding background
%files,  plots these normalized data sets on the same graph.
%NOTES:  The input file should contain a column of data filenames, a column of
%background filenames, and a column of delay times.
%atomic physics see page xxx of notebook
close all
clear all
%clc

%% Initialization and Options
alpha=40;%the of counts per photon
lambda=461*10^(-9);
detuning=-0.72*10^6; % image beam detuning in Hz
naturalwidth=32*10^6; %ion fwhm in Hz
textsize=16;
axisfontsize=14;

%All frequencies here ONLY go into alldata(k,1); frequencies.
dzeeman=-251.43;    %constant overnight
dPAS=-250;          %constant overnight
dsatabs=-41;        %constant overnight
dshift=0;           % 0 1st run, 260.5 for 2nd run
catseyesign=0;     %depends on peak
%dcat is defined above alldata !!!!!!!!!!!!

%SPECIFICS OF FITTING ROUTINES
printlinewidths=0;%1 for printing linewidth plot, 0 for not
printspectrum=1; %1 for printing spectrum plot, 0 for not
smoothplotsize=4; %number of points to average for plotting
sampleinterleave=1; %for big matrices form a matrix using only one point for every sampleinterleave points in a row/column, Always=1 for final analysis and error calc
pixelsize=9.2*10^(-6);%7.93*10^(-6);%0.00000793; %size of a pixel in meters; see p. 71w of Lab Notebook H for latest calibration (size prior to 10.10 (?) was 10.4 microns/pixel; between 10.10 and now is 9.2 microns/pixel
pixelconv=pixelsize*10^6;%converts pixels to microns
guesssmoothsize=2; %This variable is used in fitting, and also smoothes out the linear density cuts 

%ROI never bigger than 1280x1024
fullroi=1;
if fullroi==1
    roirowstart=1;      %Starting row for region of interest (min = 1)
    roirowstop=1280;    %Ending row for region of interest (max = 1280)
    roicolstart=1;      %Starting column for region of interest (min = 1)
    roicolstop=1024;    %Ending column for region of interest (max = 1024)
elseif fullroi==0
    roicolstart=390;      %Starting column for region of interest (min = 1)
    roicolstop=540;    %Ending column for region of interest (max = 1024)
    roirowstart=800;      %Starting row for region of interest (min = 1) 600
    roirowstop=950;    %Ending row for region of interest (max = 1280)
    
    %%%For variable region of interest; make sure variableROI flag is 1
    rowstartvector=[655 655 655 655 665 665 675 675 685 695 700];
    rowstopvector=[755 755 755 755 745 745 735 735 725 725 720];
    colstartvector=[300 300 300 300 305 305 315 315 325 325 340];
    colstopvector=[400 400 400 400 390 390 380 380 370 370 360];
end
absorptionsign=-1; %equals -1 or +1, respectively, depending on whether atoms are in the "atoms" image or the "background" image; normal use of imagefit has -1 value.
variableROI=0;  %1 to vary ROI dynamically for falling cloud
softwarebinsize=8; %Number of pixels square to bin (i.e. area of bin, in # of pixels = softwarebinsize*softwarebinsize
CCDbinning=1; %number of pixels binned when first recording data
prebinning=softwarebinsize*CCDbinning*sampleinterleave;
cutborders=0; %number of column and rows to trim on either side of the data matrix; same as roi parameters, 'cept inverse way of doing it.
limit=5; %number of column and rows at the corners to find the uncertainties/noise
scalebaselevel=0;% 1 to scale the base level of atom and background image, 0 for not
numberforscale=5;% number of column and rows to scale the base level of atom and background images
r1=0;       %r1 and r2 used for annular ring stuff in plasma experiment.
r2=50;
matrixsize=[1280 1024]; %matrix size for binary reading.

singlefit=0;        %0 to plot entire basename list, 1 to plot just the specified basename
plotfits=0;         %1 to plot fits/data, 0 to not plot fits; 2D Gaussian
plotevolution=0;    %1 to show surface plots for each data point on one plot
plotdensity=0;      %1 to show plots of density as function of x position (slices along y); helps to spot BEC
standalone=1;       %1 to plot data by itself (not fits, residuals, etc...); 0 to produce normal plots
fit=1;              %1 to fit, 0 just to plot; 2D Gaussian
spectrumfit=0;      %1 to fit spectrum, 0 to not fit; 1D
fittingfunction=0;  %2 for Lorentzian, 1 for Voigt, 0 for Gaussian
binaryread=1;       %0 if the files are ASCII; 1 if the files are binary
voltagetodetuning=0; %1 if we want detuning, 0 if we want voltage
save2places=0;      %1 means we save output files to both primary and secondary directories for backup purposes

vtodconversion = 6.4; %this factor sets the linear conversion rate from voltage to detuning when plotting atom size against detuning using voltage data

%% Directory and Files
%singlebasename=('PASscan1940MHz');
directory=char('/Users/Work/Documents/Analysis/Image Beam Frequency/');
basehead=char('Files_imaging');

if singlefit==0;
basenamelist=char(strcat(deblank(directory),deblank(basehead),deblank('.txt')))

[directoryvector basenamevector timevector dummy01 dummy02 dummy03 dummy04 dummy05 dummy06]=...
    textread(basenamelist, '%s%s%f%f%f%f%f%f%f','commentstyle','matlab'); %read in data file for atoms
temp=size(basenamevector);
numberofbasenames=temp(1);%max of ten for figure 200 to work
end

if singlefit==1;
    numberofbasenames=1;
end

%% Loop Through Basenames
filecounter=0   %this will keep track of all the files analyzed in all the batches
for basenamenumber=1:numberofbasenames
    filecounter=filecounter+1
    
    clear basename;
    clear file ;
    clear imagevco;
    clear vcovoltage;
    clear motdetuning;
    clear d;
    clear delay;
    clear smoothopticaldepthimagesingle;
    clear cutimageatom;
    clear cutimageback;
    clear sizes;
    
    if singlefit==0;
        basename=char(basenamevector(basenamenumber))
        batchfile=char(strcat(deblank(directoryvector(basenamenumber)),deblank(basename),deblank('.batch')));    %array of batchfile names

        [file imagevco vcovoltage motdetuning d delay dummy01 dummy02 dummy03 zeemanpower samplehold wavemeter]=textread(batchfile, '%q%f%f%s%f%f%f%f%f%f%f%s','commentstyle','matlab'); %read in all files, no limit
        temp=size(imagevco);
      
        basenamesize=temp(1);%automatically determine how many files in batch associated with this basename
        q=[basenamesize]% array of number of data files in each batch

        batchdelay(1:q,basenamenumber)=delay(1:q); %make array of all delays
        batchvcovoltage(1:q,basenamenumber)=vcovoltage(1:q); %make array of all vco voltages
        batchimagevco(1:q,basenamenumber)=imagevco(1:q); %make array of all image vco voltages

    elseif singlefit==1
        file=char(strcat(deblank(directory),deblank(singlebasename),deblank('.dat')));
        dataatom=char('atoms.bny');
        databack=char('back.bny');
    end

    %OUTPUT FILES
    outputdatafile=char(strcat(deblank(directoryvector(basenamenumber)),deblank(basename),deblank('2Dfitparams.txt')));
    %array of 2D fit parameters
    %(freq,voltage,odintegralerror,odintegral,etc...)
    %deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2)),
    outputfilesfile=char(strcat(deblank(directoryvector(basenamenumber)),deblank(basename),deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2)),deblank('.out')));
    %array of filenames from input file; currently disabled: see below
    outputfiledatamatrix=(strcat(deblank(directoryvector(basenamenumber)),deblank(basename),deblank('2Dfitparams.xls')));
    %array of data similar to outputdatafile (opticaldepths,freqs,odsdev,chi3D,filenumbers,residue1D,#atoms)
    %,deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2))
    fitdata1D=(strcat(deblank(directoryvector(basenamenumber)),deblank(basename),deblank('1Dfitparams.txt')));
    %array of parameters from fitting routine(delay,density,Nions,amplitude,amplitudedev,centerfreq,centerfreq(dev),linewidth,linewidth(dev),bg,bgdev,bgslope,bgslopedev)
    %,deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2))
    basenamespectrumplot=(strcat(deblank(directoryvector(basenamenumber)),deblank(basename),deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2)),deblank('spec.jpg')));
    %jpg file of spectrum and fit
    basenamespectrumstatsplot=(strcat(deblank(directoryvector(basenamenumber)),deblank(basename),deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2)),deblank('specstats.jpg')));
    %jpg file of residuals, chisquared, etc...
    linewidthvsdelayplot=(strcat(deblank(directoryvector(basenamenumber)),deblank(basehead),deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2)),deblank('linewidths.jpg')));
    %jpg file of linewidths plot
    if save2places==1
        outputdatafile2=char(strcat(deblank(directory2),deblank(basename),deblank('2Dfitparams.txt')));
        outputfilesfile2=char(strcat(deblank(directory2),deblank(basename),deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2)),deblank('.out')));
        outputfiledatamatrix2=(strcat(deblank(directory2),deblank(basename),deblank('2Dfitparams.xls')));
        fitdata1D2=(strcat(deblank(directory2),deblank(basename),deblank('1Dfitparams.txt')));
        basenamespectrumplot2=(strcat(deblank(directory2),deblank(basename),deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2)),deblank('spec.jpg')));
        basenamespectrumstatsplot2=(strcat(deblank(directory2),deblank(basename),deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2)),deblank('specstats.jpg')));
        linewidthvsdelayplot2=(strcat(deblank(directory2),deblank(basehead),deblank('-'),deblank(num2str(r1)),deblank('-'),deblank(num2str(r2)),deblank('linewidths.jpg')));
    end
    if binaryread==1
        dataatom=char('atoms.bny');
        databack=char('back.bny');
    elseif binaryread==0
        dataatom=char('atoms.dat');
        databack=char('back.dat');     
    end

    clear opticaldepthsfrom2Dfits
    clear alldata
    clear odsdev
    clear sdev
    clear chi3D
    clear sigx3D
    clear sigy3D
    clear residue1D
    clear filenumbers
    clear gradmatrix
    clear nogradoddata
    clear nogradodmatrix
    sumsigx=0;
    sumsigy=0;
    avgsigx=0;
    avgsigy=0;
    avgsigz=0;
    Nions=0;
    % the columns of this matrix are the OD, OD error, chi square 3D, 1D
    % residue no errors

    %this loop processes all the data files in this batch
    counter=0;
    for k=1:q
        s=char(strcat(deblank(directoryvector(basenamenumber)),deblank(file(k)),deblank(dataatom(1,:))));
        t=char(strcat(deblank(directoryvector(basenamenumber)),deblank(file(k)),deblank(databack(1,:))));
        counter=counter+1
        if binaryread==1;
            %Read in binary file created in LabView
            fid=fopen(s,'rb','ieee-be'); %Open atoms image
            %LabView saves binary in "big endian" ('be') format: most significant bit in
            %lowest memory address. Matlab needs this info to import binary file correctly.
            fullrawimageatom=fread(fid,matrixsize,'*int16');
            fclose(fid);
            %LabView and Matlab treat matrix coordinates differently.
            %(0,0) for LabView is lower left corner; for Matlab is upper left corner.
            %Therefore, flip about horizontal axis.
            %fullrawimageatom=flipud(fullrawimageatom); %I think this makes it
            %opposite LabView. 2007/01/15

            fid=fopen(t,'rb','ieee-be'); %Open background image
            fullrawimageback=fread(fid,matrixsize,'*int16');
            fclose(fid);
            %fullrawimageback=flipud(fullrawimageback); %I think this makes it
            %opposite LabView. 2007/01/15

        elseif binaryread==0;
            fullrawimageatominter=dlmread(s); %read in data file for atoms
            fullrawimagebackinter=dlmread(t); %read in data file for back
            fullrawimageatom=fullrawimageatominter';
            fullrawimageback=fullrawimagebackinter';
        end
        
            bgsize=size(fullrawimageback);
            fullsizefull=size(fullrawimageatom);    %Size of atoms matrix
            fullcolumnfull=fullsizefull(1);
            fullrowfull=fullsizefull(2); %no. of rows 
        if variableROI==1
            roirowstart=rowstartvector(k);
            roirowstop=rowstopvector(k);
            roicolstart=colstartvector(k);
            roicolstop=colstopvector(k);
        end
        %%%%Region of interest; done before any other array modification.
        %%%%Define in terms of original coordinates (full size array = 1280x1024)
        roiimageatom=fullrawimageatom(roirowstart:roirowstop,roicolstart:roicolstop);
        roiimageback=fullrawimageback(roirowstart:roirowstop,roicolstart:roicolstop);
        roisizes=size(roiimageatom);
        roicol=roisizes(1);
        roirow=roisizes(2);
        
        %%%%true software binning
        if softwarebinsize > 1
            for i=1:roicol/softwarebinsize
                for j=1:roirow/softwarebinsize
                    rawimageatom(i,j)=sum(sum(roiimageatom(i*softwarebinsize-(softwarebinsize-1):i*softwarebinsize,j*softwarebinsize-(softwarebinsize-1):j*softwarebinsize)));
                    rawimageback(i,j)=sum(sum(roiimageback(i*softwarebinsize-(softwarebinsize-1):i*softwarebinsize,j*softwarebinsize-(softwarebinsize-1):j*softwarebinsize)));
                end
            end
        else
            rawimageatom = roiimageatom;
            rawimageback = roiimageback;
        end

        sizefull=size(rawimageatom);    %Size of atoms matrix post-binning
        columnfull=sizefull(1);
        rowfull=sizefull(2);

        if scalebaselevel==1
            sumylowxlow=sum(sum(rawimageatom(1:numberforscale,1:numberforscale)));;
            sumylowxhigh=sum(sum(rawimageatom(1:numberforscale,rowfull-numberforscale:rowfull)));
            sumyhighxlow=sum(sum(rawimageatom(columnfull-numberforscale:columnfull,1:numberforscale)));
            sumyhighxhigh=sum(sum(rawimageatom(columnfull-numberforscale:columnfull,rowfull-numberforscale:rowfull)));
            baselevelatom=sumylowxlow+sumylowxhigh+sumyhighxlow+sumyhighxhigh;
            averageylowxlow=sum(sum(rawimageback(1:numberforscale,1:numberforscale)));
            sumylowxhigh=sum(sum(rawimageback(1:numberforscale,rowfull-numberforscale:rowfull)));
            sumyhighxlow=sum(sum(rawimageback(columnfull-numberforscale:columnfull,1:numberforscale)));
            sumyhighxhigh=sum(sum(rawimageback(columnfull-numberforscale:columnfull,rowfull-numberforscale:rowfull)));
            baselevelback=sumylowxlow+sumylowxhigh+sumyhighxlow+sumyhighxhigh;
            rawimageatom=rawimageatom*baselevelback/baselevelatom;
        end %of scalebaselevel

        %%%%Form matrix from only part of the data
        %%%%Sample interleave section; removes every sampleinterleave row or column
        for i=1:columnfull/sampleinterleave
            for j=1:rowfull/sampleinterleave
                for v=1:sampleinterleave
                    tempatom(v)=sum(rawimageatom(((i-1)*sampleinterleave+1):i*sampleinterleave,((j-1)*sampleinterleave+v)));
                    tempback(v)=sum(rawimageback(((i-1)*sampleinterleave+1):i*sampleinterleave,((j-1)*sampleinterleave+v)));
                end
                rawimageatom2(i,j)=sum(tempatom);
                rawimageback2(i,j)=sum(tempback);
            end
        end
        rawimageatom=rawimageatom2; clear rawimageatom2 %cleaning up array
        rawimageback=rawimageback2; clear rawimageback2 %cleaning up array

        sizefull=size(rawimageback);
        columnfull=sizefull(1); %no. of cols after sampleinterleave
        rowfull=sizefull(2); %no. of rows after sampleinterleave
        for i=1+cutborders:columnfull-cutborders,
            for j=1+cutborders:rowfull-cutborders,
                cutimageback(i-cutborders,j-cutborders)=rawimageback(i,j);
            end
        end

        sizefull=size(rawimageatom);
        columnfull=sizefull(1);
        rowfull=sizefull(2); %no. of rows

        %%%%Cut edges of data by cutborders
        for i=1+cutborders:columnfull-cutborders,
            for j=1+cutborders:rowfull-cutborders,
                cutimageatom(i-cutborders,j-cutborders)=rawimageatom(i,j);
            end
        end
        sizes=size(cutimageatom);
        column=sizes(1);
        row=sizes(2); %no. of rows
        S=zeros(column,row,q);

        %cutimageback=cutimageback(1:column,1:row);
        %sigma=sqrt(cutimageback(:,:));%expected noise, the factor of two acounts for the shot noise of the atom signal and background
        %sigma2=cutimageback(:,:);   ck
        cutimageatom=abs(cutimageatom);
        cutimageback=abs(cutimageback);
        opticaldepthimagesingle=(log(cutimageatom)-log(cutimageback))*absorptionsign;%makes absorption positive
        opticaldepthimagebatch(1:column,1:row,k)=opticaldepthimagesingle; %array of log(data) -log(back)
        
%         figure(777)
%         surf(cutimageatom)
%         hold on
%         set(gca,'View',[0 90])
%         shading flat
%         grid off
%         axis off
%         xlabel('x');
%         ylabel('y');
%         hold off
       
        %Fit Data
        %Generate guesses and bounds
        %smoothed optical density for finding guesses
        Zguess=filter2(ones(guesssmoothsize,guesssmoothsize),opticaldepthimagesingle(:,:))/guesssmoothsize^2; %atom-back, smoothed lots for finding guesses
        amplitude=max(max(Zguess));%%this will give trouble if slope is too big
        %guesses for x and y position of cloud; bin by 4 and choose max
        %bin; multiply by 4 to get coordinates in full size array
        if softwarebinsize > 5
            centerfindbinsize=1; %bin by 1 is sufficient in this case because binning was already pretty high.
        elseif softwarebinsize <= 5
            centerfindbinsize=4; %bin by at least 4 because binning was low to begin with
        end
        for i=1:size(Zguess,1)/centerfindbinsize
            for j=1:size(Zguess,2)/centerfindbinsize
                Zguessbinned(i,j)=sum(sum(Zguess(i*centerfindbinsize-(centerfindbinsize-1):i*centerfindbinsize,j*centerfindbinsize-(centerfindbinsize-1):j*centerfindbinsize)));
            end
        end
        [C,I] = max(Zguessbinned); [C2,I2]=max(max(Zguessbinned));
        xoffset=I2*centerfindbinsize;
        yoffset=I(I2)*centerfindbinsize;
        sigx=row/10;
        sigy=column/10;%%%%ONE OF THESE SHOULD BE COLUMN
        averageylowxlow=sum(sum(Zguess(1:5,1:5)))/25;
        averageylowxhigh=sum(sum(Zguess(1:5,row-5:row)))/25;
        averageyhighxlow=sum(sum(Zguess(column-5:column,1:5)))/25;
        averageyhighxhigh=sum(sum(Zguess(column-5:column,row-5:row)))/25;
        offset=(averageylowxlow+averageylowxhigh+averageyhighxlow+averageyhighxhigh)/4 ;
        slopex=((averageylowxhigh+averageyhighxhigh)-(averageylowxlow+averageyhighxlow))/(2*row);
        slopey=((averageyhighxlow+averageyhighxhigh)-(averageylowxlow+averageylowxhigh))/(2*column);

        InitialConditions=[amplitude, sigx, sigy, xoffset, yoffset,offset,slopex,slopey];%
        %pause
        lowerbound=[0, 1, 1, 10, 10, -3, 0, 0,-.01,-.01];%
        upperbound=[1, 50, 50, 90, 90, 3, 1, 1,.01,.01];%
        temp=size(InitialConditions);
        numberofparameters=temp(2);%number of 2d fit parameters

        %Statistics of Image
        %Weird things happen when you calculate these with smoothed images!!!
        %Get estimate of noise, you need to subtract image from background because of fringes and variations in intensity
        average1=sum(sum(cutimageatom(1:limit,1:limit)-cutimageback(1:limit,1:limit)))/limit^2;
        rms1=sqrt(sum(sum((cutimageatom(1:limit,1:limit)-cutimageback(1:limit,1:limit)).^2))/limit^2);
        rmsdeviation1=sqrt(rms1^2-average1^2);

        average2=sum(sum(cutimageatom(1:limit,row-limit:row)-cutimageback(1:limit,row-limit:row)))/limit^2;
        rms2=sqrt(sum(sum((cutimageatom(1:limit,row-limit:row)-cutimageback(1:limit,row-limit:row)).^2))/limit^2);
        rmsdeviation2=sqrt(rms2^2-average2^2);

        average3=sum(sum(cutimageatom(column-limit:column,1:limit)-cutimageback(column-limit:column,1:limit)))/limit^2;
        rms3=sqrt(sum(sum((cutimageatom(column-limit:column,1:limit)-cutimageback(column-limit:column,1:limit)).^2))/limit^2);
        rmsdeviation3=sqrt(rms3^2-average3^2);

        average4=sum(sum(cutimageatom(column-limit:column,row-limit:row)-cutimageback(column-limit:column,row-limit:row)))/limit^2;
        rms4=sqrt(sum(sum((cutimageatom(column-limit:column,row-limit:row)-cutimageback(column-limit:column,row-limit:row)).^2))/limit^2);
        rmsdeviation4=sqrt(rms4^2-average4^2);

        error=abs((rmsdeviation1+rmsdeviation2+rmsdeviation3+rmsdeviation4)/(4))%
        %FIX IMAGINARY ERRORS....THIS IS A PROBLEM!!!!!

        %Create fitting matrices
        %(sum(sum(cutimageback))/(row*column)) is average value of pixel in cutimageback
        errorod(1:row*column)=error/(sum(sum(cutimageback))/(row*column)); %normalized error - relative to average pixel value
        xvec=1:row;
        yvec=1:column;%%%%ONE OF THESE SHOULD BE COLUMN
        [xdata,ydata]=meshgrid(xvec,yvec);

        numberofpoints=column*row*ones(column,row);%matrix in which each entry is the number of points
        coordmatrix=[xdata(:),ydata(:),errorod(:),numberofpoints(:)];

        if fit==1 && column*row<100000; % skip fit if fit==0 or matrix is larger than 100,000 pixels (very slow otherwise)
            %actual fit
            Z2=(opticaldepthimagesingle(:)./(errorod(:)*sqrt(column*row)));  %flattened(1D) list of image weighted by uncertainties, not smoothed

            [P,r,J]=nlinfit(coordmatrix,Z2(:),@gaussfitfunctionerrorsslope,InitialConditions); %fitting statement - expecting data divided by weight and sqrt(number of points)

            ci=nlparci(P,r,J); %error in the fit parameters 95% confidence interval
            %%saveci=zeros(1,6,2);
            saveci(k,:,:)=ci;
            sdev=(saveci(:,:,2)-saveci(:,:,1))/4;%turn confidence interval into 1 sigma
            %%%%%%THIS MUST BE RECORDED SOMEWHERE< ALONG WITH ALL FIT PARAMETERS
            %odsdev=sdev(:,1);
            %second derivative of chisquared wrt amplitude at fit paramters
            %open 'Dgaussfitfunctionrf.m'%derivative of chi squared, not yet summed over all points
            %presumdervchi3D=Dgaussfitfunctionrf(P, coordmatrix); % this is the second dervative of chi3d with respect to amplitude, not yet summed over all points
            %dervchi3(counter)=sum(presumdervchi3D);
            %odsdev2(counter)= sqrt(2)*sqrt(1/((column*row)*dervchi3(counter)));
            %error=.01
            %temp=size(smoothopticaldepthimagesingle)
            %Z3 = smoothopticaldepthimagesingle(:)./errorod/8100

            %Below are the calculations of residues and chi sqares considering and not considering errors for the 3D fits
            weightedresidue3D=(Z2(:)-gaussfitfunctionerrorsslope(P, coordmatrix)); %weighted residues
            chi3D(k)=(transpose(weightedresidue3D)*weightedresidue3D);
            %weightedresidue3Dnoerrors=(opticaldepthimagesingle(:)-gaussfitfunctionnoerrorsslope(P, coordmatrix));
            %chi3Dnoerrors(counter)=(transpose(weightedresidue3Dnoerrors)*weightedresidue3Dnoerrors)*(1/(column*row))*(1/errorod(1))^2

            %residue=cumsum(((smoothopticaldepthimagesingle(:)-gaussfitfunctionrf(P, coordmatrix))./sigma(:)).^2);
            %sqrt(residue(8100))/(8100)  %chisquared
            %form output parameters
            %fraction absorption is p(1)

            %removing slope and background
            gradmatrix=(xdata(:)-P(4))*P(7)+(ydata(:)-P(5))*P(8)+P(6);
            nogradoddata=opticaldepthimagesingle(:)-gradmatrix(:);
            nogradodmatrix = reshape(nogradoddata, size(xdata));

            pixelcount=0;
            odsum=0;
            odaverage=0;
            totalod=0;
            odintegral=0;
            %Sum the optical depth to get total optical depth
            for i=1:row
                for j=1:column
                    tempindex=(i-1)*column+j;
                    totalod=totalod+nogradoddata(tempindex);
                end
            end
            odintegralerror=errorod(1)*sqrt(row*column)*(prebinning*pixelsize)^2;
            %pixelcount
            %odaverage=odsum/pixelcount
            odintegral=totalod*(prebinning*pixelsize)^2;

            %odsdev(k,1)=error/(sum(sum(cutimageback))/(row*column))/sqrt(pixelcount);
            odsdev(k,1)=odintegralerror;

            P; %fit parameters
            onesigma = (ci(:,2)-ci(:,1))./4; %1-sigma "error" for fit parameters
            output(1:1,k,basenamenumber)=totalod;%total od
            output(2:2,k,basenamenumber)=P(1);%od average for selected region (peak od?)
            output(3:3,k,basenamenumber)=pixelsize*prebinning*P(4);%positionx in meters
            output(4:4,k,basenamenumber)=pixelsize*prebinning*P(5);%positiony in meters
            sigx3D(k)=pixelsize*prebinning*abs(P(2));%sigmax in meters
            sigy3D(k)=pixelsize*prebinning*abs(P(3));%sigmay in meters
            sigxerror3D(k)=pixelsize*prebinning*abs(sdev(k,2));%sigmaxerror in meters
            sigyerror3D(k)=pixelsize*prebinning*abs(sdev(k,3));%sigmayerror in meters
            sumsigx=sumsigx+pixelsize*prebinning*abs(P(2));
            sumsigy=sumsigy+pixelsize*prebinning*abs(P(3));
            output(5:5,k,basenamenumber)= odintegral;%absorption fraction, or more accurately - optical depth
            
            for j=1:numberofparameters
                output(5+j,k,basenamenumber)=P(j);%write fit information into output
            end
            for jj=1:numberofparameters
                output(5+numberofparameters+jj,k,basenamenumber)=(ci(jj,2)-ci(jj,1))./4; %95% ci for all fit parameters will go into output file
            end
            output(1:1,k,basenamenumber);
            output(2:2,k,basenamenumber);
            output(3:3,k,basenamenumber);
            output(4:4,k,basenamenumber);
            output(5:5,k,basenamenumber);
            %output(6:(6+2*numberofparameters),k,basenamenumber);
%             for j=1:numberofparameters
%                 lowererror(j,k,basenamenumber)=ci(j);
%             end
%             for j=numberofparameters+1:2*numberofparameters
%                 uppererror(j,k,basenamenumber)=ci(j); %this seems wrong
%             end
            %for j=1:5

            %    output(13+j,k,basenamenumber)=P(j);%write fit information into output

            %      lowererror(j,k,basenamenumber)=ci(j);
            %  end
        elseif column*row >= 100000
            error('analysis stopped because array is larger than 100,000 pixels')
        end% ends the if fit==0/1 , skip fit if fit==0
        %%%%%%%%WHY NO TEXT FOR K>6 It seems if it
        %has a figure number over 7 => no text!????????
        b=ones(smoothplotsize,smoothplotsize);
        smoothopticaldepthimagesingle=filter2(b,opticaldepthimagesingle(:,:))/smoothplotsize^2; %atom-back, smoothed a bit for fitting
        if plotevolution==1
            figure(2222)
            subplot(2,2,k)
            surf(smoothopticaldepthimagesingle)
            hold on
            set(gca,'View',[0 90])
            %set(gca,'View',[0 0])
            shading flat
            grid off
            %fontsize=12;
           % ticks=[0 50 100 150]; ticklabelstemp=ticks.*prebinning.*pixelsize.*1000;
            %ticklabels=strcat(num2str(ticklabelstemp(1)),'|',num2str(ticklabelstemp(2)),'|',...
            %    num2str(ticklabelstemp(3)),'|',num2str(ticklabelstemp(4)));
           % set(gca,'XTick',ticks,'XTickLabel',ticklabels,'FontSize',fontsize,'FontWeight','bold')
            %set(gca,'YTick',ticks,'YTickLabel',ticklabels,'FontSize',fontsize,'FontWeight','bold')
            %if k==3
            %    xlabel('x (mm)','FontSize',fontsize,'FontWeight','bold');
            %end
            %ylabel('y (mm)','FontSize',fontsize,'FontWeight','bold');
            %timelabel=[3 4.5 6];
            %text(130,75,strcat(num2str(timelabel(k)),' ms'),'FontSize',fontsize,'FontWeight','bold');
            %xlim([0 row]);
            %ylim([0 column]);
            title(strcat(num2str(imagevco(k)),'ms'));
            %hold off
        end
        if plotfits==1;
            %%%%%smooth data forplotting
            b=ones(smoothplotsize,smoothplotsize);
            %opticaldepthimagesingle is optical depth log(cutimageatom/cutimageback)
            smoothopticaldepthimagesingle=filter2(b,opticaldepthimagesingle(:,:))/smoothplotsize^2; %atom-back, smoothed a bit for fitting
            
            if standalone==1
                figure(400+k) %Plot raw data (atoms - back)
                surf(smoothopticaldepthimagesingle)
                hold on
                set(gca,'View',[0 90])
                shading flat
                grid off
                axis off
                %colormap('Hot')
                xlabel('x');
                ylabel('y');
                title(file(k));
                hold off
                
                %%%Sarah's thesis plot - absorption image - appropriate formatting
%                 figure(500+k) %Plot raw data (atoms - back)
%                 surf(smoothopticaldepthimagesingle)
%                 hold on
%                 set(gca,'View',[0 90],'FontSize',10,'FontWeight','bold')
%                 shading flat
%                 colorbar('FontSize',10,'FontWeight','bold')
%                 %grid off
%                 %axis off
%                 %xlabel('x');
%                 %ylabel('y');
%                 xlim([0 row]);
%                 ylim([0 column]);
%                 set(gca,'XTickLabel',{num2str(0);num2str(1.66);num2str(3.33);num2str(4.99);num2str(6.66);num2str(8.32);num2str(9.98)})
%                 set(gca,'YTickLabel',{num2str(0);num2str(1.66);num2str(3.33);num2str(4.99);num2str(6.66);num2str(8.32);num2str(9.98);num2str(11.65);num2str(13.31)})
%                 xlabel('Size (mm)','FontSize',14,'FontWeight','bold');
%                 ylabel('Size (mm)','FontSize',14,'FontWeight','bold');
%                 zlabel('Optical Depth','FontSize',14,'FontWeight','bold');
%                 title('Typical Absorption Image','FontSize',14,'FontWeight','bold');
%                 xxloc=128;yyloc=-4;
%                 text(xxloc,yyloc,char('Optical'),'FontSize',14,'FontWeight','bold');
%                 text(xxloc,yyloc-7,char('Depth'),'FontSize',14,'FontWeight','bold');
%                % title(file(k));
%                 hold off
                
%                 figure(888) %FIGURE for Tom's talk at Laser Science
%                 Conference (for triplet P zero temperature evolution
%                 from 2007/08/01 data)
%                 subplot(3,3,k*3-2)
%                 surf(smoothopticaldepthimagesingle)
%                 hold on
%                 set(gca,'View',[0 90])
%                 shading flat
%                 grid off
%                 %axis off
%                 fontsize=12;
%                 ticks=[0 50 100 150]; ticklabelstemp=ticks.*prebinning.*pixelsize.*1000;
%                 ticklabels=strcat(num2str(ticklabelstemp(1)),'|',num2str(ticklabelstemp(2)),'|',...
%                     num2str(ticklabelstemp(3)),'|',num2str(ticklabelstemp(4)));
%                 set(gca,'XTick',ticks,'XTickLabel',ticklabels,'FontSize',fontsize,'FontWeight','bold')
%                 set(gca,'YTick',ticks,'YTickLabel',ticklabels,'FontSize',fontsize,'FontWeight','bold')
%                 if k==3
%                     xlabel('x (mm)','FontSize',fontsize,'FontWeight','bold');
%                 end
%                 ylabel('y (mm)','FontSize',fontsize,'FontWeight','bold');
%                 timelabel=[3 4.5 6];
%                 text(130,75,strcat(num2str(timelabel(k)),' ms'),'FontSize',fontsize,'FontWeight','bold');
%                 xlim([0 row]);
%                 ylim([0 column]);
%                % title(file(k));
%                 hold off
%                 %%%ROI was col:  525 650  row:  850 950

            else
                figure(k) % Create a new figure ...
                %plot smoothed data, even though we fit unsmoothed
                legend('name(k,8:namelength)')
                subplot(2, 2, 1) % plotdata
                surf(smoothopticaldepthimagesingle, 'Marker', '.'); % Plot data
                %imagesc(smoothopticaldepthimagesingle, [0,1]); % Plot data
                %%fullrawimageatom-fullrawimageback
                %colormap(summer);
                %view([0 90])
                %cax=caxis ;
                colorbar
                shading flat
                axis off;
                grid off;
                xlabel('x')
                ylabel('y')
                title(file(k))
                legend('',num2str(batchvcovoltage(k,basenamenumber)),0)
            end %end of standalone case statement

            %USED THIS FOR ODT EVOLUTION PLOT: MAY DELETE LATER
%             if k<=4
%                 figure(999)
%                 subplot(2,2,k)
%                 surf(smoothopticaldepthimagesingle)
%                 hold on
%                 set(gca,'View',[0 90])
%                 shading flat
%                 grid off
%                 xlim([0 row]);
%                 ylim([0 column]);
%                 if k==1
%                 title('1ms');
%                 elseif k==2
%                     title('50ms');
%                 elseif k==3
%                     title('100ms');
%                 elseif k==4
%                     title('150ms');
%                 end
%             end

            if fit==1; % skip plot of fit if fit==0
                figure(500+k) %Plot raw data (atoms - back)
                subplot(2,2,1)
                surf(smoothopticaldepthimagesingle)
                hold on
                set(gca,'View',[0 90])
                shading flat
                grid off
                axis off
                title(basename,'FontSize',axisfontsize)
                xlabel('x');
                ylabel('y');
                xlim([0 row]);
                ylim([0 column]);
                % title(file(k));
                hold off

                %%%%%%%%%PLACEMENT OF CONSTANTS FOR RESIDUAL PLOTS
                xresidualconstants=0;
                yresidualconstants=column*pixelconv*2;
                zresidualconstants=1.5*(P(1));
                %%%%%%%%%%%%%%
                subplot(2, 2, 2) % plot fit

                %%%Generate plot of optical density based on fitted parameters, remove distortion of data due to dividing by weight and sqrt(number of points)
                errorod(1:row*column)=1; %
                numberofpoints=ones(column,row);%
                coordmatrix=[xdata(:),ydata(:),errorod(:),numberofpoints(:)];
                opticaldepthcalcvector=gaussfitfunctionerrorsslope(P, coordmatrix);%optical depth calculated from fit parameters (NOT WEIGHTED)
                opticaldepthcalcmatrix = reshape(opticaldepthcalcvector, size(xdata));

                surf(pixelconv*xvec,pixelconv*yvec,opticaldepthcalcmatrix, 'Marker', '.'); % Plot fit
                zlabel('Optical Depth','FontSize',axisfontsize);
                xlabel('X (microns)','FontSize',axisfontsize);
                ylabel('Y (microns)','FontSize',axisfontsize);
                set(gca,'FontSize', axisfontsize);
                text(xresidualconstants,yresidualconstants,zresidualconstants, strcat(char('Peak OD ='), num2str(P(1)), ' \pm ', num2str(sdev(k,1))),'FontSize',textsize);
                text(xresidualconstants,yresidualconstants,zresidualconstants -P(1)/4, strcat(char('\sigma_x ='), num2str(sigx3D(k)), ' \pm ', num2str(sigxerror3D(k))),'FontSize',textsize);
                text(xresidualconstants,yresidualconstants,zresidualconstants -P(1)/2, strcat(char('\sigma_y ='), num2str(sigy3D(k)), ' \pm ', num2str(sigyerror3D(k))),'FontSize',textsize);
                %surf(ypred, 'Marker', '.'); % Plot fit
                %view([0 90])
                %caxis(cax);
                colorbar
                shading flat
                %axis off;
                grid off;
                subplot(2, 2, 3) %plot residuals
                %form smoothed residuals
                residuals=smoothopticaldepthimagesingle-opticaldepthcalcmatrix;
                temp=reshape(residuals,size(xdata)); %put in matrix form
                smoothresiduals=filter2(b,temp(:,:))/smoothplotsize^2; %weighted residue smoothed a bit for fitting
                surf(pixelconv*xvec,pixelconv*yvec,smoothresiduals, 'Marker', '.');
                zlabel('Optical Depth','FontSize',axisfontsize);
                xlabel('X (microns)','FontSize',axisfontsize);
                ylabel('Y (microns)','FontSize',axisfontsize);
                set(gca,'FontSize', axisfontsize);
                % view([0 90])
                caxis('auto');
                title('true residues (NOT weighted by uncertainties)','FontSize',axisfontsize)
                colorbar
                shading flat
                %axis off;
                grid off;

                subplot(2, 2, 4) % plot weighted residuals
                %form smoothed residuals
                temp=reshape(weightedresidue3D,size(xdata)); %put in matrix form
                smoothweightedresidue3D=filter2(b,temp(:,:))/smoothplotsize^2; %weighted residue smoothed a bit for fitting

                surf(pixelconv*xvec,pixelconv*yvec,smoothweightedresidue3D, 'Marker', '.');
                zlabel('Optical Depth','FontSize',axisfontsize);
                xlabel('X (microns)','FontSize',axisfontsize);
                ylabel('Y (microns)','FontSize',axisfontsize);
                set(gca,'FontSize', axisfontsize);
                % view([0 90])
                caxis('auto');
                title('residues weighted by uncertainties','FontSize',axisfontsize)
                colorbar
                shading flat
                %axis off;
                grid off;
            end % end the if for: skip plot of fit if fit==0
        end %end plot if loop
        
        if plotdensity==1
%             figure(1000+k) %Plot raw data (atoms - back)
%             subplot(2,1,1)
%             plot(Zguess(floor(P(5))-1,:)+Zguess(floor(P(5)),:)+Zguess(floor(P(5))+1,:))
%             hold on
%             %plot(Zguess(floor(P(5))-1,:)+Zguess(floor(P(5)),:)+Zguess(floor(P(5))+1,:))
%             grid off
%             title(basename,'FontSize',axisfontsize)
%             xlabel('x');
%             ylabel('Normalized Density');
%             hold off
%             subplot(2,1,2)
%             plot(Zguess(:,floor(P(4))-1)+Zguess(:,floor(P(4)))+Zguess(:,floor(P(4))+1))
%             hold on
%             %plot(Zguess(:,floor(P(4))-1)+Zguess(:,floor(P(4)))+Zguess(:,floor(P(4))+1))
%             grid off
%             xlabel('y');
%             ylabel('Normalized Density');
%             hold off
            
            figure(2000)
            subplot(4,4,k)
            plot(Zguess(floor(P(5))-1,:)+Zguess(floor(P(5)),:)+Zguess(floor(P(5))+1,:),'-b')
            hold on
            plot(Zguess(:,floor(P(4))-1)+Zguess(:,floor(P(4)))+Zguess(:,floor(P(4))+1),'--r')
            xlabel('x and y position');
            ylabel('Normalized Density');
            if k==8
                legend('x','y','Position','EastOutside');
            end
            title(strcat(num2str(imagevco(k)),'ms'));
            hold off
            
        end %end of plotdensity case structure
                
        %%%%Collect various data and write to outputdatafile, but only if fit is on.
        if fit==1 && column*row<100000; % skip fit if fit==0
            %add data for this file to the alldata matrix
            %detuning = dimage(basenamenumber)*10^6;      % image beam detuning in Hz
            absorptioncross = 6*pi*(lambda/(2*pi))^2/(1+(2*detuning/naturalwidth)^2);     % absorption cross section; used to calculate #atoms
            numatoms = odintegral./absorptioncross;     % number of atoms obtained from the optical depth
            numatomserror = odintegralerror.*numatoms./odintegral;  % dN = N*dOD/OD
            %numatomserror = odintegralerror./absorptioncross;      % dN = dOD/absorptioncross

            %frequencies not important to non-PAS experiments now.
            %dshift = aomvector(basenamenumber)*1000;
            dcat(k)=(-.02288*(vcovoltage(k)^4)+0.497*(vcovoltage(k)^3)-3.436*(vcovoltage(k)^2)+23.93*vcovoltage(k)+169);
            %alldata(k,1)=(dshift(basenamenumber)-dsatabs(basenamenumber)+2*(dzeeman(basenamenumber)+dPAS+(catseyesign(basenamenumber)*dcat(k))))
            alldata(k,1)=(dshift+2*(dzeeman+2*dPAS-dsatabs+(catseyesign*dcat(k))));
            alldata(k,2)=batchvcovoltage(k,basenamenumber);
            alldata(k,3)=odintegralerror;
            alldata(k,4)=odintegral;
            for index=1:(size(output,1));
                alldata(k,index+4)=output(index,k,basenamenumber);
            end; %fit parameters and 1-sigma error for each parameter
            alldata(k,4+size(output,1)+1)=numatomserror;
            alldata(k,4+size(output,1)+2)=numatoms; %Number of atoms is useful
            alldata(k,4+size(output,1)+3)=chi3D(k);
            alldata(k,4+size(output,1)+4)=pixelsize*prebinning; %size parameter saved so hardware and software binning settings only set once in imagefit.
            temp=size(char(file(k,:)));
            allfiles(k)=file(k);

            %write all 2D fit parameters to outputdatafile
            dlmwrite(outputdatafile, alldata,'\t');
            if save2places==1
                dlmwrite(outputdatafile2, alldata,'\t');
            end
            %Columns for alldata
            %1	frequencies (MHz)
            %2	batchvcovoltage
            %3	odintegralerror
            %4	odintegral
            %5	totalod
            %6	od average for selected region
            %7	position x in meters
            %8	position y in meters
            %9	odintegral
            %10	raw 2d output P(1) peak OD
            %11	raw 2d output P(2)sig x in pixels
            %12	raw 2d output P(3)sig y in pixels
            %13	raw 2d output P(4)
            %14	raw 2d output P(5)
            %15 raw 2d output P(6)
            %16 raw 2d output P(7)
            %17 raw 2d output P(8)
            %18 raw 2d output 1-sigma for P(1)
            %19 raw 2d output 1-sigma for P(2)
            %20 raw 2d output 1-sigma for P(3)
            %21 raw 2d output 1-sigma for P(4)
            %22 raw 2d output 1-sigma for P(5)
            %23 raw 2d output 1-sigma for P(6)
            %24 raw 2d output 1-sigma for P(7)
            %25 raw 2d output 1-sigma for P(8)
            %26	number of atoms error
            %27	number of atoms
            %28 chi3D
            %29 size parameter (pixelsize*prebinning)
            %dlmwrite(outputfilesfile, allfiles,'\t')

            if k==q
%                 %Plot number of atoms, peak OD, and sizes each versus
%                 changed variable.
%                 figure(5000+filecounter)
%                 subplot(2,2,1)
%                 hold on
%                 if voltagetodetuning == 1
%                     imagevcod = imagevco * 6.5;
%                     plot(imagevcod(1:counter),alldata(:,27)./10^6,'*b')
%                     xlabel('Detuning in MHz')
%                     ylabel('Number of atoms (10^6)');
%                 end
            else
                figure(1)
                subplot(2,2,1)
                plot(imagevco(1:counter),alldata(:,27)./10^6,'*b')
                hold on
                %xlabel('DAC voltage');
                set(gca, 'FontSize', axisfontsize);
                ylabel('Number of atoms (10^6)');
                title(char(basename));
                subplot(2,2,2)
                hold on
                plot(imagevco(1:counter),alldata(:,10),'*b')
                set(gca, 'FontSize', axisfontsize);
                ylabel('Peak OD');                
                subplot(2,2,3)
                hold on
                plot(imagevco(1:counter),sigx3D(1:counter),'*b')
                set(gca, 'FontSize', axisfontsize);
                %xlabel('DAC voltage');
                ylabel('RMS width X');
                xlabel('DAC voltage');
                %title('frequency calibration');
                subplot(2,2,4)
                hold on
                plot(imagevco(1:counter),sigy3D(1:counter),'*b')
                set(gca, 'FontSize', axisfontsize);
                xlabel('DAC voltage');
                ylabel('RMS width Y');
                hold off
            end
        end % end for skip fit if fit ==0
        %allfiles(k,1:temp(2))=char(file(k,:));
        
    end %loop through files in batch
    %get data from batch file

    avgsigx=sumsigx/k;
    avgsigy=sumsigy/k;
    avgsigz=sqrt(avgsigx*avgsigy);



%     figure(5001)
%     hold on
%     subplot(2,1,1)
%     plot(imagevco(1:counter),sigx3D(1:counter),'*b')
%     %xlabel('DAC voltage');
%     ylabel('RMS width X');
%     %title('frequency calibration');
%     subplot(2,1,2)
%     plot(imagevco(1:counter),sigy3D(1:counter),'*b')
%     xlabel('DAC voltage');
%     ylabel('RMS width Y');
%     hold off
    
    %%%%Spectrum fit
    if spectrumfit==1; % skip fit of spectrum if spectrumfit==0
        %%%%%%%fit spectrum
        %Must produce indevariable, calculatedspectrumcurve, residue1D,
        %fitparameters1D
        if fittingfunction==0
            %%%%%GAUSSIAN FIT
            info='fit is gaussian';
            %guesses
            amplitude=.1;
            center=0;
            FWEM=80;
            InitialConditions1D=[amplitude, center, FWEM];
            lowerbound=[0, -20, 20];
            upperbound=[1, 20, 120];

            alldata1=(alldata(1:q,9)./odsdev(1:q))*(1/sqrt(k));%vector of weighted optical depths
            opticaldepthsfrom2Dfits=alldata(:,9);%vector of optical depths , DATA IS NOT WEIGHTED BY UNCERTAINITIES

            temp=size(alldata(:,1));
            numberofdatapointvector=k*ones(temp(1),1);
            R=[alldata(1:q,1),odsdev(1:q),numberofdatapointvector(1:q)];
            [G,r,J] = nlinfit(R, alldata1(1:q),@gaussian1Dfitfunctionrfch,InitialConditions1D);
            ci2=nlparci(G,r,J);
            linewidth=G(3);
            linewidthsdev=(ci2(3,2)-ci2(3,1))/4;  %turn confidence interval into 1 sigma
            amplitude=G(1);
            amplitudedev=(ci2(1,2)-ci2(1,1))/4; %turn confidence interval for amplitude into 1 sigma
            centerfreq=G(2);
            centerfreqdev=(ci2(2,2)-ci2(2,1))/4;

            %Below are the calculations of residues and chi sqares considering and not considering errors for the 3D fits
            residue1D=(alldata1(:)-gaussian1Dfitfunctionrfch(G, R));%residuals weighted by uncertainties to normalize chi
            gaussian1Dfitfunctionrfch(G, R);
            chi1D=(transpose(residue1D)*residue1D);
            %residue1Dnoerrors=(alldata(:,9)-gaussian1Dfitfunctionrfchnoerro(G, R));%residuals
            %chi1Dnoerrors=(transpose(residue1Dnoerrors)*residue1Dnoerrors)
            fitparameters1D=[1000*delay(k);amplitude;amplitudedev;centerfreq;centerfreqdev;linewidth;linewidthsdev;chi1D;avgsigx;avgsigy;avgsigz];
            alllinewidths(1,basenamenumber)=linewidth;
            alllinewidths(2,basenamenumber)=linewidthsdev;
            alllinewidths(3,basenamenumber)=1000*delay(k);

            %make calculated spectrum, alldata(:,1)=detunings
            points=100;
            for t=1:points
                indevariable(t)=2*min(alldata(:,1))+2*t*(max(alldata(:,1))-min(alldata(:,1)))/points;
                calculatedspectrumcurve(t)=G(1)*exp(-((indevariable(t)-G(2))^2/(G(3))^2)); %Gaussian fitting function
            end
        end   %if fittingfunction==0 ; Gaussian fit loop

        if fittingfunction==1
            %%%%%VOIGT FIT
            info='fit is voigt';
            %guesses
            amplitude=.0001;
            center=0;
            rmsgaussianwidthguess=25;
            scanfreqlimit=10.0*(21.5*2*pi)/(2*pi);
            InitialConditions1D=[amplitude,center,rmsgaussianwidthguess];
            lowerbound=[0, -20, 1];
            upperbound=[1, 20, 40];

            alldata1=(alldata(1:q,9)./odsdev(1:q))*(1/sqrt(k));%vector of weighted optical depths, this was good for gauss1Dfitfunctionrfch
            opticaldepthsfrom2Dfits=alldata(:,9);%vector of optical depths , DATA IS NOT WEIGHTED BY UNCERTAINITIES

            %numberofdatapointvector=(alldata(:,1)./alldata(:,1)*counter);%%%%%%%
            temp=size(alldata(:,1));
            numberofdatapointvector=k*ones(temp(1),1);

            R=[alldata(1:q,1),odsdev(1:q),numberofdatapointvector(1:q)];

            [G,r,J] = nlinfit(R, alldata1,@voigt1Dfitfunctionrfchdata,InitialConditions1D);%fit including weights
            ci2=nlparci(G,r,J);
            linewidth=G(3);
            linewidthsdev=(ci2(3,2)-ci2(3,1))/4;%turn confidence interval into 1 sigma
            amplitudepeak=G(1);
            amplitudedev=(ci2(1,2)-ci2(1,1))/4; %turn confidence interval for amplitude into 1 sigma
            centerfreq=G(2);
            centerfreqdev=(ci2(2,2)-ci2(2,1))/4;
            lambda=422*(10)^(-9);
            temperature=(.088/8.315)*(G(3)*(10)^(6)*lambda)^2;  %temperature=((G(3)*lambda)^2)*molarmass/R
            errortemp1=2*(G(3)-ci2(3,1))*temperature/G(3);
            errortemp2=2*(G(3)-ci2(3,2))*temperature/G(3);
            temperaturebound=[temperature-errortemp1,temperature-errortemp2];
            %Below are the calculations of residues and chi sqares
            residue1D=(alldata1(:)-voigt1Dfitfunctionrfchdata(G, R));%residuals including weight
            %calculatedspectrum=voigt1Dfitfunctionrfchdata(G, R);
            chi1D=(transpose(residue1D)*residue1D) ;

            Nions=2*pi*avgsigx*avgsigy*amplitudepeak*4/(2*pi*21.5)*1e18*3
            density=4*amplitudepeak/(sqrt(2*pi))/(2*pi*21.5)/avgsigz*1e18*3 %m^-3
            %residue1Dnoerrors=(alldata(:,9)-voigt1Dfitfunctionrfchdatanoerro(G, R));%residuals
            %chi1Dnoerrors=(transpose(residue1Dnoerrors)*residue1Dnoerrors)
            fitparameters1D=[1000*delay(k);timezero;gatewidth;density;Nions;amplitudepeak;amplitudedev;centerfreq;centerfreqdev;linewidth;linewidthsdev;chi1D;temperature;errortemp1; errortemp2];
            alllinewidths(1,basenamenumber)=linewidth;
            alllinewidths(2,basenamenumber)=linewidthsdev;
            alllinewidths(3,basenamenumber)=1000*delay(k);
            %CALCULATE VOIGT SPECTRUM BASED ON FIT PARAMETERS
            points=100;
            for t=1:points
                indevariable(t)=2*min(alldata(:,1))+2*t*(max(alldata(:,1))-min(alldata(:,1)))/points;
                Rnoweighting=[indevariable(t),ones(1,1),ones(1,1)];
                calculatedspectrumcurve(t)=voigt1Dfitfunctionrfchdata(G,Rnoweighting);
                %   fitcurve(t)=G(1)*exp(-((indevariable(t)-G(2))^2/(G(3))^2)); %Gaussian fitting function
            end
        end % end  if fittingfunction==1 , do voigt fit

        if fittingfunction==2
            %%%%%LORENTZIAN FIT
            info='fit is lorentzian';
            %guesses
            baselineguess = alldata(1,19);
            amplitudeguess = 0.2;
            centerguess = (max(alldata(:,1))+min(alldata(:,1)))/2;
            linewidthguess = 60;
            bgslopeguess = (alldata(q,19)-alldata(1,19))/(max(alldata(:,1))-min(alldata(:,1)));
            %0.1;
            InitialConditions1D=[baselineguess,amplitudeguess,centerguess,linewidthguess,bgslopeguess];
            lowerbound=[10^5,0,min(alldata(:,1)),1,0];
            upperbound=[10^7,1,max(alldata(:,1)),1000,10];

            alldata1=(alldata(1:q,19)./alldata(1:q,18))*(1/sqrt(k));%vector of weighted optical depths, this was good for gauss1Dfitfunctionrfch
            opticaldepthsfrom2Dfits=alldata(:,9);%vector of optical depths , DATA IS NOT WEIGHTED BY UNCERTAINTIES

            temp=size(alldata(:,1));
            numberofdatapointvector=k*ones(temp(1),1);

            R=[alldata(1:q,1),alldata(1:q,18),numberofdatapointvector(1:q)];

            [G,r,J] = nlinfit(R, alldata1,@lorentzian1Dfitfunctionrfch,InitialConditions1D);%fit including weights
            ci2=nlparci(G,r,J);
            %linewidth=G(3);
            %linewidthsdev=(ci2(3,2)-ci2(3,1))/4;%turn confidence interval into 1 sigma
            bg=G(1);
            bgdev=(ci2(1,2)-ci2(1,1))/4; %turn confidence interval for background into 1 sigma
            amplitude=G(2);
            amplitudedev=(ci2(2,2)-ci2(2,1))/4; %turn confidence interval for amplitude into 1 sigma
            centerfreq=G(3);
            centerfreqdev=(ci2(3,2)-ci2(3,1))/4; %turn confidence interval for center frequency into 1 sigma
            linewidth=G(4);
            linewidthdev=(ci2(4,2)-ci2(4,1))/4; %turn confidence interval for width into 1 sigma
            bgslope=G(5);
            bgslopedev=(ci2(5,2)-ci2(5,1))/4; %turn confidence interval for bgslope into 1 sigma
            %lambda=422*(10)^(-9);

            %Below are the calculations of residues and chi sqares
            residue1D=(alldata1(:)-lorentzian1Dfitfunctionrfch(G, R));%residuals including weight
            %calculatedspectrum=lorentzian1Dfitfunctionrfchdata(G, R);
            chi1D=(transpose(residue1D)*residue1D) ;
            
            Nions=2*pi*avgsigx*avgsigy*amplitude*4/(2*pi*21.5)*1e18*3
            density=4*amplitude/(sqrt(2*pi))/(2*pi*21.5)/avgsigz*1e18*3 %m^-3
            %residue1Dnoerrors=(alldata(:,9)-lorentzian1Dfitfunctionrfchdatanoerro(G, R));%residuals
            %chi1Dnoerrors=(transpose(residue1Dnoerrors)*residue1Dnoerrors)
            fitparameters1D=[1000*delay(k);density;Nions;amplitude;amplitudedev;centerfreq;centerfreqdev;linewidth;linewidthdev;bg;bgdev;bgslope;bgslopedev;];
            alllinewidths(1,basenamenumber)=linewidth;
            alllinewidths(2,basenamenumber)=linewidthdev;
            alllinewidths(3,basenamenumber)=1000*delay(k);
            %CALCULATE LORENTZIAN SPECTRUM BASED ON FIT PARAMETERS
            points=100;
            for t=1:points
                indevariable(t)=min(alldata(:,1))+t*(max(alldata(:,1))-min(alldata(:,1)))/points;
                calculatedspectrumcurve(t) = (bg+bgslope*(indevariable(t)-centerfreq)).*(1-amplitude*((linewidth.^2)/4)./((indevariable(t)-centerfreq).^2+(linewidth.^2)/4));
            end
        end % end  if fittingfunction==2 , do lorentzian fit

        filenumbers=transpose([1:k]);
        datamatrix=[opticaldepthsfrom2Dfits(1:q),alldata(1:q,1),odsdev(1:q),transpose(chi3D(1:q)),residue1D(1:q),filenumbers(1:q),alldata(1:q,3),alldata(1:q,4)];
        % the columns of this matrix are the OD(2D fits), OD error (2D fits), chi
        % square 3D, 1D residue no errors, file numbers, odintegralerror,
        % odintegral
        %dlmwrite(outputfiledatamatrix,datamatrix,'\t',0,0);
        dlmwrite(outputfiledatamatrix,datamatrix,'\t');
        dlmwrite(fitdata1D,transpose(fitparameters1D),'\t');
        if save2places==1
            dlmwrite(outputfiledatamatrix2,datamatrix,'\t');
            dlmwrite(fitdata1D2,transpose(fitparameters1D),'\t');
        end
        
        %%%%%MAKE PLOTS OF SPECTRUM
        spectrumplotnumber=400%+basenamenumber
        spectrumstatsplotnumber=600%+basenamenumber
        figure(spectrumplotnumber)
        %alldata(:,1) is peak optical density from 2D fits for each file in batch
        h=errorbar(alldata(1:q,1),alldata(1:q,19),alldata(1:q,18), 'rx'); % plot freq vs OD and ODstdev
        hold on
        % This sub-section exists to plot different scans with different colors.
        % The division by 40 only works if you have exactly 40 data points; it isn't
        % universal.  I will find a more elegant solution later.
        % If you don't get the # of points correct, the color defaults to red which is
        % an acceptable degradation.
        ccounter=counter/40;
        if ccounter==1
            set(h,{'Color'},{[0,0,1]}) % Blue
        end
        if ccounter==2
            set(h,{'Color'},{[0,1,0]}) % Green
        end
        if ccounter==3
            set(h,{'Color'},{[1,0,0]}) % Red
        end
        if ccounter==4
            set(h,{'Color'},{[1,0,1]}) % Pink
        end
        if ccounter==5
            set(h,{'Color'},{[0,1,1]}) % Aqua
        end
        h=plot(indevariable, calculatedspectrumcurve,'-b');
        text(min(alldata(:,1)),min(opticaldepthsfrom2Dfits),num2str(fitparameters1D));
        text(max(alldata(:,1)),max(opticaldepthsfrom2Dfits),info);
        text(max(alldata(:,1)),max(opticaldepthsfrom2Dfits)/2,num2str(chi1D));
        if scalebaselevel==1
            text(max(alldata(:,1)),max(opticaldepthsfrom2Dfits)*0.75,'with base level scaling')
        end
        title(fitdata1D)
        hold off

        figure(spectrumstatsplotnumber)
        subplot(2,2,1)
        h=plot(alldata(1:q,1),residue1D(1:q), 'rx');
        text(min(alldata(:,1))-(1.5*k),min(opticaldepthsfrom2Dfits)+.1*(max(residue1D)-min(residue1D)),outputdatafile)
        title('1D residuals')

        subplot(2,2,2)
        h=plot(alldata(1:q,1),transpose(chi3D(1:q)), '+b');
        title('3D chi square values')

        subplot(2,2,3)
        h=plot(alldata(1:q,1),sigx3D(1:q),'+b');
        title('Sigma x')

        subplot(2,2,4)
        h=plot(alldata(1:q,1),sigy3D(1:q),'+b');
        title('Sigma y')

    end % end for: skip fit of spectrum if spectrumfit==0
end %end loop through batches

%% Wrap Up
status = 'done'