%% Overview of Air Force Test Target Image System Resolution
%%% Fits groups/elements from Air Force Test Target to tell us the
%%% resolution of our imaging system
% http://www.edmundoptics.com/technical-support/technical-library/eo-tech-tools/index.cfm?techToolID=6
%%% Uses 1D Gaussian with 3 peaks for fitting

close all
clear all

%% Import Image File Names
%%% Data from Basler Scout is in .bmp format. Batch file contains list of
%%% all file names.
directory=char('/Users/Work/Documents/Free Space Coupling/Data/2012.02.07/');
imagefilename=char('O17_FORT_Cloud_60MHz_TS_1000us.gif');

%% Define Options
%%%Plot options
fontsize=14; %Set font size for plots
map=colormap(gray); %Set color map for plots
%%%Region of interest options
fullroi=0; %1=use the full image dimensions; 0=use some portion of the image
variableROI=0;  %1 to vary ROI dynamically for falling cloud
%%%Fitting Options
fit=0; %1 to fit the 2D gaussian; 0 to just display images
%%%Camera Parameters
pixelsize=20.41e-6; %[m] size of one camera pixel in meters (depends on magnification); from calibration images Enrique took 2011.12.05

%% Open Image File
imagefile=char(strcat(deblank(directory),deblank(imagefilename)));
imagedata=imread(imagefile);
imagedata=double(imagedata);
ydim=size(imagedata,1);
xdim=size(imagedata,2);

%% Prepare Image Data
%%% Define a region of interest; done before any other array modification
if fullroi==1
    roirowstart=1;     %Starting row for region of interest (min = 1)
    roirowstop=ydim;    %Ending row for region of interest (max = 494)
    roicolstart=1;     %Starting column for region of interest (min = 1)
    roicolstop=xdim;    %Ending column for region of interest (max = 656)
elseif fullroi==0
    if variableROI==0
        %%% Use these numbers if all files should have the same ROI
        roirowstart=13;   %Starting row for region of interest (min = 1)
        roirowstop=439;    %Ending row for region of interest (max = 494)
        roicolstart=80;   %Starting column for region of interest (min = 1)
        roicolstop=505;    %Ending column for region of interest (max = 656)
        
    elseif variableROI==1
        %%%For variable region of interest; to use, set variableROI flag to 1
        rowstartvector=[227 247 257 267 277 287 297 307]; %Define the row number at the top of the region of interest for each file in the batch list
        rowstopvector=[328 338 348 388 428 488 428 448]; %Define the row number at the bottom of the region of interest for each file in the batch list
        colstartvector=[148 148 128 128 128 128 128 128]; %Define the column number at the left of the region of interest for each file in the batch list
        colstopvector=[359 359 379 379 379 379 379 379]; %Define the column number at the right of the region of interest for each file in the batch list
        
        roirowstart=rowstartvector(filenamenumber);
        roirowstop=rowstopvector(filenamenumber);
        roicolstart=colstartvector(filenamenumber);
        roicolstop=colstopvector(filenamenumber);
    end
end

roiimage=imagedata(roirowstart:roirowstop,roicolstart:roicolstop);
roisizes=size(roiimage);
column=roisizes(1);
row=roisizes(2);

imageforanalysis=roiimage;

%% Fit Each Image
if fit==1
    %%%Generate guesses and bounds
    amplitudeguess=100;widthguess=1.5;
    amplitude1=amplitudeguess;amplitude2=amplitudeguess;amplitude3=amplitudeguess;
    center1=4;center2=9;center3=14;
    %InitialConditions=[amplitude1,center1,width1,amplitude2,center2,width2,amplitude3,center3,width3]; %Initial guesses for fitting of atom cloud with 2D gaussian
    InitialConditions=[center1,center2,center3,amplitudeguess];
    xvec=1:row;
    yvec=1:column;
    [xdata,ydata]=meshgrid(xvec,yvec);
    numberofpoints=column*row*ones(column,row); %Matrix in which each entry is the number of points
    barwidth=widthguess*ones(column,row); %Matrix in which each entry is the number of points
    coordmatrix=[xdata(:),ydata(:),numberofpoints(:),barwidth(:)]; %This matrix contains the x and y data points, the normalized error for each pixel, and the total number of points
    imagelist=(imageforanalysis(:).*sqrt(column*row));  %Flattened(1D) list of image weighted by uncertainties, not smoothed
    
    %%% Actual fitting
    [P,r,J]=nlinfit(coordmatrix,imagelist(:),@GaussianMultiPeakFitFunction,InitialConditions); %Fitting - expecting data divided by weight and sqrt(number of points)
    
    center1=P(1);
    center2=P(2);
    center3=P(3);
    width=widthguess;
    amplitude=P(4);
    
    x=0:0.1:50;
    Z=amplitude*1./(1+(2*(x-center1)/width).^2)+amplitude*1./(1+(2*(x-center2)/width).^2)+amplitude*1./(1+(2*(x-center3)/width).^2);
    
end %end of fit case structure

if fit==1
    figure(1)
    h1=image(imageforanalysis);
    set(h1,'CDataMapping','scaled');
    hold on
    %axis off
    plot(x,Z./250,'-b')
    %ylim([0,10000]);
elseif fit==0
    figure(1)
    image(imageforanalysis)
    h1=image(imageforanalysis);
    colormap('bone')
    set(h1,'CDataMapping','scaled');
    hold on
    axis off
end %of fit plotting case structure

%DONE