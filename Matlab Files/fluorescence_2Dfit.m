%% Overview of 2D Fitting Program
%%% Fits Atom Cloud Sizes
%%% Uses 2D gaussian fit; based on imagefit_binaryread from Rice

close all
clear all

%% Import Image File Names
%%% Data from Basler Scout is in .bmp format. Batch file contains list of
%%% all file names.
directory=char('/Users/Work/Documents/Spin Squeezing/Data/2011.12.13/');
batchfilename=char('20111213_FORT_Lifetime_Fluorescence');

batchfile=char(strcat(deblank(directory),deblank(batchfilename),deblank('.batch'))); %Create a text file with a list of filenames to be analyzed

[filenamevector delaytimevector]=textread(batchfile,'%s%f','commentstyle','matlab'); %Read in list of image file names
numberoffilenames=size(filenamevector,1);

%% Define Options
%%%Plot options
fontsize=14; %Set font size for plots
map=colormap(jet); %Set color map for plots
%%%Region of interest options
fullroi=1; %1=use the full image dimensions; 0=use some portion of the image
variableROI=0;  %1 to vary ROI dynamically for falling cloud
%%%Image size reduction options
softwarebinsize=1; %Number of pixels square to sum/bin (i.e. area of bin, in # of pixels, is softwarebinsize*softwarebinsize)
%%%Fitting Options
fit=0; %1 to fit the 2D gaussian; 0 to just display images
numberforscale=5; %Size of corner area used to determine background light level is numberforscale*numberforscale.  Also used to determine error/pixel (noise).
guesssmoothsize=2; %This variable is used to smooth out image data to generate image fit guesses
%%%Camera Parameters
pixelsize=20.41e-6; %[m] size of one camera pixel in meters (depends on magnification); from calibration images Enrique took 2011.12.05

%% Loop Through All Image Files
for filenamenumber=1:numberoffilenames
    clear coordmatrix xdata ydata errorod errorodplot numberofpoints
    filenamenumber
    
    %% Open Image File
    imagefile=char(strcat(deblank(directory),deblank(filenamevector(filenamenumber)),deblank('.bmp')));
    
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
            roirowstart=250;   %Starting row for region of interest (min = 1)
            roirowstop=350;    %Ending row for region of interest (max = 494)
            roicolstart=1;   %Starting column for region of interest (min = 1)
            roicolstop=656;    %Ending column for region of interest (max = 656)
            
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
    roicol=roisizes(1);
    roirow=roisizes(2);
    
    %%% Bin adjacent pixels together
    numcolbins=roicol/softwarebinsize;numrowbins=roirow/softwarebinsize;
    binnedimage=zeros(numcolbins,numrowbins);
    if softwarebinsize > 1
        for i=1:numcolbins
            for j=1:numrowbins
                binnedimage(i,j)=sum(sum(roiimage(i*softwarebinsize-(softwarebinsize-1):i*softwarebinsize,j*softwarebinsize-(softwarebinsize-1):j*softwarebinsize)));
            end
        end
    else
        binnedimage = roiimage;
    end
    
    sizefull=size(binnedimage);    %Size of matrix post-binning
    column=sizefull(1);
    row=sizefull(2);
    
    imageforanalysis=binnedimage;
    
    %% Fit Each Image
    if fit==1
        if column*row>=100000 %Condition insures that Matlab doesn't get an array so large that it will take forever to process
            error('Analysis stopped because array is larger than 100,000 pixels.')
        else
            %%%Generate guesses and bounds
            guessimage=filter2(ones(guesssmoothsize,guesssmoothsize),imageforanalysis(:,:))/guesssmoothsize^2; %Image smoothed lots for finding guesses
            amplitude=max(max(guessimage)); %Guess at initial peak amplitude of fluorescence signal
            sigx=row/10; sigy=column/10; %Guess at the initial x and y sizes of atom cloud
            %To determine guesses for x and y position of cloud, bin by 4 and
            %choose max signal bin; multiply by 4 to get coordinates in full
            %size array
            if softwarebinsize > 5
                centerfindbinsize=1; %Bin by 1 is sufficient in this case because binning is already pretty high
            elseif softwarebinsize <= 5
                centerfindbinsize=4; %Bin by at least 4 because binning is low to begin with
            end
            for i=1:size(guessimage,1)/centerfindbinsize
                for j=1:size(guessimage,2)/centerfindbinsize
                    guessimagebinned(i,j)=sum(sum(guessimage(i*centerfindbinsize-(centerfindbinsize-1):i*centerfindbinsize,j*centerfindbinsize-(centerfindbinsize-1):j*centerfindbinsize)));
                end
            end
            [C,I] = max(guessimagebinned); [C2,I2]=max(max(guessimagebinned));
            xpos=I2*centerfindbinsize; ypos=I(I2)*centerfindbinsize; %Guess at initial x and y coordinates of atom cloud
            
            sumylowxlow=sum(sum(guessimage(1:numberforscale,1:numberforscale))); %sum of pixels in top lefthand corner
            sumylowxhigh=sum(sum(guessimage(1:numberforscale,row-numberforscale:row))); %sum of pixels in top righthand corner
            sumyhighxlow=sum(sum(guessimage(column-numberforscale:column,1:numberforscale))); %sum of pixels in bottom lefthand corner
            sumyhighxhigh=sum(sum(guessimage(column-numberforscale:column,row-numberforscale:row))); %sum of pixels in bottom righthand corner
            offset=(sumylowxlow+sumylowxhigh+sumyhighxlow+sumyhighxhigh)/(4*numberforscale^2); %Mean signal per pixel that lives in the corner regions. Guess at the offset present in the image due to background light
            
            slopex=((sumylowxhigh+sumyhighxhigh)-(sumylowxlow+sumyhighxlow))/(numberforscale.^2.*row); %Guess at initial x slope (if image has an overall x slope in signal intensity)
            slopey=((sumyhighxlow+sumyhighxhigh)-(sumylowxlow+sumylowxhigh))/(numberforscale.^2.*column); %Guessa t initial y slope (if image has an overall y slope in signal intensity)
            
            InitialConditions=[amplitude,sigx,sigy,xpos,ypos,offset,slopex,slopey]; %Initial guesses for fitting of atom cloud with 2D gaussian
            
            %%%Determine error per pixel
            average1=sum(sum(imageforanalysis(1:numberforscale,1:numberforscale)))/numberforscale^2;
            rms1=sqrt(sum(sum((imageforanalysis(1:numberforscale,1:numberforscale)).^2))/numberforscale^2);
            rmsdeviation1=sqrt(rms1^2-average1^2);
            
            average2=sum(sum(imageforanalysis(1:numberforscale,row-numberforscale:row)))/numberforscale^2;
            rms2=sqrt(sum(sum((imageforanalysis(1:numberforscale,row-numberforscale:row)).^2))/numberforscale^2);
            rmsdeviation2=sqrt(rms2^2-average2^2);
            
            average3=sum(sum(imageforanalysis(column-numberforscale:column,1:numberforscale)))/numberforscale^2;
            rms3=sqrt(sum(sum((imageforanalysis(column-numberforscale:column,1:numberforscale)).^2))/numberforscale^2);
            rmsdeviation3=sqrt(rms3^2-average3^2);
            
            average4=sum(sum(imageforanalysis(column-numberforscale:column,row-numberforscale:row)))/numberforscale^2;
            rms4=sqrt(sum(sum((imageforanalysis(column-numberforscale:column,row-numberforscale:row)).^2))/numberforscale^2);
            rmsdeviation4=sqrt(rms4^2-average4^2);
            
            error=abs((rmsdeviation1+rmsdeviation2+rmsdeviation3+rmsdeviation4)/4);
            
            errorod(1:row*column)=error/(sum(sum(imageforanalysis))/(row*column)); %Normalized error - relative to average pixel value
            xvec=1:row;
            yvec=1:column;
            [xdata,ydata]=meshgrid(xvec,yvec);
            numberofpoints=column*row*ones(column,row); %Matrix in which each entry is the number of points
            coordmatrix=[xdata(:),ydata(:),errorod(:),numberofpoints(:)]; %This matrix contains the x and y data points, the normalized error for each pixel, and the total number of points
            imagelist=(imageforanalysis(:)./(errorod(:).*sqrt(column*row)));  %Flattened(1D) list of image weighted by uncertainties, not smoothed
            
            %%% Actual fitting
            [P,r,J]=nlinfit(coordmatrix,imagelist(:),@gaussfitfunctionerrorsslope,InitialConditions); %Fitting - expecting data divided by weight and sqrt(number of points)
            
            ci=nlparci(P,r,J); %error in the fit parameters 95% confidence interval
            saveci(filenamenumber,:,:)=ci;
            sdev=(saveci(:,:,2)-saveci(:,:,1))/4; %turn confidence interval into 1 sigma
            weightedresidue3D=(imagelist(:)-gaussfitfunctionerrorsslope(P,coordmatrix)); %weighted residues (considers errors)
            chi3D(filenamenumber)=(transpose(weightedresidue3D)*weightedresidue3D); %chi squared value for this fit (considers errors)
            
            %%% Remove slope and background from fitted image
            gradmatrix=(xdata(:)-P(4))*P(7)+(ydata(:)-P(5))*P(8)+P(6);
            nogradoddata=imagelist(:)-gradmatrix(:);
            nogradodmatrix=reshape(nogradoddata,size(xdata));
            
            %%% Prepare fit parameters and derived quantities for output
            onesigma = (ci(:,2)-ci(:,1))./4; %1-sigma "error" for fit parameters
            onesigma = transpose(onesigma);
            sigx3D(filenamenumber)=pixelsize*abs(P(2));%x cloud size in meters
            sigy3D(filenamenumber)=pixelsize*abs(P(3));%y cloud size in meters
            sigxerror3D(filenamenumber)=pixelsize*abs(sdev(filenamenumber,2));%error in x cloud size in meters
            sigyerror3D(filenamenumber)=pixelsize*abs(sdev(filenamenumber,3));%error in y cloud size in meters
            peakpositionx=pixelsize*P(4);%positionx in meters
            peakpositiony=pixelsize*P(5);%positiony in meters
            alldata=[P,onesigma,pixelsize]; %alldata=transpose(alldata);
            
            %%% Write to output file
            outputdatafile=char(strcat(deblank(directory),deblank(filenamevector(filenamenumber)),deblank('_2Dfitparams.txt')));
            dlmwrite(outputdatafile,alldata,'\t');
            
            %%%Generate plot of optical density based on fitted parameters, remove distortion of data due to dividing by weight and sqrt(number of points)
            errorodplot(1:row*column)=1; %
            numberofpoints=ones(column,row);%
            coordmatrix=[xdata(:),ydata(:),errorodplot(:),numberofpoints(:)];
            gaussiancalcvector=gaussfitfunctionerrorsslope(P, coordmatrix); %2D Gaussian calculated from fit parameters (unweighted by errors)
            gaussiancalcmatrix = flipud(reshape(gaussiancalcvector, size(xdata)));
            
        end %end of matrix size case structure
    end %end of fit case structure
    
    %% Display Each image
    if fit==1
        figure(filenamenumber)
        subplot(2,1,1) %Plot the data
        image(imageforanalysis)
        hold on
        axis off
        %xlabel('x','FontSize',fontsize,'FontWeight','bold');
        %ylabel('y (gravity)','FontSize',fontsize,'FontWeight','bold');
        title(char(strcat(num2str(delaytimevector(filenamenumber)),' ms Drop Time')),'FontSize',fontsize,'FontWeight','bold');
        
        subplot(2,1,2) %Plot the fit
        surf(xvec,yvec,gaussiancalcmatrix,'Marker', '.');
        hold on
        caxis([0 35])
        zlabel('Amplitude (arb.)','FontSize',fontsize);
        xlabel('X (microns)','FontSize',fontsize);
        ylabel('Y (microns)','FontSize',fontsize);
        set(gca,'FontSize', fontsize);
        view([0 90])
        %colorbar
        shading flat
        grid off
        axis off
    elseif fit==0
        figure(filenamenumber)
        image(imageforanalysis)
        hold on
        axis off
        title(char(strcat(num2str(delaytimevector(filenamenumber)),' ms Drop Time')),'FontSize',fontsize,'FontWeight','bold');
    end %of fit plotting case structure
    
    %% Create an image showing a drop time sequence
    figure(1000)
    subplot(7,7,filenamenumber*7)
    surf(xvec,yvec,gaussiancalcmatrix,'Marker', '.');
    caxis([0 35])
    view([0 90])
    shading flat
    grid off
    subimage(imageforanalysis,map)
    hold on
    axis off
    
end %end of loop through all image files