%% Gaussian Beam Fitting
%%%For list of 1/e^2 sizes and positions, fit to beam propagation equation
%%%to get Rayleigh range, waist position, and waist size of beam
%%%Future additions:  fitting in both x and y transverse dimensions

close all
clear all

%% Initialization
fontsize = 16;
lambda = 690e-9; %[m] wavelength of the laser beam
errorflag = 1; %1 to include error bars as weights in fitting; 0 otherwise

%% Import data
%%%Data should be in 2 column format: longitudinal position (z) and 1/e^2
%%%radius of the beam at that longitudinal position
filepath = char('/Users/work/Documents/Nanofibers/Laser Profiling/Blue-detuned Laser Profiling/Telescope 2 Output Side/');
inputfilename = char('ProbeLaser_NanofiberOutput_Telescope1_20130906_redcoord.txt');
filename = strcat(filepath,inputfilename);
[z beamsize]=textread(filename,'%f%f','commentstyle','matlab'); %read in data file, z in cm and beamsize in mm

%% Data preparation
z = z.*1e-2; %[m] converts z from cm to m
beamsize = beamsize.*1e-3; %[m] converts beamsize from mm to m
lambdavector = zeros(1,length(z));
lambdavector(1) = lambda;

%% Fitting
%%%Fit to beam propagation equation
InitialGuess=[1e-4 3]; %[m m] initial guess for waist value and position of waist
datamatrix = [z(:) lambdavector(:)]; %[m m] data matrix contains both the values to be fitted and the fixed parameter (laser wavelength)
if errorflag==1
    errorvector = beamsize./beamsize; %weight all points equally in absence of better information
    %errorvector=[0.8 0.4 0.4 0.2 0.1];
    [P,r,J]=nlinfitweight(datamatrix,beamsize(:),@fittogaussianbeampropagation,InitialGuess,errorvector(:));
else
    [P,r,J]=nlinfit(datamatrix,beamsize(:),@fittogaussianbeampropagation,InitialGuess);
end

fitwaist = P(1) %[m] waist size; minimum beam size
fitwaistposition = P(2) %[m] longitudinal position of minimum beam size

%% Generate vector that represents the fit to the points
stepsize = (max(z)-min(z))/100;
fitz=-3:.001:3; %fitz = min(z):stepsize:max(z);
zR = pi*fitwaist^2/lambda; %[m] Rayleigh range of the beam
fitbeamsize = fitwaist.*sqrt(1+((fitz-fitwaistposition)./zR).^2);%[m] 1/e^2 beam size as a function of z position

%% Plot the data, the fit, and the fit parameters
figure(1)
if errorflag==1
    errorbar(z(:).*1e2,beamsize(:).*1e3,errorvector(:)./10,'sr','MarkerFaceColor','r','MarkerSize',8)
else
    plot(z(:).*1e2,beamsize(:).*1e3,'sr','MarkerFaceColor','r','MarkerSize',8)
end
hold on
plot(fitz(:).*1e2,fitbeamsize(:).*1e3,'-k','LineWidth',1.2)
%ylim([0 (max(beamsize(:))+0.1.*max(beamsize(:))).*1e3]);
ylim([0 1.9]);
%xlim([50 150]);
set(gca,'FontSize',fontsize,'FontWeight','bold');
xlabel('Longitudinal Position [cm]','FontSize',fontsize,'FontWeight','bold');
ylabel('1/e^2 Beam Radius [mm]','FontSize',fontsize,'FontWeight','bold');
text(min(fitz.*1e2),max(fitbeamsize.*1e3),strcat('w0 = ',num2str(fitwaist.*1e3,3),' mm'),'FontSize',fontsize,'FontWeight','bold');
text(min(fitz.*1e2),max(fitbeamsize.*1e3)-max(fitbeamsize.*1e3)./50,strcat('z0 = ',num2str(fitwaistposition.*1e2,3),' cm'),'FontSize',fontsize,'FontWeight','bold');
%text(min(fitz.*1e2),1.2,strcat('w_0 = ',num2str(fitwaist.*1e3,3),' mm'),'FontSize',fontsize,'FontWeight','bold');
%text(min(fitz.*1e2),1,strcat('z_0 = ',num2str(fitwaistposition.*1e2,3),' cm'),'FontSize',fontsize,'FontWeight','bold');
%title('894 nm input to Nanofiber','FontSize',fontsize,'FontWeight','bold');