 %% Gaussian Beam Propagation
%%%For some beam size and wavelength of light, determine the beam's
%%%propagation characteristics.
%%%If there is a lens, predict the changed path.
%%%Future additions: propagation of astigmatic beams (x and y different)

close all
clear all

%% Initialization
fontsize = 16;
lambda = 1064e-9; %[m] wavelength of the beam
z = -1:0.0001:2; %[m] longitudinal position

%% Incoming beam characteristics
%%%Incoming beam characteristics could be derived from measurement or they
%%%could be arbitrarily chosen if one were designing a system from scratch.
w0 = 1.828e-3; %[m] 1/e^2 radius of beam at its maximum intensity
z0 = 0; %[m] position of the minimum beam size in lab coordinate system
zR = pi*w0^2/lambda; %[m] Rayleigh range of the beam

qinitialinverse = -1i.*lambda./(pi.*w0.^2); %[1/m] value of the inverse of the complex beam parameter at the waist of the incoming beam
qinitialinversevector = ABCDmatrix_FreeSpace(qinitialinverse,z-z0); %[1/m] vector of inverse of the complex beam parameter as a function of position for incoming beam

incomingbeamsize = BeamSizeFunction(lambda,qinitialinversevector); %[m] beam size as a function of z position

%% Lens 1 Beam Characteristics
f1 = 150e-3; %[m] focal length of lens 1
z1 = 0.010; %[m] position of lens 1 relative to waist of incoming beam
zlens1 = z0+z1; %[m] position of lens 1 relative to lab coordinate system

qbeforelens1inverse=ABCDmatrix_FreeSpace(qinitialinverse,z1); %[1/m] value of the inverse of the complex beam parameter just before lens 1
qafterlens1inverse = ABCDmatrix_ThinLens(qbeforelens1inverse,f1); %[1/m] value of the inverse of the complex beam parameter just after lens 1
qafterlens1inversevector = ABCDmatrix_FreeSpace(qafterlens1inverse,z-zlens1); %[1/m] vector of the inverse of the complex beam parameter after lens 1

afterlens1beamsize = BeamSizeFunction(lambda,qafterlens1inversevector); %[m] 1/e^2 beam radius as a function of z position

%% Lens 2 Beam Characteristics
f2 = 38.1e-3; %[m] focal length of lens 2
z2 = 0.187; %[m] position of lens 2 relative to lens 1 waist
zlens2 = zlens1 + z2; %[m] position of lens 2 relative to lab coordinate system

qbeforelens2inverse=ABCDmatrix_FreeSpace(qafterlens1inverse,z2); %[1/m] value of the inverse of the complex beam parameter just before lens 2
qafterlens2inverse = ABCDmatrix_ThinLens(qbeforelens2inverse,f2); %[1/m] value of the inverse of the complex beam parameter just after lens 2
qafterlens2inversevector = ABCDmatrix_FreeSpace(qafterlens2inverse,z-zlens2); %[1/m] vector of the inverse of the complex beam parameter after lens 2

afterlens2beamsize = BeamSizeFunction(lambda,qafterlens2inversevector); %[m] 1/e^2 beam radius as a function of z position

afterlens2beamsize(30000)
%% Plot beam propagation
zeropoint=300; %position at z vector where zero lives
figure(1)
hold on
plot(z(1:end).*1e2,incomingbeamsize(1:end).*1e3,'-k','LineWidth',1.5) %plots size of beam after lens 2 over entire z range of interest
plot(z(1:end).*1e2,afterlens1beamsize(1:end).*1e3,'-g','LineWidth',1.5) %plots size of beam after lens 3 over entire z range of interest
plot(z(1:end).*1e2,afterlens2beamsize(1:end).*1e3,'-r','LineWidth',1.5) %plots size of beam after lens 3 over entire z range of interest
set(gca,'FontSize',fontsize,'FontWeight','bold');
xlabel('Longitudinal Position [cm]','FontSize',fontsize,'FontWeight','bold');
ylabel('1/e^2 Beam Radius [mm]','FontSize',fontsize,'FontWeight','bold');
%xlim([0 40]);
ylim([0 4]);
legend('Input Beam','After Lens 1','After Lens 2','Location','NorthEast');