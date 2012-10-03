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
z = -2:0.01:2.0; %[m] longitudinal position

%% Import empirically measured data
%%%This section isn't necessary for the general program.  It's here to
%%%double check against the analytically calculated values
filename = char('/Users/work/Documents/Spin Squeezing/Equipment/Amoco Laser PAA0146/BeamProfile_AmocoPAA0146_20110715.txt'); %%Macintosh path format; slashes in opposite direction!
[zexp beamsizeexp]=textread(filename, '%f%f','commentstyle','matlab'); %read in data file, z in cm and beamsize in mm
zoffset = 0; %[cm] distance from laser output to second lens
zexp=zexp+zoffset; %[cm] distance of measurement from output of fiber laser

fitwaist = 1.07e-4; %[m] waist determined from the fit to the experimental data
fitwaistposition = -0.0507; %[m] position of the waist relative to the output of the laser
zRfit = pi*fitwaist^2/lambda; %[m] Rayleigh range of the beam based on the fit parameters
fitbeamsize = fitwaist.*sqrt(1+((z-fitwaistposition)./zRfit).^2);%[m] 1/e^2 beam size as a function of z position

%% Incoming beam characteristics
%%%Incoming beam characteristics could be derived from measurement or they
%%%could be arbitrarily chosen if one were designing a system from scratch.
w0 = 1.07e-4; %[m] 1/e^2 radius of beam at its maximum intensity
z0 = -0.0507; %[m] position of the minimum beam size in lab coordinate system
zR = pi*w0^2/lambda; %[m] Rayleigh range of the beam

incomingbeamsize = w0.*sqrt(1+((z-z0)./zR).^2);%[m] 1/e^2 beam size as a function of z position

%% Define lens 1
f1 = 350e-3; %[m] focal length of lens 1
z1 = 0.400104; %[m] position of lens 1 relative to waist of incoming beam
zlens1 = z0+z1; %[m] position of lens 1 relative to lab coordinate system

incomingR = z1.*(1+zR^2./z1.^2); %[m] radius of curvature of the incoming beam
incomingW = w0.*sqrt(1+(z1./zR).^2);%[m] 1/e^2 beam size of incoming beam at the position of the lens
qbeforelens1inverse = 1./incomingR-1i.*lambda./(pi.*incomingW.^2); %[1/m] value of the inverse of the complex beam parameter just before lens 1

%% Calculate beam characteristics after lens 1
lens1R = 1./((1/incomingR)-(1/f1)); %[m] radius of curvature of beam after lens 1
lens1W = incomingW; %[m] 1/e^2 beam size of beam just after lens 1
qafterlens1inverse = qbeforelens1inverse-1./f1; %[1/m] value of the inverse of the complex beam parameter just after lens 1

afterlens1z0 = f1./(1+(f1/zR).^2); %[m] position of beam waist relative to lens 1
afterlens1z0 = 202.3497e-2;
qwaist1 = 1/qafterlens1inverse+afterlens1z0; %[m] value of the complex beam parameter at the waist after lens 1
w1 = sqrt(lambda.*imag(qwaist1)./pi); %[m] minimum 1/e^2 size of the beam after lens 1
afterlens1zR = (pi.*w1.^2)./lambda; %[m] Rayleigh range of beam after lens 1

afterlens1beamsize = w1.*sqrt(1+((z-(afterlens1z0+zlens1))./afterlens1zR).^2); %[m] 1/e^2 beam size as a function of z position

%% Define lens 2
f2 = 100e-3; %[m] focal length of lens 2
z2 = 0.11; %[m] position of lens 2 relative to lens 1 waist
zlens2 = zlens1 + afterlens1z0 + z2; %[m] position of lens 2 relative to lab coordinate system

incomingR1 = z2.*(1+afterlens1zR^2./z2.^2); %[m] radius of curvature of the beam just before lens 2 (after going through lens 1)
incomingW1 = w1.*sqrt(1+(z2./afterlens1zR).^2); %[m] 1/e^2 beam size of the beam at the position of lens 2 (after going through lens 1)
qlens2 = qwaist1 + z2; %[m] complex beam parameter before lens 2
qlens2inverse = 1./qlens2; %[1/m] value of the inverse of the complex beam parameter before lens 2

%% Calculate beam characteristics after lens 2
lens2R = 1./((1/incomingR1)-(1/f2)); %[m] radius of curvature of beam after lens 2
lens2W = incomingW1; %[m] 1/e^2 beam size of beam just after lens 2
qafterlens2inverse = qlens2inverse-1./f2; %[1/m] value of the inverse of the complex beam parameter just after lens 2

afterlens2z0 = f2./(1+(f2./afterlens1zR).^2); %[m] relative position of beam waist after lens 2
qwaist2 = 1/qafterlens2inverse+afterlens2z0; %[m] complex beam parameter at the position of the waist after lens 2
w2 = sqrt(lambda.*imag(qwaist2)./pi); %[m] minimum 1/e^2 size of the beam after lens 2
afterlens2z0 = abs(real(qwaist2)); %[m] relative position of beam waist after lens 2
afterlens2zR = pi.*w2.^2./lambda; %[m] Rayleigh range of beam after lens 2

afterlens2beamsize = w2.*sqrt(1+((z-(zlens2+afterlens2z0))./afterlens2zR).^2); %[m] 1/e^2 beam size as a function of z position
%afterlens2beamsize(end) %[m] 1/e^2 beam size at the maximum propagation length (z value) defined by this program

% Define lens 3
%%Definition of lenses 3 and 4 is not essential for any calculations done
%%in this program thus far.  Nor have I trouble-shot the code, so I'm not
%%sure it's right.  The principles are the same as in lens 1 and 2 above.
f3 = 1000e-3; %[m] focal length of lens 3
z3 = 0.49; %[m] position of lens 3 relative to lens 2 waist
zlens3 = zlens2 + afterlens2z0 + z3; %[m] position of lens 3 relative to lab coordinate system

incomingR2 = z3.*(1+afterlens2zR^2./z3.^2); %[m] radius of curvature of the beam just before lens 3 (after going through lens 2)
incomingW2 = w2.*sqrt(1+(z3./afterlens2zR).^2); %[m] 1/e^2 beam size of the beam at the position of lens 3 (after going through lens 2)
qlens3 = qwaist2 + z3; %[m] complex beam parameter before lens 3
qlens3inverse = 1./qlens3; %[1/m] value of the inverse of the complex beam parameter before lens 3

%numbers below obtained from ABCD matrix program; my third lens
%calculations aren't working for some reason.
afterlens3z0 = 47.4e-2;
w3 = 524e-6;
afterlens3zR = 81.1e-2;
afterlens3beamsize = w3.*sqrt(1+((z-(zlens3+afterlens3z0))./afterlens3zR).^2); %[m] 1/e^2 beam size as a function of z position

%% Calculate beam characteristics after lens 3
% lens3R = 1./((1./incomingR2)-(1./f3)); %[m] radius of curvature of beam after lens 3
% lens3W = incomingW2; %[m] 1/e^2 beam size of beam just after lens 3
% qafterlens3inverse = qlens3inverse-1./f3; %[1/m] value of the inverse of the complex beam parameter just after lens 3
% 
% %afterlens3z0 = f3./(1+(f3./afterlens2zR).^2); %[m] relative position of beam waist after lens 3
% %afterlens3z0=-afterlens3z0;
% afterlens3z0 = -0.62;
% qwaist3 = 1/qafterlens3inverse+afterlens3z0; %[m] complex beam parameter at the position of the waist after lens 3
% w3 = sqrt(lambda.*imag(qwaist3)./pi); %[m] minimum 1/e^2 size of the beam after lens 3
% afterlens3z0 = abs(real(qwaist3)); %[m] relative position of beam waist after lens 3
% %afterlens3z0=-afterlens3z0;
% afterlens3zR = pi.*w3.^2./lambda; %[m] Rayleigh range of beam after lens 3

% afterlens3z0 = -0.155;
% w3 = 520.0e-6;
% afterlens3zR = 0.797;
% afterlens3beamsize = w3.*sqrt(1+((z-(zlens3+afterlens3z0))./afterlens3zR).^2); %[m] 1/e^2 beam size as a function of z position
%min(afterlens3beamsize);

% %% Define lens 4
% f4 = 300e-3; %[m] focal length of lens 4
% z4 = 0.305; %[m] position of lens 4 relative to lens 3 waist
% zlens4 = zlens3 + afterlens3z0 + z4; %[m] position of lens 4 relative to lab coordinate system
% 
% incomingR3 = z4.*(1+afterlens3zR^2./z4.^2); %[m] radius of curvature of the beam just before lens 4 (after going through lens 3)
% incomingW3 = w3.*sqrt(1+(z4./afterlens3zR).^2); %[m] 1/e^2 beam size of the beam at the position of lens 4 (after going through lens 3)
% qlens4 = qwaist3 + z4; %[m] complex beam parameter before lens 4
% qlens4inverse = 1./qlens4; %[1/m] value of the inverse of the complex beam parameter before lens 4
% 
% %% Calculate beam characteristics after lens 3
% lens4R = 1./((1/incomingR3)-(1/f4)); %[m] radius of curvature of beam after lens 4
% lens4W = incomingW3; %[m] 1/e^2 beam size of beam just after lens 4
% qafterlens4inverse = qlens4inverse-1./f4; %[1/m] value of the inverse of the complex beam parameter just after lens 4
% 
% afterlens4z0 = f3./(1+(f3./afterlens3zR).^2); %[m] relative position of beam waist after lens 4
% qwaist4 = 1/qafterlens4inverse+afterlens4z0; %[m] complex beam parameter at the position of the waist after lens 4
% w4 = sqrt(lambda.*imag(qwaist4)./pi); %[m] minimum 1/e^2 size of the beam after lens 4
% afterlens4z0 = abs(real(qwaist4)); %[m] relative position of beam waist after lens 4
% afterlens4zR = pi.*w4.^2./lambda; %[m] Rayleigh range of beam after lens 4
% 
% afterlens4beamsize = w4.*sqrt(1+((z-(zlens4+afterlens4z0))./afterlens4zR).^2); %[m] 1/e^2 beam size as a function of z position

%% Plot beam propagation
figure(1)
plot(z(:).*1e2,incomingbeamsize(:).*1e3,'-k','LineWidth',1.5) %plots size of incoming beam over entire z range of interest
hold on
plot(z(:).*1e2,afterlens1beamsize(:).*1e3,'--r','LineWidth',1.5) %plots size of beam after lens 1 over entire z range of interest
%plot(z(:).*1e2,afterlens2beamsize(:).*1e3,':b','LineWidth',1.5) %plots size of beam after lens 2 over entire z range of interest
plot(zexp,beamsizeexp,'sg','MarkerFaceColor','g'); %plot experimentally measured beam sizes to compare with analytical beam propagation
%plot(z(:).*1e2,fitbeamsize(:).*1e3,':c','LineWidth',1.5) %plot size of beam derived from experimentally measured beam sizes to compare with analytical beam propagation
set(gca,'FontSize',fontsize,'FontWeight','bold');
xlabel('Longitudinal Position [cm]','FontSize',fontsize,'FontWeight','bold');
ylabel('1/e^2 Beam Radius [mm]','FontSize',fontsize,'FontWeight','bold');
xlim([100*min(z) 100*max(z)]);