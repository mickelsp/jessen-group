%INCLUDES 2-BODY LOSSES - NUMBER/FLUORESCENCE LOSS EQUATION
function Number= Fluoresce2body(coeff,data)

%fit parameter(s)
Gamma=coeff(1);
betap=coeff(2);
Ln=coeff(3);

%data
t=data(:,1);
%v=data(:,2);

%COMPLETE EQUATION
% Nmax=(-Gamma+sqrt(Gamma.^2+4.*beta.*Ln))./(2.*beta);
% Fluo=Nmax.*(1-exp(-(Gamma+(Nmax.*beta./(sqrt(2).*V))).*t))./(1+(beta.*Nmax./(beta.*Nmax+2.*sqrt(2).*V.*Gamma)).*exp(-(Gamma+(Nmax.*beta./(sqrt(2).*V))).*t));

%INCOMPLETE EQUATION
Nss=(-Gamma+sqrt(Gamma.^2+4.*betap.*Ln))./(2.*betap);
CHI=betap.*Nss./(betap.*Nss+Gamma);
Number=Nss.*(1-exp(-(Gamma+2.*betap.*Nss).*t))./(1+CHI.*exp(-(Gamma+2.*betap.*Nss).*t));
%%% If Nss is in millions, then L is millions/second and betap is
%%% m^3*millions/second assuming that volume is in m^3.



