%FOR 1-BODY LOSSES - NUMBER LOSS EQUATION
function N= Lifetime(coeff,data);

%fit parameter(s)
N0=coeff(1);
beta=coeff(2);

%x axis data
t=data(:,1);

%equation to fit to
%N = N0.*10^6.*exp(-gamma.*t)./(1+(N0.*10^6.*beta./gamma).*(1-exp(-gamma.*t)));
N = N0.*exp(-beta.*t);