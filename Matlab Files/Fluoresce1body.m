%INCLUDES 2-BODY LOSSES - NUMBER/FLUORESCENCE LOSS EQUATION
function N= Fluoresce1body(coeff,data);

%fit parameter(s)
gamma=coeff(1);
Ln=coeff(2);

%data
t=data(:,1);



N=(Ln/gamma).*(1-exp(-gamma.*t));




