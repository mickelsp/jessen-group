function [beta,r,J] = nlinfitweight(X,y,model,beta0,weight)
global model1 weight1  % these are defined as global variables to be used in the fuction newmodel
model1=model; weight1=weight;
% This function calls the nlinfit.m, but you can also provide uncertainties to 
% your Y values. It modifies Y and Yfit so that the weighted chi^2 is
% minimized. 
% Call the nlinfitweight as nlinfitweight(X,y,model,beta0,weight)
% Given that nlinfit can be called as nlinfit(X,y,model,beta0)
% Typically, X is a design matrix of predictor
% (independent variable) values, with one row for each value in Y.
% However, X may be any array that FUN is prepared to accept.  FUN is
% a function that accepts two arguments, a coefficient vector and the
% array X, and returns a vector of fitted Y values.  BETA0 is a vector
% containing initial values for the coefficients. (For more details
% please see the description of nlinfit. Apart from the above, weight is 
% a vector of the same length as Y, and it is defined as below:-
% chi^2=Sum((Yfiti-Yi)/weighti)^2


ynew=y./weight;
[beta,r,J]=nlinfit(X,ynew,@newmodel,beta0);

%-- the following defines the new function to be used in nlinfit
function yfitnew=newmodel(beta0, X)
global model1 weight1
yfit=feval(model1,beta0,X);
yfitnew=yfit(:)./weight1(:);