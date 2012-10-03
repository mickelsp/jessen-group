%BOXCARINFLOW   Inflow into a boxcartrain
%   This simple function is used to add biomass to a boxcartrain in a model definition. It returns a vector 
%   with all zeros except for the first element, which is filled with the inflow.
%
%
%   Usage:
%   res=BOXCARINFLOW(VAR, INFLOW) - returns a zero vector of sizeof(VAR) with inflow in the first element.
%
%   See also boxcartrain, model, insim

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function res=boxcarinflow(A,inflow)
res=zeros(size(A));
res(1)=inflow;
