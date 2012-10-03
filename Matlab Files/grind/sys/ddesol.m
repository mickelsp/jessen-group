function [tout,yout,varargout] = ddesol(odefile,tspan,y0,options,varargin)
global g_grind;
%DDE23('F',LAGS,HISTORY,TSPAN)
lags=zeros(size(g_grind.lags));
varargout=varargin;
for i=1:length(g_grind.lags)
   lags(i)=abs(evalin('base',g_grind.lags{i}));
end;
sol=dde23(odefile,lags,y0,[tspan(1),tspan(length(tspan))],options);
if length(tspan)>2
   [yout]=deval(sol,tspan');
   tout=tspan';
   yout=yout';
else
   yout=sol.y';
   tout=sol.x';
end;

