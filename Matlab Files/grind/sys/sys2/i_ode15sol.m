function [tout,yout] = i_ode15sol(odefile,tspan,y0,options,varargin)
%yp0=zeros(size(y0));
yp0 = feval(odefile,tspan(1),y0);
[y0,yp0] = decic(@i_currode,tspan(1),y0,[],yp0,[]);
%[TOUT,YOUT] = ODE15I(ODEFUN,TSPAN,Y0,YP0,OPTIONS)
[tout,yout] = ode15i(@i_currode,tspan,y0,yp0,options);

