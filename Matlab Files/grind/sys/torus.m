%TORUS   Polar coordinates plot 
%  Create a time plot using polar coordinates. The first two axis of the phase plane.
%
%   Usage:
%   TORUS - create a torus plot with a period of 365 time steps. The x axis starts at -1.
%   TORUS PERIOD XSTART - create torus plot with a period of PERIOD steps. The x axis 
%   starts at -XSTART.
%
%   See also ru, ax    

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function torus(period,increment)
global g_t g_Y t g_grind;
i_parcheck;
if nargin<1
   period=365;
else
   period=i_checkstr(period);
end;
if nargin<2
   increment=1;
else
   increment=i_checkstr(increment);
end;
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
end;
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
if ~iX.isvar|| ~iY.isvar
   error('GRIND:torus:NoStatevars','Error: there are no state variables on the axes, use ax to set the first 2 axes');
end;
H=i_makefig('torus');
hld=ishold;
plot3(sin(g_t/period*2*pi).*(increment+g_Y(:,iX.no)),cos(g_t/period*2*pi).*(increment+ ...
   g_Y(:,iX.no)),g_Y(:,iY.no));
set(gca,'DrawMode', 'fast');
hold on;
box off;
lims=max(abs([get(gca,'xlim'),get(gca,'ylim')]));
zlim=get(gca,'zlim');
plot3([0;0],[0;0],get(gca,'zlim'),'k')
plot3([0;0],[-lims,lims],[zlim(1),zlim(1)],'k');
plot3([-lims,lims],[0;0],[zlim(1),zlim(1)],'k');
xlabel(['sin(t)*' g_grind.xaxis.var]);
ylabel(['cos(t)*' g_grind.xaxis.var]);
zlabel(g_grind.yaxis.var);
t1=0:0.05:6.5;
plot3(sin(t1)*lims,cos(t1)*lims,zlim(1)*ones(1,length(t1)),'k');
plot3(sin(t1)*increment,cos(t1)*increment,zlim(1)*ones(1,length(t1)),'k');
if ~hld
   hold off;
end;
us1.period=period;
us1.increment=increment;
set(H, 'userdata', us1);
