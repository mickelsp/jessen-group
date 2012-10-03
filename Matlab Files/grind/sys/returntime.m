%RETURNTIME   Estimates the time necessary to reach a stable node
%   From the current initial conditions, the time is estimated till the change in the 
%   model is negligible. This command analyses the results of the last run. (if there 
%   is no last run or if parameters have changed, it calls RU). 
%   
%   Usage:
%   RETURNTIME - estimates the returntime based on the current SIMTIME settings and a default
%   value for the maximum change in equilibrium (1E-8). 
%   RET=RETURNTIME - store the returntime in the variable RET.
%   RETURNTIME ERR - use the maximum change of ERR.
%   RETURNTIME ERR MAXT - use ERR and simulate MAXT timesteps to find an equilibrium.
%   
%   
%   See also returntime2d, ru, simtime, findeq

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function rettime = returntime(err1, maxt, methodopt)
global g_Y t g_t g_grind;
i_parcheck;
err = 1E-8;
if nargin == 0
   methodopt='-both';
else
   if ischar(err1)&&(err1(1)=='-')
      methodopt=err1;
   else
      err =  i_checkstr(err1);
      if nargin==1
         methodopt='-both';
      end;
   end;
end;
if (nargin  == 2)
   if ischar(maxt)&&(maxt(1)=='-')
      methodopt=maxt;
    elseif ~isnan(maxt)
      g_grind.ndays =  i_checkstr(maxt);
      methodopt='-both';
    else
      methodopt='-both';
   end;
end;
switch lower(methodopt)
    case '-change'
       method=0;
    case '-equil'
       method=1;
    otherwise
       method=2;
end 
maxt=g_grind.ndays;
if method==2
   err2=err*2;
elseif method==1
   err2=err;
   err=1E40; %very larger number: will never been reached
elseif method==0
   err2=1E40;
end;
N0 = i_initvar;
if ~isnan(maxt)||i_settingschanged(N0)
   i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
end;
if length(g_t)-g_grind.tstep<2
   tt=g_t;
   YY=g_Y;
else
   tt=g_t(1):1:g_t(size(g_t,1));
   YY=interp1(g_t,g_Y,tt);
   if size(g_Y,2)==1
      YY=YY';
   end;   
end;
iret = 1;
conv=0;
for i = 1:size(YY, 1) - 1
   diff = mean(abs(YY(i, :) - YY(i + 1, :)));
   diff2 = mean(abs(YY(i, :) - YY(size(YY,1), :)));
   conv= (diff<err) & (diff2<err2);
   if conv
 %     disp(sprintf('iret %0.5g diff %0.5g diff2 %0.5g',[iret,diff,diff2]));
      break;
   end;
   iret = i;
end;
if ~conv
   rettime=tt(length(tt))-tt(1);
   if (nargout == 0)
      disp('no equilibrium reached (cannot detect cycles)');
   end;
else
   rettime =tt(iret) - tt(1);
end;

