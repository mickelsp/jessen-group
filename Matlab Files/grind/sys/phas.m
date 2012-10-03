%PHAS   Phase space 
%   Create or open a phase space plot, showing the last run. If
%   there is no run, or if parameters/state variables have been
%   changed, the model is run first.
%
%   Usage:
%   PHAS - if there is a variable for the z-axis selected, create
%   a 3D phase space,
%   else if there is a variable for the y-axis, create a 2D phase plane 
%   else create a 1D phase space. 
%   PHAS 2 - create a 2D phase plane if there is something on the
%   y-axis, else create
%   a 1D phase space.
%   PHAS 3 - create a 3D phase space
%
%   See also ax, null, null3, ru

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function phas(nDim)
global g_grind;
i_parcheck;
if nargin == 0
   if isempty(g_grind.zaxis.var)
      nDim = 2;
   else
      nDim = 3;
   end;
else
   nDim = i_checkstr(nDim);
end;
if nDim == 3
   if isempty(g_grind.zaxis.var)
      ax('?');
      errordlg('Cannot create 3D plot if there is nothing on the Z axis (see AX)');
      error('GRIND:phas:NoZAxis','Cannot create 3D plot if there is nothing on the Z axis (see AX)');
   end;
   [h,fignew] = i_makefig('phase3');
   if fignew
        set(h, 'WindowButtonDown', 'i_callb(''mdown'')');
        set(h, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
     set(gca, 'View', [322.5, 30]);
   end;
elseif isempty(g_grind.yaxis.var)
   h = i_figno('phase1');
else
   h = i_makefig('phase2');
end;
if i_settingschanged
   ru;
else
   i_phas(h, 0);
end;

