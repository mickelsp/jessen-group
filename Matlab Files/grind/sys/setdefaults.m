%SETDEFAULTS   Defaults settings
%   Save/load default settings
%   
%   Usage:
%   SETDEFAULTS saves the default settings
%   SETDEFAULTS 1 reads the default settings
%   SETDEFAULTS -1 resets the default settings to the system defaults

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function def = setdefaults(doload)
global g_grind;
if nargin < 1
   doload = 0;
end;
doload = i_checkstr(doload);
filenam = fullfile(grindpath , 'defaults.mat');
if doload < 0
   delete(filenam);
   setsystemdefaults;
   disp('default settings reset');
elseif doload
   g_grind.odefile=get_odefile('curr_ode');
   g_grind.xaxis.lim = [0, 10];
   g_grind.yaxis.lim = [0, 10];
   g_grind.zaxis.lim = [0, 10];
   id = fopen(filenam, 'r');
   if id > 0
      fclose(id);
      load(filenam);
      if exist('defaults','var')
         g_grind.drawnow = defaults.drawnow; %#ok
         g_grind.slowdown = defaults.slowdown;
         g_grind.tstep = defaults.g_tstep;
         g_grind.ndays = defaults.g_ndays;
         g_grind.pen.nextpen = defaults.g_pen.nextpen;
         g_grind.pen.colormap = defaults.g_pen.colormap;
         g_grind.pen.linewidth = defaults.g_pen.linewidth;
         g_grind.pen.fontsize = defaults.g_pen.fontsize;
         g_grind.pen.tickdir = defaults.g_pen.tickdir;
      end;
   else
      setsystemdefaults;
   end;
else
   defaults.version = 1;
   defaults.drawnow = g_grind.drawnow;
   defaults.slowdown = g_grind.slowdown;
   defaults.g_tstep = g_grind.tstep;
   defaults.g_ndays = g_grind.ndays;
   defaults.g_pen.nextpen = g_grind.pen.nextpen;
   defaults.g_pen.colormap = g_grind.pen.colormap;
   defaults.g_pen.linewidth = g_grind.pen.linewidth;
   defaults.g_pen.fontsize = g_grind.pen.fontsize;
   defaults.g_pen.tickdir = g_grind.pen.tickdir;
   disp(defaults)
   save(filenam, 'defaults');
end;
if nargout == 1
   def = defaults;
end

function setsystemdefaults
global g_grind;
g_grind.slowdown = 0;
g_grind.drawnow = 1;
g_grind.tstep = NaN;
g_grind.ndays = 1000;
g_grind.pen = nextpen([]);
return

function odef = get_odefile(odefile)
j = length(odefile);
while (j > 0) && ~isempty(strfind('1234567890', odefile(j)))
   j = j - 1;
end;
odefile1 = [odefile(1:j) '%d.m'];
i = 0;
while exist(fullfile(grindpath , sprintf(odefile1, i)), 'file') && (i < 200)
   i = i + 1;
end;
odef = sprintf([odefile(1:j) '%d'], i);

