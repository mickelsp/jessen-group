%VIEWCELLS   view a matrix variable
%  The values of a matrix variable are shown as a movie of 
%  changing colors. The figure includes tools to move forward of backwards.
%  This command is commonly used for spatial explicit models and 
%  cellular automats.
%
%   Usage:
%   VIEWCELLS - shows all vector/matrix variables.
%   VIEWCELLS X- wait X seconds between frames.
%   VIEWCELLS X CMAP - wait X seconds between frames and use the colormap CMAP 
%   (see the MATLAB command COLORMAP for names of colormaps).
%   VIEWCELLS -out (-o) - change the variables to be plotted.
%   VIEWCELLS -out plotno  var [siz1 siz2] - sets the output in a command line:
%   plotno  = number of plot, var = the state or auxiliary variable, [siz1 siz2] 
%   is the size of the variable.
%   VIEWCELLS -default (-d) - reset the default output.
%
%
%   See also vectplot, setcolormap

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function viewcells(t_pause, colmap, fun, siz)
global t g_grind;
i_parcheck;
tframe = [];
makeavi = [];
if (nargin == 0) || strncmpi(t_pause, '-d', 2) %-defaults
   if ~isfield(g_grind, 'viewcells')
      g_grind.viewcells.currno = 1;
      for i = 1:length(g_grind.statevars.vectnames)
         if g_grind.statevars.dims{i}.dim1 * g_grind.statevars.dims{i}.dim2 > 1
            viewcells('-out', i, g_grind.statevars.vectnames{i}, [g_grind.statevars.dims{i}.dim1, g_grind.statevars.dims{i}.dim2]);
         end;
      end;
   end;
   if nargin ~= 0
      return;
   end;
end;
if (nargin == 1) && strncmpi(t_pause, '-l', 2)
   displayout;
   return;
elseif (nargin > 0) && ischar(t_pause) && strncmpi(t_pause, '-o', 2)
   if nargin == 1
      i_viewcellsdlg('init');
      return;
   end;
   %     if isfield(g_grind, 'viewcells') & ~isempty(g_grind.viewcells.vars)
   %       answer={'1', g_grind.viewcells.vars{1}.name,sprintf('[%d %d]',g_grind.viewcells.vars{1}.dims)};
   %    elseif g_grind.statevars.vector
   %        answer={'1', g_grind.statevars.vectnames{1},sprintf('[%d %d]',g_grind.statevars.dims{1},...
   %           dim1, g_grind.statevars.dims{1}.dim2)};
   %     else
   %       error('viewcells applies only for vector variables');
   %    end;
   %    prompt={'The number of the plot:','The state variable/auxiliary variable:','The size of the variable:'};
   %    dlgTitle = 'viewcells -out';
   %    answer = inputdlg(prompt, dlgTitle, 1, answer);
   %   if isempty(answer)
   %        return;
   %    end;
   %    colmap = answer{1};
   %    fun = answer{2};
   %   siz = i_checkstr(answer{3});
   %   end;
   No = i_checkstr(colmap);
   if isempty(No) || (No < 1)
      No = 1;
   end;
   g_grind.viewcells.currno = No;
   if nargin == 4
      siz = i_checkstr(siz);
   elseif nargin == 3
      siz = evalin('base', sprintf('size(%s);', fun));
   end;
   g_grind.viewcells.vars{No}.name = fun;
   g_grind.viewcells.vars{No}.dims = siz;
   return;
elseif nargin < 1
   t_pause = 0.01;
else
   if ischar(t_pause)
      % option: t=xx
      if ~isempty(strfind(t_pause,'t'))&&~isempty(strfind(t_pause, '='))
         f=strfind(t_pause, '=');
         tframe = i_checkstr(t_pause(f(1) + 1:length(t_pause)));
         t_pause = 0.01;
      end;
      %option avi
   elseif isstruct(t_pause)
      % makeavi.filename
      % makeavi.pathname
      % makeavi.fps frames per second (15)
      % makeavi.from t = start
      % makeavi.to  t= end
      % makeavi.axisonly 1/0?
      makeavi= struct('filename','test.avi','pathname','','fps',15,'from',t,'to',t+g_grind.ndays,'axisonly',1);
      f = fieldnames(t_pause);
      for i = 1:length(f)
%         makeavi = setfield(makeavi, f{i}, getfield(t_pause, f{i}));
         makeavi.(f{i}) = t_pause.(f{i});
      end;
      colmap = 'i_grindcolormap';
      t_pause = 0.01;
   end;
end;
if nargin < 2
   colmap = 'i_grindcolormap';
end;
t_pause = i_checkstr(t_pause);
if ~g_grind.statevars.vector
   errordlg('The current model has no matrix state variables');
   error('GRIND:viewcells:NoMatrixVar','The current model has no matrix state variables');
end;

N0 = i_initvar;
if i_settingschanged(N0, g_grind.ndays)
   %  disp('running');
   i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
end;
i_movie('varcontour', t_pause, tframe, colmap, makeavi);

function displayout
global g_grind;
if ~isfield(g_grind, 'viewcells')
   viewcells('-d');
end;
for no = 1:length(g_grind.viewcells.vars);
   if ~isempty(g_grind.viewcells.vars{no}.name)
      fprintf('viewcells -out %d %s [%d %d]\n', no, g_grind.viewcells.vars{no}.name, ...
         g_grind.viewcells.vars{no}.dims);
   end;
end;

