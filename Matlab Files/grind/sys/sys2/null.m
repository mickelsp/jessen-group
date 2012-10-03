%NULL   2D phase plane with nullclines
%   Creates nullclines in the 2D phase plane (it is also allowed to have a
%   parameter on the x or y-axes).
%   If there is no variable on the y-axis, a one-dimensional variant is used
%   (itermap or plotdiff)
%    
%   Usage: 
%   NULL - creates nullclines with default accuracy
%   NULL N - creates nullclines using a grid of N x N points (default=50)
%
%   Remark:
%   NULL is also a MATLAB function to determine the null space of a matrix, for help on this function see:
%   help(fullfile(matlabroot,'toolbox','matlab','matfun','null'))
%   You can still use this function (if the first argument is a matrix).
%   
%   
%   See also itermap, plotdiff, null3, ax, growths

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function Z=null(npoints, opt)
if (nargout>0)||((nargin>0)&&(~ischar(npoints)) && (length(npoints)>1))
    %there is a problem that NULL was already defined in MATLAB,
    %if the first argument is a matrix or if there is output
    %it is assumed that the matlab NULL is meant.
    cur=pwd;
    cd(fullfile(matlabroot,'toolbox','matlab','matfun'));
    if nargin==2
        Z=null(npoints,opt);
    elseif nargin==1
        Z=null(npoints);
    else
        Z=null;
    end;
    cd(cur);
    return;
end;
global g_grind ;
i_parcheck;
findnull = 0;
hasdata = 0; %if data are entered i_vector is not evaluated
if nargin == 0
   npoints = 50;
else
   if ischar(npoints) && strncmpi(npoints, '-f', 2)
      if nargin>1
          npoints=i_checkstr(opt);
      else
          npoints = 50;
      end;
      findnull = 1;
   elseif ischar(npoints) && strncmpi(npoints, '-v', 2)
      %plug in data from a user-defined function (used for R* evaluations)
      Vect=i_checkstr(opt);
      hasdata=1;
      npoints=size(Vect,1);
   end;
   npoints = i_checkstr(npoints);
   if (nargin == 2) && ischar(opt) && strncmpi(opt, '-f', 2)
      findnull = 1;
   end;
end;
if findnull
   g_grind.findnull.ndx = []; %option to find the nullclines of other variables
   for i = 1:g_grind.statevars.dim
      nam = i_statevars_names(i);
      if ~strcmp(g_grind.xaxis.var, nam) && ~strcmp(g_grind.yaxis.var, nam)
         g_grind.findnull.ndx = [g_grind.findnull.ndx(:); i];
      end;
   end;
elseif isfield(g_grind, 'findnull')
   g_grind = rmfield(g_grind, 'findnull');
end;
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
if ~(iX.isvar||iX.ispar||iX.isext) || ~(iY.isvar||iY.ispar||iY.isext)
   if isempty(g_grind) || g_grind.solver.isdiffer
      itermap(1);
   else
      plotdiff;
   end;
   return;
end;
if (isempty(iX.no) || isempty(iY.no))
   ax('?');
   errordlg('Cannot create null-isoclines if there are no state variables on the axes, use "phas 2" instead.');
   error('GRIND:null:NoStatevars','null: Cannot create null-isoclines if there are no state variables on the axes, use "phas 2" instead.');
end
if hasdata
    [X,Y] = meshgrid(linspace(g_grind.xaxis.lim(1),g_grind.xaxis.lim(2),npoints),linspace(g_grind.yaxis.lim(1),g_grind.yaxis.lim(2),npoints));
else
    [X, Y, Vect] = i_vector(npoints, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim, [],0);
end;
[H, new] = i_makefig('phase2');
if new
   set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
   set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
end;
set(H, 'Name', 'Phase plane');
oldhold = ishold;
hold('on');
leg = cell(g_grind.statevars.dim, 1);
%Hs = zeros(g_grind.statevars.dim, 1);
colors = {'m', 'b', 'r', 'g','k'};
lines = {':',':.','--','-'};
maxi = size(colors, 2) * size(lines, 2);
nleg = 0;
for i = 1:size(Vect,3)
   j = mod(i, maxi);
   pen = [char(lines{floor(j / size(colors, 2) + 1)}) char(colors{mod(j, size(colors, 2)) + 1})];
%    if g_grind.version.isoctave %not compatible
%        [c] = contour(X, Y, Vect(:, :, i), [0, 0]);
%    %    set(h, 'linewidth', g_grind.pen.linewidth, 'linestyle',pen(1),'edgecolor',pen(2),'facecolor',pen(2)  );
%    else
       set(gca, 'drawmode','fast');
       [c, h] = contour(X, Y, Vect(:, :, i), [0, 0], pen);
       set(h, 'linewidth', g_grind.pen.linewidth);
%   end;
   if ~isempty(c)
      nleg = nleg + 1;
      if g_grind.solver.isdiffer
         leg{nleg} = i_disptext([i_statevars_names(i) '_{t} = ' i_statevars_names(i) '_{t+1}']);
      else
         leg{nleg} = i_disptext([i_statevars_names(i) '''' '=0']);
      end;
      %Hs(nleg) = h(1);
   end;
end;
if new
   i_plotdefaults(gcf);
   set(gca, 'drawmode','fast');
   if nleg > 0
      leg = leg(1:nleg);
      legend(leg{:});
   end;
   h = xlabel(i_disptext(g_grind.xaxis.var));
   set(h, 'fontsize', g_grind.pen.fontsize);
   h = ylabel(i_disptext(g_grind.yaxis.var));
   set(h, 'fontsize', g_grind.pen.fontsize);
end;
if ~hasdata && (g_grind.statevars.dim > 2)
   title(i_disptext(['Valid for ' i_othervars(i_initvar, iX.no, iY.no)]));
end;
set(gca, 'XLim', g_grind.xaxis.lim);
set(gca, 'YLim', g_grind.yaxis.lim);
if ~oldhold
   hold off;
end;
