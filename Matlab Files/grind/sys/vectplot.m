%VECTPLOT   special plots for vector state variables
%   Special function for vector state variable. Different ways of 
%   visualizing vector state variables (contour or surface plot 
%   with on the y-axis the reference number of the vector and on the x axis time,
%   or a movie. It is not possible to use this function for scalar variables.
%
%
%   Usage:
%   VECTPLOT - makes the current plots.
%   VECTPLOT -out (-o) opens dialog box to select the plot.
%   VECTPLOT -out N0 [XAXIS YAXIS] ATYPE - make a 2D plot/movie (if there 
%   is no time on the axes a movie is made). The axes can be any function of the
%   state variables or time t. ATYPE is type of plot (see below).
%   VECTPLOT -out NO [XAXIS YAXIS ZAXIS] ATYPE - the same but now with 3 axes.
%
%   ATYPE can be (see MATLAB help on the functions):
%   surface 
%   pcolor 
%   bar (or bar3)
%   stem (or stem3)
%   scatter
%   line or plot (or plot3)
%
%   VECTPLOT -default (-d) - reset the default output.
%
%   See also time, model, viewcells

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function vectplot(flag, no,xyzax, atype)
global g_grind t;
%colno = 1;
i_parcheck;
%type = 0;
%log1 = 0;
%xa = '';
%ya = '';
if (nargin == 0) || strncmpi(flag, '-d', 2) %-defaults
   if ~isempty(g_grind.statevars.vectnames)
      if ~isfield(g_grind, 'vectplot')  
         g_grind.vectplot.currno = 1;
         no = 1;
         for i = 1:length(g_grind.statevars.vectnames)
            if g_grind.statevars.dims{i}.dim1 * g_grind.statevars.dims{i}.dim2 > 1
               g_grind.vectplot.vars{no}.xax = 't';
               g_grind.vectplot.vars{no}.yax = g_grind.statevars.vectnames{i};
               g_grind.vectplot.vars{no}.zax = '';
               g_grind.vectplot.vars{no}.type = 'pcolor';
               no = no + 1;
            end;
         end;
      end;
   else
      error('GRIND:vectplot:NoVectors','Vectplot can only be used for vector state variables');
   end;
   if nargin ~= 0
      return;
   end;
end;
if (nargin == 1) && strncmpi(flag, '-l', 2)
   displayout;
   return;
elseif (nargin == 1) && strncmpi(flag, '-o', 2)
   %  prompt = {'Number of plot','Variable/function on x-axis (or t)', ...
   %     'Variable/function on y-axis','Variable/function on z-axis',...
   %     'Graph type (scatter/stem/stem3/bar/surface/pcolor)'};
   %     if isfield(g_grind, 'vectplot')
   %        no = g_grind.vectplot.currno;
   %        answer = {sprintf('%d', no), g_grind.vectplot.vars{no}.xax, g_grind.vectplot.vars{no}.yax, ...
   %          g_grind.vectplot.vars{no}.zax, g_grind.vectplot.vars{no}.type};
   %     elseif ~isempty(g_grind.statevars.vectnames)
   %        answer={'1','t',g_grind.statevars.vectnames{1},'','pcolor'};
   %     else
   %        answer={'1','','','','line'};
   %     end;
   %  answer = inputdlg(prompt, 'vectplot', 1, answer);
   %  if isempty(answer)
   %     error('Cancelled');
   %  end;
   %  no = i_checkstr(answer{1});
   %  g_grind.vectplot.vars{no}.xax = answer{2};
   %  if isempty(g_grind.vectplot.vars{no}.xax)
   %     g_grind.vectplot.vars{no}.xax = sprintf('cellno(%s)', g_grind.statevars.vectnames{1});
   %  end;
   %  g_grind.vectplot.vars{no}.yax = answer{3};
   %  g_grind.vectplot.vars{no}.zax = answer{4};
   %  g_grind.vectplot.vars{no}.type = answer{5};
   %  g_grind.vectplot.currno = no;
   i_vectplotdlg('init');
   return;
elseif (nargin > 1) &&  strncmpi(flag, '-o', 2)
   no = i_checkstr(no);
   g_grind.vectplot.currno = no;
   f=[strfind(xyzax, '[') strfind(xyzax,']')];
   if length(f) == 2
      xyzax = xyzax(f(1) + 1:f(2) - 1);
   end;
   f = strfind(xyzax,',');
   if isempty(f)
      f = strfind(xyzax,' ');
   end;
   if isempty(f)
      g_grind.vectplot.vars{no}.xax = sprintf('cellno(%s)', xyzax);
      g_grind.vectplot.vars{no}.yax = xyzax;
      g_grind.vectplot.vars{no}.zax = '';
   elseif length(f) == 1
      g_grind.vectplot.vars{no}.xax = xyzax(1:f - 1);
      g_grind.vectplot.vars{no}.yax = xyzax(f + 1:end);
      g_grind.vectplot.vars{no}.zax = '';
   else
      g_grind.vectplot.vars{no}.xax = xyzax(1:f(1) - 1);
      g_grind.vectplot.vars{no}.yax = xyzax(f(1) + 1:f(2) - 1);
      g_grind.vectplot.vars{no}.zax = xyzax(f(2) + 1:end);
   end;
   if nargin > 3
      if strncmpi(atype, 'su', 2)
         g_grind.vectplot.vars{no}.type = 'surface';
      elseif strncmpi(atype, 'pc', 2)
         g_grind.vectplot.vars{no}.type = 'pcolor';
      elseif strncmpi(atype, 'ba', 2)
         g_grind.vectplot.vars{no}.type = 'bar';
      elseif strncmpi(atype, 'st', 2)
         g_grind.vectplot.vars{no}.type = 'stem';
      elseif strncmpi(atype, 'sc', 2)
         g_grind.vectplot.vars{no}.type = 'scatter';
      elseif strncmpi(atype, 'l', 1)||strncmpi(atype, 'pl', 2)
         g_grind.vectplot.vars{no}.type = 'plot';
      else
         g_grind.vectplot.vars{no}.type = lower(atype);
      end;
      if ~isempty(strfind(atype, '3'))&&isempty(strfind(g_grind.vectplot.vars{no}.type,'3'))
         g_grind.vectplot.vars{no}.type = [g_grind.vectplot.vars{no}.type '3'];
      end;
   else
      g_grind.vectplot.vars{no}.type = 'pcolor';
   end
   return;
end;
for no = 1:length(g_grind.vectplot.vars)
   xax = i_getno(g_grind.vectplot.vars{no}.xax);
   xax.name = g_grind.vectplot.vars{no}.xax;
   xax.t = strcmp(xax.name, 't');
   yax = i_getno(g_grind.vectplot.vars{no}.yax);
   yax.name = g_grind.vectplot.vars{no}.yax;
   yax.t = strcmp(yax.name, 't');
   zax = i_getno(g_grind.vectplot.vars{no}.zax);
   zax.name = g_grind.vectplot.vars{no}.zax;
   zax.t = strcmp(zax.name, 't');
   type = g_grind.vectplot.vars{no}.type;
   if strcmpi(xax.name, sprintf('cellno(%s)', zax.name))
      xax = zax;
      zax = i_getno('');
      zax.name = '';
   end;
   if strcmpi(yax.name, sprintf('cellno(%s)', zax.name))
      yax = zax;
      zax = i_getno('');
      zax.name = '';
   end;
   
   N0 = i_initvar;
   if i_settingschanged(N0, g_grind.ndays)
      i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
   end;
   X = getvalues(xax);
   Y = getvalues(yax);
   Z = getvalues(zax);
   %i = 0;
   %H = i_figno('vectplot') + i;
   %while ishandle(H + i)
   %   i = i + 1;
   %end;
   H = i_makefig('vectplot', no);
   set(H, 'Name', 'Vector plot');
   if (xax.t || yax.t) && isempty(Z)   %2D or 3D with cellno()
      if (size(Y, 2) > 1) && (size(X ,2)==1)
         X = repmat(X, 1, size(Y, 2));
      end;
      if (size(X, 2) > 1) && (size(Y ,2)==1)
         Y = repmat(Y, 1, size(X, 2));
      end;
      
      if yax.t
         H1 = X;
         X = Y';
         Y = H1';
      end;
      switch type
       case 'stem3'
         [x, y] = meshgrid(1:size(Y, 2), 1:size(Y, 1));
         stem3(x,y, Y);
         if yax.t
            xlabel(sprintf('Cell number %s', xax.name));
            ylabel(yax.name);
            zlabel(xax.name);
         else
            xlabel(sprintf('Cell number %s', yax.name));
            xlabel(xax.name);
            zlabel(yax.name);
         end;
       case 'scatter'
         scatter(X(:), Y(:));
       case 'bar'
         bar(X(:), Y(:));
       case 'bar3'
         bar(X, Y); %works only if Y is a matrix
       case 'stem'
         stem(X(:), Y(:));
       case 'plot'
         plot(X(:), Y(:));
       case 'surface'
         surf(Y');
         shading flat;
         if yax.t
            xlabel('Cell number');
            ylabel(yax.name);
            zlabel(xax.name);
         else
            ylabel('Cell number');
            xlabel(xax.name);
            zlabel(yax.name);
         end;
       case 'pcolor'
         pcolor(Y');
         shading flat;
         if yax.t
            xlabel(sprintf('Cell number %s', xax.name));
            ylabel(yax.name);
         else
            ylabel(sprintf('Cell number %s', yax.name));
            xlabel(xax.name);
         end;
         colorbar;
         i_plotdefaults;
       otherwise
         error('GRIND:vectplot:NotImplemented','Not yet implemented');
      end;
   elseif (xax.t || yax.t || zax.t)  && ~isempty(Z)   %not a movie but a 3D figure
      if (size(Y, 2) > 1) && (size(X ,2)==1)
         X = repmat(X, 1, size(Y, 2));
      end;
      if (size(X, 2) > 1) && (size(Y ,2)==1)
         Y = repmat(Y, 1, size(X, 2));
      end;
      switch type
       case 'plot3'
         plot3(X, Y, Z);
       case 'scatter3'
         plot3(X, Y, Z, '.');
       case 'stem3'
         stem3(X, Y, Z);
       case 'scatter'
         error('GRIND:vectplot:NoScatter3','Cannot use scatter in 3D, use scatter3 instead');
       case 'stem'
         error('GRIND:vectplot:NoStem3','Cannot use stem in 3D, use stem3 instead');
       case 'plot'
         error('GRIND:vectplot:NoPlot3','Cannot use plot in 3D, use plot3 instead');
       otherwise
         error('GRIND:vectplot:NotImplemented','Not implemented');
      end;
      xlabel(xax.name);
      ylabel(yax.name);
      zlabel(zax.name);
   elseif ~xax.t && ~yax.t && ~zax.t
      opt.xlabel = xax.name;
      opt.ylabel = yax.name;
      if ~isempty(zax.name)
         opt.zlabel = zax.name;
      end;
      movplot(type, X, Y, Z, opt);
   else
      error('GRIND:vectplot:NotImplemented','Not implemented');
   end;
end;
disp('use <a href="matlab:vectplot -out">vectplot -out</a> to change the output/add plots');
function X = getvalues(ax)
if isempty(ax.name)
   X = [];
   return
else
   X = i_getoutfun(ax.name);
end;
function displayout
global g_grind;
for no = 1:length(g_grind.vectplot.vars);
   if isempty(g_grind.vectplot.vars{no}.zax)
      fprintf('vectplot -out %d [%s %s] %s\n', no, g_grind.vectplot.vars{no}.xax, ...
         g_grind.vectplot.vars{no}.yax, g_grind.vectplot.vars{no}.type);
   else
      fprintf('vectplot -out %d [%s %s %s] %s\n', no, g_grind.vectplot.vars{no}.xax, ...
         g_grind.vectplot.vars{no}.yax, g_grind.vectplot.vars{no}.zax, g_grind.vectplot.vars{no}.type);
   end;
end;
