%GROWTHS   Growth of each state variable
%   Make filled contour plots with the growth of each state
%   variable in the 2D phase plane.
%
%   Usage:
%   GROWTHS - determine the growth in a grid of 40x40 points.
%   GROWTHS N - determine the growth in a grid of NxN points.
%
%   See also null

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function growths(igrid)
global g_grind;
i_parcheck;
if nargin == 0
   igrid = 40;
else
   igrid = i_checkstr(igrid);
end;
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
if g_grind.statevars.dim>20
   disp('Too many state variables only the first 20 shown');
end;
N0 = i_initvar;
if ~(isempty(iX.no) || isempty(iY.no))
   [X, Y, Vect] = i_vector(igrid, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim);
   %  H = figure(i_figno('growths') + 1);
   for i = 1:min(20,g_grind.statevars.dim)
      H = i_makefig('growths',i);
      i_plotdefaults(H);
      set(H, 'Name', ['Growth of ', char(i_statevars_names(i))]);
      surf(X, Y, Vect(:, :, i))
      hold('on');
      contour3(X, Y, Vect(:, :, i), [0, 0], 'k');
      hold off;
      xlabel(i_disptext(g_grind.xaxis.var));
      ylabel(i_disptext(g_grind.yaxis.var));
      zlabel('growth');
      set(gca, 'XLim', g_grind.xaxis.lim);
      set(gca, 'YLim', g_grind.yaxis.lim);
      %shading flat;
      shading('interp');   % plot a surface
      if getrelease<14 % for some reason matlab 7.5 cannot do this lightning
          light;
          lighting('gouraud'); 
          material('dull');%/shiny/metal
      end;
      colorbar;
      if g_grind.statevars.dim > 2
         title(i_disptext(['Valid for ' i_othervars(N0, iX.no, iY.no)]));
      end;
   end;
else
   errordlg('No state variable on the y axis');
   error('GRIND:growths:NoStateVar','growths: No state variable on the y axis');
end;
