%NULL3   3D phase space with nullclines planes
%   Create nullclines in the 3D state space (it is also allowed to have a
%   parameter on one of the axes).
%   Note that there should be state variables or parameters on all 
%   three axes 
%
%   Usage:
%   NULL3  - creates nullcline planes using default settings
%   NULL3 X1 - creates nullclines based on a 3D grid of X1*X1*X1 points (default=25). 
%   
%
%   See also ax, null, phas

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function null3(npoints, opt)
global g_grind;
i_parcheck;
hasdata = 0;
if nargin == 0
   npoints = 25;
   nlines = npoints;
elseif ischar(npoints) && strncmpi(npoints, '-v', 2)
   %plug in data from a user-defined function (used for R* evaluations)
   Vs = i_checkstr(opt);
   hasdata = 1;
   npoints = size(Vs{1}, 1);
   nlines = size(Vs{1}, 3);
else
   npoints = i_checkstr(npoints);
   nlines = npoints;
end;

iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
iZ = i_getno(g_grind.zaxis.var);
if ~(iX.isvar || iX.ispar) || ~(iY.isvar||iY.ispar) ||~(iZ.isvar || iZ.ispar)
   ax('?');
   errordlg('Cannot create 3D null-isoclines if there are no state variables on the axes.')
   error('GRIND:null3:NoStatevariables','Cannot create 3D null-isoclines if there are no state variables on the axes');
end;
N0 = i_initvar;
Zaxis = g_grind.zaxis.lim;
if abs(Zaxis(1)) < 0.0001
   Zaxis(1) = 0.0001;
end
oldZ = i_getparvar(iZ);
Xs = zeros(npoints + 1, npoints + 1, nlines);
Ys = Xs;
Zs = Xs;
%Vs = Xs;
try
   if ~hasdata
      Vs = cell(1,g_grind.statevars.dim);
      for i = 1:g_grind.statevars.dim
         Vs{i} = Xs;
      end;
      d = cell(1, g_grind.statevars.dim);
      for i = 1:size(d, 2)
         d{i} = zeros(3, 1000);
      end;
      %    nd = ones(1, g_grind.statevars.dim);
      %    maxn = repmat(1000,1, g_grind.statevars.dim);
      %    nanrow = repmat(NaN,3, 1);
      wb = waitbar(0, 'Calculating...');
      set(wb,'name','null3')
      for z1 = 1:nlines
         waitbar(z1 / nlines, wb);
         pz = (z1 - 1) * (Zaxis(2) - Zaxis(1)) / nlines + Zaxis(1);
         N0 = i_initvar;
         N0 =  i_setparvar(iZ, N0, pz);
         [X, Y, Vect] = i_vector(npoints, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim, N0);
         Xs(:, :, z1) = X;
         Ys(:, :, z1) = Y;
         for i = 1:size(Vect, 3)
            V = Vs{i};
            V(:, :, z1) = Vect(:, :, i);
            Vs{i} = V;
         end;
         Zs(:, :, z1) = pz * ones(npoints + 1);
      end;
      close(wb);
   else
       [Xs,Ys,Zs]=meshgrid(linspace(g_grind.xaxis.lim(1),g_grind.xaxis.lim(2),npoints),...
           linspace(g_grind.yaxis.lim(1),g_grind.yaxis.lim(2),npoints),...
           linspace(g_grind.zaxis.lim(1),g_grind.zaxis.lim(2),nlines));
   end;
   [H, new] = i_makefig('phase3');
   if new
      set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
      set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
   end;
   set(H,'name', '3D phase space');
   oldhold = ishold;
   hold('on');
   leg = cell(length(Vs), 1);
   col = [0 0 1; 1 0 0;0.5 1 0; 0 1 0; 0 0 0.5; 0.5 0 0; 0.5 1 0; 0 0.5 0];
   for i = 1:length(Vs)
      h = patch(isosurface(Xs, Ys, Zs, Vs{i}, 0));
      camlight('left'); lighting('phong');
      alpha(0.85); %opaquenesss
      material('metal');
      set(h,'FaceColor', col(mod(i, size(col, 1)), :), 'EdgeColor','none');
      leg{i} = [i_statevars_names(i) '''' '=0'];
   end;
   if new
      legend(leg{:});
   end;
   set(gca, 'View', [322.5, 30]);
   box('on');
   xlabel(i_disptext(g_grind.xaxis.var));
   ylabel(i_disptext(g_grind.yaxis.var));
   zlabel(i_disptext(g_grind.zaxis.var));
   set(gca,'XLim', g_grind.xaxis.lim);
   set(gca,'YLim', g_grind.yaxis.lim);
   set(gca,'ZLim', g_grind.zaxis.lim);
   if g_grind.statevars.dim > 3
      title(i_disptext(['Valid for ' i_othervars(N0, iX.no, iY.no, iZ.no)]));
   end;
   i_setparvar(iZ, N0, oldZ);
   if oldhold
      hold on;
   end;
catch err
   %   err=lasterror;
   N0 = i_setparvar(iZ, N0, oldZ);
   i_keep(N0);
   rethrow(err);
end;



