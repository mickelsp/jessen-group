%ATTRBASIN   Find basins of attraction by simulation
%   Vary the initial states of two state variables and run till an attractor
%   is reached. Plot the mean state ifn the attractor. Used to visualize the 
%   sensitive dependence of the initial state of systems with transient chaos. 
%   Use null and perturb or manifolds to find a simple separatrix.
%
%   Usage:
%   ATTRBASIN - creates a 20x20 plot with the equilibrium state of the x-axis
%   ATTRBASIN N VAR1 - creates a NxN plot with the equilibrium state of variable VAR.
%   RES=ATTRBASIN(N,VAR1) - save the simulation results in the 3 dimensional matrix RES.
%   ATTRBASIN N VAR1 RES - recreate the plot with the result matrix RES.
%
%
%   See also ax, null, perturb, manifolds

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:25 $
function EndsVar = attrbasin(n, var2, EndsX)
global g_Y t g_grind g_t;
i_parcheck;
if nargin < 1
   n = 20;
else
   n = i_checkstr(n);
end;
if nargin < 2
   var2 = g_grind.xaxis.var;
end;
Xs = zeros(n);
Ys = zeros(n);
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
oldX = i_getparvar(iX);
oldY = i_getparvar(iY);
try
   iRes = i_getno(var2);
   if nargin < 3
      EndsX = zeros(n, n, g_grind.statevars.dim);
      if ~g_grind.solver.isdiffer
         stabstep = 2;
      else
         stabstep = NaN;
      end;
      oldstep = g_grind.tstep;
      wb = waitbar(0, 'Calculating...');
      set(wb,'name','attrbasin')
      N0 = i_initvar;
      for j = 1:n
         waitbar(j / n, wb)
         px = g_grind.xaxis.lim(1) + (g_grind.xaxis.lim(2) - g_grind.xaxis.lim(1)) * j / n;
         N0 = i_setparvar(iX, N0, px);
         for k = 1:n
            py = g_grind.yaxis.lim(1) + (g_grind.yaxis.lim(2) - g_grind.yaxis.lim(1)) * k / n;
            N0 = i_setparvar(iY, N0, py);
            Xs(k, j) = px;
            Ys(k, j) = py;
            t1 = t;
            g_grind.tstep = stabstep;
            i_ru(g_grind.odefile,t1, g_grind.ndays, N0, 0);
            N1 = transpose(g_Y(size(g_Y, 1), :));
            t1 = g_t(size(g_t, 1));
            g_grind.tstep = oldstep;
            i_ru(g_grind.odefile,t1, g_grind.ndays / 10, N1, 0);
            mY = mean(g_Y);
            for l = 1:g_grind.statevars.dim
               EndsX(k, j, l) = mY(l);
            end;
         end;
      end;
      delete(wb);
   else
      N0 = i_initvar;
      for j = 1:n
         px = g_grind.xaxis.lim(1) + (g_grind.xaxis.lim(2) - g_grind.xaxis.lim(1)) * j / n;
         N0 = i_setparvar(iX, N0, px);
         for k = 1:n
            py = g_grind.yaxis.lim(1) + (g_grind.yaxis.lim(2) - g_grind.yaxis.lim(1)) * k / n;
            N0 = i_setparvar(iY, N0, py);
            Xs(k, j) = px;
            Ys(k, j) = py;
         end;
      end;
   end;
   colmap = 'i_grindcolormap';
   H = i_makefig('attrbasin');
   set(H,'name','Basin of attraction');
   pcolor(Xs, Ys, EndsX(:, :, iRes.no));
   ylabel(['initial ' g_grind.yaxis.var]);
   xlabel(['initial ' g_grind.xaxis.var]);
   shading flat;
   colormap(colmap);
   cb = colorbar;
   ylab = get(cb, 'ylabel');
   set(ylab,'string',['Mean equilibrium ' var2]);
   if nargout > 0
      EndsVar = EndsX;
   end;
   N0 = i_setparvar(iX, N0, oldX);
   i_setparvar(iY, N0, oldY);
catch err
%   err=lasterror;
   N0 = i_setparvar(iX, N0, oldX);
   i_setparvar(iY, N0, oldY);
   rethrow(err);
end;
