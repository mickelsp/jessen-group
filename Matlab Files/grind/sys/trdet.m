%TRDET   Plot of trace and deteminant of the Jacobian 
%  The trace and determinant of the Jacobian matrix determine the stability of equilibria. 
%  Plotting them help to determine the eigenvalues of difference and 
%  differential equations. This function works only for 2 dimensional equations.
%
%  See also eigen, findeq

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function trdet(donumerical)
global g_grind;
i_parcheck;
if nargin == 0
   donumerical = 0;
else
   donumerical = i_checkstr(donumerical);
end;
niters = g_grind.solver.iters;
ud = get(get(0,'currentfigure'), 'userdata');
if isfield(ud, 'iters')
   niters = ud.iters;
end;
N0 = i_initvar;
[Jacobian, eigenvalues] = i_eigen(donumerical,niters,N0);
if length(eigenvalues) == 2
   tr = trace(Jacobian);
   determ = det(Jacobian);
   disp('Jacobian (J) =');
   disp(Jacobian);
   fprintf('trace(J) = %g\n\n', tr);
   fprintf('det(J)   = %g\n', determ);
   h = i_makefig('trdet');
   i_plotdefaults(h); 
   hold off;
   plot(tr, determ, 'r+', 'markersize',10);

   % daspect([1 1 1]);
   %  oldhold = ishold;
   hold on;
   set(h,'Name','Trace and determinant of Jacobian (J)');
   title('Trace and determinant of Jacobian (J)');
   xlabel('trace(J)');
   ylabel('det(J)');
   H1 = get(h, 'CurrentAxes');
   max1 = max(abs([tr * 1.1 determ * 1.1]));
   if max1 < 0.18
      lim = [ - 0.2 0.2];
   elseif max1 < 0.35
      lim = [ - 0.4 0.4];
   elseif max1 < 0.45
      lim = [ - 0.5 0.5];
   elseif max1 < 0.95
      lim = [ - 1 1];
   elseif max1 < 1.4
      lim = [ - 1.5, 1.5];
   elseif max1 < 10
      lim = [ - (round(max1) + 1), round(max1) + 1];
   else
      lim = [ - max1, max1];
   end;
   set(H1, 'Xlim', lim);
   set(H1, 'Ylim', lim);
   x = lim(1):(lim(2) - lim(1)) / 100:lim(2);
   if g_grind.solver.isdiffer
      plot(x, 0.25 * x.^2,':k');
      x = [lim(1) lim(2)];
      plot(x, x - 1,':r');
      plot(x,  -x - 1,':r');
      x = [-2 2 0  -2];
      y = [1 1  -1 1];
      plot(x, y,'-k');
   else
      plot([lim(1),lim(2)],[0,0], 'k');
      plot([0,0],[lim(1),lim(2)], 'k');
      plot(x, 0.25 * x.^2,':k');
   end;
   Nres1 = N0;
   for k = 1:niters
      Nres1 = feval(g_grind.odefile, 1, Nres1);
   end;
   if g_grind.solver.isdiffer
      Nres1 = Nres1 - N0;
   end;
   if sum(Nres1) > 0.001
      i_warningdlg('GRIND:trdet:noequilibrium','The model is not in an equilibrium\nUse <a href="matlab:findeq">findeq</a> to find the equilibrium');
   end;
else
   error('GRIND:trdet:No2dim','This function can only be applied for two dimensional systems');
end;

