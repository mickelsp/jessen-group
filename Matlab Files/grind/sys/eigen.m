%EIGEN   Calculate eigenvalues of currently selected equilibrium
%   Display eigenvalues and eigenvectors of the current
%   initial conditions. The eigenvalues are plotted in
%   the complex plane.
%
%   Usage:
%   EIGEN - Calculate eigenvalues with analytical Jacobian
%   (see enterjac). If there is no Jacobian entered, the Jacobian
%   matrix is approximated numerically. 
%   EIGEN -numerical (-n) Jacobian matrix is approximated numerically, even
%          if an analytical Jacobian is entered.
%   EIGEN -analytical (-a) Jacobian matrix is calculated numerically, only
%          if an analytical Jacobian is entered.
%   Eigs=EIGEN - calculates all eigenvalues during the last run.
%   Eigs=EIGEN('-paranal') - calculates all eigenvalues from the last paranal run.
%
%
%   See also findeq, enterjac, paranal, trdet

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function Eigs=eigen(opt)
global g_grind g_Y t;
i_parcheck;
if g_grind.solver.nonautonomous
   if strcmp(questdlg('Might not find correct eigenvalues in nonautonomous equation, continue?','Error','Yes','No','No'),'No')
      error('GRIND:eigen:NonAutonomous','Make equation autonomous by adding tau''=1 and replacing all t by tau');
   end;
end;
doparanal=0;
donumerical = isempty(g_grind.Jacobian);
if (nargin==1)
   if strncmpi(opt,'-n',2)
      donumerical=1;
   elseif strncmpi(opt,'-a',2)
      donumerical=0;
   elseif strncmpi(opt,'-p',2)
      doparanal=1;
   else
      donumerical = i_checkstr(opt);
   end;
end
niters = g_grind.solver.iters;
ud = get(get(0,'currentfigure'), 'userdata');
if isfield(ud, 'iters')
   niters = ud.iters;
end;
if nargout>0 
   if ~doparanal
      N0 = i_initvar;
      if isempty(g_Y)
         i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
      end;
      Eigs=zeros(size(g_Y));
      for i=1:size(g_Y,1)
         N0=g_Y(i,:)';
         [J,eigenvalues]=i_eigen(donumerical,niters,N0);
         Eigs(i,:)=eigenvalues';
      end;
   else
      Eigs=outfun('eigen','-p');
   end;
   return;   
elseif donumerical
   disp('Jacobian numerically approximated');
else
   disp('Uses Jacobian entered by user');
end;
N0 = i_initvar;
[Jacobian, eigenvalues, eigenvect] = i_eigen(donumerical,niters,N0);
%eigenvalues = diag(eigenvalues);
disp('Jacobian:');disp(Jacobian);
disp('Eigen vectors:');disp(eigenvect);
disp('Eigenvalues:');disp(eigenvalues);
%eigenvect
%eigenvalues
%N0 = i_initvar;
Nres1 = N0;
for k = 1:niters
   Nres1 = feval(g_grind.odefile, 1, Nres1);
end;
if g_grind.solver.isdiffer
   Nres1 = Nres1 - N0;
end;
h = i_makefig('eigen');
plot(real(eigenvalues),imag(eigenvalues),'rx','MarkerSize',12);
if sum(Nres1) > 0.001  
   i_warningdlg('GRIND:eigen:noequilibrium', 'The model is not in an equilibrium\nUse <a href="matlab:findeq">findeq</a> to find the equilibrium');
end;
daspect([1 1 1]);
oldhold = ishold;
hold on;
if g_grind.solver.isdiffer
   x = -pi:.1:pi + 0.1;
   plot(sin(x), cos(x), 'k');
end;
set(h,'Name','Plot of eigenvalues in the complex plane');
xlabel('Real part of eigenvalue');
ylabel('Imaginary part of eigenvalue');
H1 = get(h, 'CurrentAxes');
if g_grind.solver.isdiffer
   max1=max(abs([get(H1, 'Xlim') get(H1, 'Ylim')]));
   lim = [ - max1, max1];
else
   max1 = max(abs([real(eigenvalues); imag(eigenvalues)]));
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
      max1=max(abs([get(H1, 'Xlim') get(H1, 'Ylim')]));
      lim = [ - max1, max1];
   end;
end;
set(H1, 'Xlim', lim);
set(H1, 'Ylim', lim);
plot([complex(0, lim(1)), complex(0, lim(2))], 'k');
plot([complex(lim(1), 0), complex(lim(2), 0)], 'k');
if ~oldhold
   hold off;
end
