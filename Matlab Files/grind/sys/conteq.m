%CONTEQ   Continue a stable or unstable equilibrium 
%   This command can do a simple 1D bifurcation analysis. First select an equilibrium 
%   with findeq, then you can continue the stable or unstable equilibrium. It can detect a 
%   fold (=saddle-node), transcritical and a Hopf bifurcation. It cannot distinguish 
%   between a transcritical and fold bifurcation. A plot of the equilibrium and the 
%   eigenvalues is made
%
%   Usage:
%   CONTEQ -  - user is prompted for information
%   CONTEQ APAR TOVAL - continue APAR from the current value till TOVAL. Take 200 steps.
%   CONTEQ APAR TOVAL NSTEP - continue APAR from the current value till TOVAL. 
%   Take NSTEP steps.
%
%   See also findeq, paranal, conteq2d

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function conteq(apar, nend, nsteps, silent)
evalin('base','global g_cont;');
global g_grind g_cont g_paranal;
if isfield(g_grind, 'pars')&&isempty(g_grind.pars)
   error('GRIND:conteq:NoParameters','No parameters to analyse');
end;
if nargin < 4
   silent = 0;
end;
i_parcheck;
g_cont.Hopf = [];
g_cont.Fold = [];
g_cont.DoubleEig = [];
hasperm = ~isempty(g_grind.permanent);
if hasperm
   Nperm = defpermanent('-p', 1);
end;
nend2 = [];
if (nargin == 0)
   answer = i_conteqdlg;
   if isempty(answer)
      return;
   end;
   apar = g_grind.pars{answer.ipar};
   nend = answer.min;
   nend2 = answer.max;
   if isempty(nend) && ~isempty(nend2)
      nend = nend2;
      nend2 = [];
   end;
   nsteps = answer.nsteps;
   %   N0 = answer.N0;
   start = evalin('base', apar);
   if isempty(g_grind.continue.N0)
      for i = 1:length(g_grind.continue.eqlist) - 1
         s = g_grind.continue.eqlist{i};
         N1 = str2num(s(strfind(s, ':') + 2:end)); %#ok
         i_keep(N1);
         if ~isempty(nend2)
            conteq(apar, [nend nend2], nsteps, silent);
         else
            conteq(apar, nend, nsteps, silent);
         end;
      end;
      return;
   end;
elseif (nargin==1) && strcmp(apar, '-r')
   if isfield(g_grind, 'continue')
      ud = g_grind.continue;
      apar = g_grind.pars{ud.ipar};
      nend = ud.min;
      nend2 = ud.max;
      nsteps = ud.nsteps;
      start = evalin('base', apar);
   else
      conteq;
      return;
   end;
else
   start = evalin('base', apar);
   if nargin < 2
      if start ~= 0
         nend = start + start * 2;
      else
         nend = 1;
      end;
   else
      nend  = i_checkstr(nend);
      if length(nend) == 2
         nend2 = nend(2);
         nend = nend(1);
      end;
   end;
   if nargin  < 3
      nsteps = 200;
   else
      nsteps  = i_checkstr(nsteps);
   end;
end;
g_cont.parname = apar;
if start==nend && ~isempty(nend2) && start~=nend2
   nend = nend2;
   fprintf('Nothing to do: start value of %s equals end value (%g)\n', apar, start);
   nend2 = [];
elseif ~isempty(nend2)
   if start < nend && start < nend2
      nend = max(nend, nend2);
      nsteps = nsteps * 2;
      nend2 = [];
   end;
   if start > nend && start > nend2
      nend = min(nend, nend2);
      nsteps = nsteps * 2;
      nend2 = [];
   end;
   if ~isempty(nend2)
      ntot = nsteps;
      nsteps = ceil(ntot * abs((start - nend) / (nend2 - nend))) + 1;
      nsteps2 = ceil(ntot * abs((start - nend2) / (nend2 - nend))) + 1;
   end;
end;
if start ~= nend
   [cont, foldexpected] = continueeq(apar, start, nend, nsteps);
   g_cont.Y = cont.Y;
   g_cont.eig = cont.eig;
   g_cont.p = cont.p;
   istart = 1;
   dim = size(g_cont.Y, 2);
   if ~isempty(nend2);
      if start ~= nend2
         %     disp(sprintf('Nothing to do: start value of %s equals end value (%g)', apar, start));
         %   else
         [cont1, foldexpected1] = continueeq(apar, start, nend2, nsteps2);
         foldexpected = foldexpected || foldexpected1;
         if nend < nend2
            istart = length(g_cont.p);
            g_cont.Y = [flipud(g_cont.Y); cont1.Y];
            g_cont.eig = [flipud(g_cont.eig); cont1.eig];
            g_cont.p = [flipud(g_cont.p); cont1.p];
         else
            istart = length(cont1.p);
            g_cont.Y = [flipud(cont1.Y); g_cont.Y];
            g_cont.eig = [flipud(cont1.eig); g_cont.eig];
            g_cont.p = [flipud(cont1.p); g_cont.p];
         end;
      end;
   end;
   if hasperm
      g_cont.perm = repmat(Nperm', length(g_cont.p), 1);
   end;
   [g_cont.isstable, g_cont.issaddle, g_cont.isspiral] =  i_stability(g_cont.eig', g_grind.solver.isdiffer);
   if ~silent
      i_makefig('conteq');
      hold on;
      if length(g_cont.p) == length(g_cont.eig)
         plot3(g_cont.p, real(g_cont.eig(:,1:dim)),imag(g_cont.eig(:,1: dim)));
         i_plotdefaults;
         set(gca,'DrawMode', 'fast');
         xlabel(apar);
         ylabel('real eigenvalue');
         zlabel('imaginary eigenvalue');
         l = get(legend, 'userdata');
         if ~isempty(l)
            oldleg = l.lstrings;
         else
            oldleg = {};
         end;
      else
          i_warningdlg('GRIND:conteq:cannotplot','Problem with dimensions, cannot plot eigenvalues');
          oldleg={};
      end;
   end;
   leg = cell(1, dim);
   for i = 1:dim
      leg{i} = sprintf('eigen(%d)', i);
   end;
   
   n = length(g_cont.p);
   ifold = zeros(20, 1);
   nfold = 0;
   k = dim + 1;
   for i = 1:n - 2
      if (g_cont.isstable(i) ~= g_cont.isstable(i+1)) && g_cont.isspiral(i)
         fprintf('Hopf bifurcation between %s=%0.5g and %0.5g\n', apar, g_cont.p(i), g_cont.p(i + 1));
         plotbif(i, i + 1, g_cont.p, g_cont.eig, dim, 'ko',silent);
         g_cont.Hopf = [g_cont.Hopf; (g_cont.p(i) + g_cont.p(i + 1)) / 2];
         leg{k} = 'Hopf';
         k = k + 1;
      end;
      if (g_cont.isspiral(i) ~= g_cont.isspiral(i + 1))
         fprintf('Double eigenvalues between %s=%0.5g and %0.5g\n', apar, g_cont.p(i), g_cont.p(i + 1));
         plotbif(i, i + 1, g_cont.p, g_cont.eig, dim, 'kx',silent);
         g_cont.DoubleEig = [g_cont.DoubleEig ; (g_cont.p(i) + g_cont.p(i + 1)) / 2];
         leg{k} = 'Double eigenvalues';
         k = k + 1;
      end;
      if (g_cont.isstable(i) ~= g_cont.isstable(i+1)) && ~g_cont.isspiral(i)
         fprintf('Fold or transcritical bifurcation between %s=%0.5g and %0.5g\n', apar, g_cont.p(i), g_cont.p(i + 1));
         plotbif(i, i + 1, g_cont.p, g_cont.eig, dim, 'k^',silent);
         g_cont.Fold = [g_cont.Fold ; (g_cont.p(i) + g_cont.p(i + 1)) / 2];
         nfold = nfold + 1;
         ifold(nfold) = i + 1;
         leg{k} = 'Fold or transcritical';
         k = k + 1;
      end;
   end;
   if foldexpected&&(n > 1)
      plotbif(n - 1, n - 1, g_cont.p, g_cont.eig,dim, 'k^',silent);
      i = n - 1;
      if i < 1
         i = 1;
      end;
      g_cont.Fold = [g_cont.Fold; g_cont.p(i)];
      leg{k} = 'Fold';
   end;
   if ~silent
      if ~isempty(oldleg)
         legend({oldleg{:} leg{:}});  %#ok
      else
         legend(leg{:});
      end;
      if ~exist('g_paranal','var')
         g_paranal = [];
      end;
      oldparanal = g_paranal;
      g_paranal = g_cont;
      g_paranal.t = zeros(size(g_cont.p));
      if ~isfield(g_grind, 'paranal')||~isfield(g_grind.paranal, 'dlg')
         paranal('-defaults');
         g_grind.paranal.dlg = i_paranaldialog('initstruct', 1);
         g_grind.paranal.dlg.par{1} = apar;
         g_grind.paranal.dlg.nend(1) = nend;
         if isempty(nend2)
            g_grind.paranal.dlg.start(1) = start;
         else
            g_grind.paranal.dlg.start(1) = nend2;
         end;
      end;
      if nfold > 0
         i = 1;
         while (ifold(i) < istart) && (i < nfold)
            i = i + 1;
         end;
         if i > 1
            sel = ifold(i - 1):ifold(i);
         elseif ifold(i) < istart
            sel = ifold(i):length(g_cont.p);
         else
            sel = 1:ifold(i);
         end;
         g_paranal.Y = g_paranal.Y(sel, :);
         g_paranal.p = g_paranal.p(sel);
         g_paranal.t = g_paranal.t(sel);
         if hasperm && ~isempty(g_paranal.perm)
            g_paranal.perm = g_paranal.perm(sel);
         end;
      end;
      oldlines = g_grind.paranal.dlg.lines;
      g_grind.paranal.dlg.lines = 1;
      if g_cont.isstable(istart)
         g_grind.pen.pen = '-';
      else
         g_grind.pen.pen = ':';
      end;
      i_paranal(1);
      g_grind.paranal.dlg.lines = oldlines;
      g_grind.pen.pen = '-';
      g_paranal = oldparanal;
      %      [h, isnew] = i_makefig('paranal');
      %      iX = i_getno(g_grind.xaxis.var);
      %      iY = i_getno(g_grind.yaxis.var);
      %      if isnew
      %         xlabel(i_disptext(apar));
      %         ylabel(i_disptext(g_grind.xaxis.var));
      %         zlabel(i_disptext(g_grind.yaxis.var));
      %         set(findobj(gcf,'type','line'),'linewidth',g_grind.pen.linewidth);
      %         i_plotdefaults(gcf);
      %      end;
      %      set(h, 'Name', 'Parameter analysis');
      %      hold on;
      %      i_plotdefaults;
      %      if dim > 1
      %         if g_cont.isstable(istart)
      %            s = '-';
      %         else
      %            s = ':';
      %         end;
      %         [ydata,zdata] =  i_conteqsimpfun(iX, iY, g_cont.Y(sel,:));
      %         plot3(g_cont.p(sel), ydata, zdata,s);
      %         set(gca,'DrawMode', 'fast');
      %     else
      %         ydata =  i_conteqsimpfun(iX, iY, g_cont.Y(sel,:));
      %         plot(g_cont.p(sel), ydata);
      %      end;
   end;
else
   fprintf('Nothing to do: start value of %s equals end value (%g)\n', apar, start);
end;

% function [ydata, zdata] = i_conteqsimpfun(iX, iY, YY)
% global g_Y g_t g_grind;
% oldY = g_Y;
% oldt = g_t;
% g_Y = YY;
% g_t = zeros(size(YY, 1), 1);
% if isempty(iX.no) & ~isempty(g_grind.xaxis.var)
%    ydata =  i_getoutfun(g_grind.xaxis.var);
% else
%    ydata = g_Y(:, iX.no);
% end;
% if isempty(g_grind.yaxis.var)
%    zdata = [];
% elseif isempty(iY.no)
%    zdata  = i_getoutfun(g_grind.yaxis.var);
% else
%    zdata = g_Y(:, iY.no);
% end;
% g_Y = oldY;
% g_t = oldt;

function [cont, foldexpected] = continueeq(apar, start, nend, nsteps)
global g_grind;
oldN0 = i_initvar;
disturb = (nend - start) / nsteps;
[N0, eqfound] = findeq(0);
if ~eqfound
   error('GRIND:conteq:NoEquilibrium','Cannot find equilibrium');
end;
%nmax = 200;
dim = length(N0);
cont.Y = NaN * ones(nsteps, dim);
cont.eig = NaN * ones(nsteps, dim);
cont.p = zeros(nsteps, 1);
p = evalin('base', apar);
oldpar = p;
%[N0, eqfound] = findeq(0);

foldexpected = 0;
[Jacobian, eigenval] = i_eigen(1,g_grind.solver.iters);
[isstable1, issaddle1, isspiral1] =  i_stability(eigenval, g_grind.solver.isdiffer);
s = 'Continuing ';
if issaddle1 && ~isspiral1
   s = [s 'a saddle equilibrium from %s=%0.5g till %0.5g\n'];
else
   if isstable1
      s = [s 'a stable '];
   else
      s = [s 'an unstable '];
   end;
   if isspiral1
      s = [s 'spiral from %s=%0.5g till %0.5g\n'];
   else
      s = [s 'node from %s=%0.5g till %0.5g\n'];
   end;
end;
fprintf(s, apar, start, nend);
cont.Y(1, :) = N0.';
N1 = N0;
cont.eig(1, :) = eigenval';
cont.p(1) = p;
n = 2;
while eqfound && (n < nsteps)
   p = p + disturb;
   multassignin('base', apar, p);
   [N0, eqfound] = findeq(0);
   if eqfound
      [Jacobian, eigenval] = i_eigen(1,g_grind.solver.iters);
      cont.Y(n, :) = N0.';
      cont.eig(n, :) = eigenval';
      cont.p(n) = p;
      n = n + 1;
      i_keep(2 * N0 - N1);
      N1 = N0;
   else
      foldexpected = 1;
      fprintf('Cannot find equilibrium anymore at %s = %0.5g, Fold bifurcation?\n',apar,p);
   end;
end;
multassignin('base', apar, oldpar);
i_keep(oldN0);
cont.Y = cont.Y(1:n - 1, :);
cont.eig = cont.eig(1:n - 1, :);
cont.p = cont.p(1:n - 1);

function plotbif(i, i2, Ps, Ns, dim, symb,silent)
if ~silent&&(length(Ps)>=i)&&(length(Ps)>=i2)
   x = ones(1, dim) * (Ps(i) + Ps(i2)) / 2;
   y = real(Ns(i:i2, 1:dim));
   z = imag(Ns(i:i2, 1:dim));
   if i2 ~= i
      y = mean(y);
      z = mean(z);
   end;
   plot3(x',y',z',symb);
end;




