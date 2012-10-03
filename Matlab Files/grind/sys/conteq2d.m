%CONTEQ2D   Continue an equilibrium in two dimensions
%   This command can do a simple 2D bifurcation analysis. First select an equilibrium 
%   with findeq, then you can continue the stable or unstable equilibrium. 
%   It can detect a fold (=saddle-node), transcritical and a Hopf bifurcation with 
%   respect to two parameters. It cannot distinguish between a transcritical and fold bifurcation.
%   A 2d plot of the bifurcations is made. It can be that only a part of the plot is
%   found (change initital conditions to find the other bifurcations).  
%
%   Usage:
%   CONTEQ2D -  - user is prompted for information
%
%   See also findeq, paranal, conteq

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function conteq2d(Opt)
i_parcheck;
evalin('base','global g_cont2d;');
global g_grind g_cont g_cont2d;
if isfield(g_grind,'pars')&&isempty(g_grind.pars)
   error('GRIND:conteq2d:NoParameters','no parameters to analyse');
end;
if nargin == 0
   if isfield(g_grind, 'conteq2d') && ~isempty(g_grind.conteq2d)
      Opt = g_grind.conteq2d;
      answer{1} = Opt.par{1};
      answer{2} = num2str(Opt.min(1));
      answer{3} = num2str(Opt.max(1));
      answer{4} = num2str(Opt.nsteps);
      answer{5} = Opt.par{2};
      answer{6} = num2str(Opt.min(2));
      answer{7} = num2str(Opt.max(2));
      answer{8} = num2str(Opt.ngrid);
   else
      answer={'','0','10','300','','0','10','10'};
   end;
   prompt={'parameter xaxis','Minimum','Maximum','Number of steps (conteq)','Parameter yaxis','Minimum','Maximum','Number of steps (grid)'};
   answer = inputdlg(prompt, 'conteq 2D', 1, answer);
   Opt.par{1} = answer{1};
   Opt.min(1) = str2double(answer{2});
   Opt.max(1) = str2double(answer{3});
   Opt.nsteps = str2double(answer{4});
   Opt.par{2} = answer{5};
   Opt.min(2) = str2double(answer{6});
   Opt.max(2) = str2double(answer{7});
   Opt.ngrid = str2double(answer{8});
else
   if ischar(Opt) && strncmp(Opt,'-p',2)
      plotresults;
      return;
   end;    
end;
g_grind.conteq2d = Opt;
par1 = Opt.par{1};
par2 = Opt.par{2};
conteq(par1, [Opt.min(1), Opt.max(1)], Opt.nsteps,1);
vcont = g_cont;
n = length(vcont.p);
oldp1 = evalin('base', par1);
oldp2 = evalin('base', par2);
try
   g_cont2d.Hopfs = [];
   g_cont2d.pars{1}=par1;
   g_cont2d.pars{2}=par2; 
   g_cont2d.Doubles = [];
   g_cont2d.Folds = [];
   for i = 1:Opt.ngrid
      disp(i);
      ip = round(n * (i - 1) / (Opt.ngrid - 1) + 1);
      if ip > n
         ip = n;
      end;
      p = vcont.p(ip);
      N1 = vcont.Y(ip, :)';
      i_keep(N1);
      assignin('base', par1, p);
      conteq(par2, [Opt.min(2), Opt.max(2)], Opt.nsteps,1);
      g_cont2d.Hopfs = [g_cont2d.Hopfs; [ones(size(g_cont.Hopf)) * p, g_cont.Hopf]];
      g_cont2d.Doubles = [g_cont2d.Doubles; [ones(size(g_cont.DoubleEig)) * p, g_cont.DoubleEig]];
      g_cont2d.Folds = [g_cont2d.Folds; [ones(size(g_cont.Fold)) * p, g_cont.Fold]];
   end;
   plotresults;
   assignin('base', par1, oldp1);
   assignin('base', par2, oldp2);
catch err
   assignin('base', par1, oldp1);
   assignin('base', par2, oldp2);
   rethrow(err);
end
function plotresults
global g_cont2d;
i_makefig('conteq2d');
   hold on;
   leg = {};
   if ~isempty(g_cont2d.Hopfs)
      plot(g_cont2d.Hopfs(:, 1), g_cont2d.Hopfs(:, 2), 'or');
      leg = {'Hopf'};
   end;
   hold on;
   if ~isempty(g_cont2d.Doubles)
      plot(g_cont2d.Doubles(:, 1), g_cont2d.Doubles(:, 2), 'xb');
      leg = [leg, {'Double eigenvalues'}];
   end;
   if ~isempty(g_cont2d.Folds)
      plot(g_cont2d.Folds(:, 1), g_cont2d.Folds(:, 2), '+k');
      leg = [leg, {'Fold/transcritical'}];
   end;
   leg=i_addlegend(leg);
   if length(leg)>1
      legend(leg);
   elseif length(leg)==1
      legend(leg{1});
   end;
   i_plotdefaults;
   xlabel(g_cont2d.pars{1});
   ylabel(g_cont2d.pars{2});
   hold off;
