%FINDEQ   Find nearest equilibrium
%   Find the closest initial condition with zero growth. FINDEQ
%   uses the simplex optimizing method to find an initial
%   condition with zero growth.
%
%   Usage:
%   FINDEQ - normal use of the FINDEQ command.
%   [N1, isfound]=FINDEQ(1) - internal use of FINDEQ, return the initial
%   conditions and a flag if an equilibrium is found. The argument (1)
%   suppresses output.
%
%   See also eigen, perturb, null, findeqs

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function [N1, found] = findeq(display,N0)
global g_Y g_grind;
if nargin == 0
   display = 1;
else
   display = i_checkstr(display);
end;
if nargin<2
  i_parcheck;
  N0 = i_initvar;
end;
if g_grind.solver.nonautonomous&&(display~=1234)
   if strcmp(questdlg('Might not find correct equilibria in nonautonomous equation, continue?','Error','Yes','No','No'),'No')
      error('GRIND:findeq:NonAutonomous','Make equation autonomous by adding tau''=1 and replacing all t by tau');
   end;
end;
opt = optimset('fminsearch');
if ~isfield(g_grind,'findeq')
   g_grind.findeq.TolFun = 1E-6;
   g_grind.findeq.TolX = 1E-10;
   g_grind.findeq.MaxFunEvals = 50000;
   g_grind.findeq.MaxIter = 50000;
end;
opt.TolFun = g_grind.findeq.TolFun;
opt.TolX = g_grind.findeq.TolFun;
opt.MaxFunEvals = length(N0) * g_grind.findeq.MaxFunEvals;
opt.MaxIter = length(N0) * g_grind.findeq.MaxIter;
if display
   opt=optimset(opt,'Display','on');
else
   opt=optimset(opt,'Display','off');
end;
%opt.Display='iter';
if g_grind.solver.isdiffer
   [NRes, res, eqfound] = fminsearch(str2func('i_totgrowth2'), N0, opt);
else
   [NRes, res, eqfound] = fminsearch(str2func('i_totgrowth'), N0, opt);
end;
eqfound=(res<opt.TolFun*2) & eqfound; %No convert
if eqfound
   g_Y = [];
   i_keep(NRes);
   if display
      disp('Equilibrium found:');
      [Jacobian, eigenval] = i_eigen(isempty(g_grind.Jacobian),g_grind.solver.iters);
      [stable,issaddle]= i_stability(eigenval,g_grind.solver.isdiffer);
      val;
      if stable
         where('sblink','fo')
      else
         if issaddle
            where('sblink','go')
         else
            where('sblink','bo');
         end;
      end;
   end;
elseif display
    errordlg('Could not find an equilibrium, try other initial conditions'); 
end
if nargout > 0
   N1 = NRes;
   found = eqfound;
end;

