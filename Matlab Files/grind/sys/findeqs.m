%FINDEQS   Find all equilibria
%   Try to find all equilibria (may fail or find double equilibria)
%   It uses the simplex optimizing method and tries various random initial conditions.
%
%   Usage:
%   FINDEQS - find equibria select one.
%   FINDEQS -r - refreshes the list of equilibria.
%   FINDEQS -a - try to add to the existing list of equilibria.
%   FINDEQS(ntrial,scale) - ntrial is the number of trials (default=100), scale is the 
%   scale of the initial conditions (default=50). 
%   eqlist=FINDEQS - saves the list of equilibria to a cell array.
%
%   See also findeq, eigen

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function eqlist = findeqs(ntrial, scal)
global g_grind;
changed = 0;
eqs = [];
if nargin < 1
   ntrial = 100;
else
   if ischar(ntrial) && ((strncmpi(ntrial, '-r', 2) ||strncmpi(ntrial, '-a', 2)))
         if strncmpi(ntrial, '-a', 2)&&isfield(g_grind,'findeqs')
            for i=1:length(g_grind.findeqs.eqlist)
               eqs(i,:)=str2num(g_grind.findeqs.eqlist{i}(strfind(g_grind.findeqs.eqlist{i})+2:end,':')); %#ok
            end;
         end;
         changed = 1;
         ntrial = 100;
   else
      ntrial = i_checkstr(ntrial);
   end;
end;
if nargin < 2
   scal = 50;
else
   scal = i_checkstr(scal);
end;
settings = getparsettings(ntrial, scal);
if isfield(g_grind, 'findeqs')&&~changed
   changed=(length(settings)~=length(g_grind.findeqs.settings)) || (sum(settings~=g_grind.findeqs.settings) > 0);
else
   changed = 1;
end;
epsil = 1E-3;
N0 = i_initvar;
oldN0 = N0;
try
if ~changed
   eqlst = g_grind.findeqs.eqlist;
else
 %  eqs = [];
   eqs = addeq(eqs, N0, epsil);
   if ntrial > 0
      wb = waitbar(0, 'Searching equilibria...');
      for i = 1:ntrial
         f = rand(size(N0)) < 0.1;
         f2 = rand(size(N0)) < 0.1;
         f3 = rand(size(N0)) < 0.05;
         N0 = rand(size(N0)) * scal;
         N0(f) = 0.001 * N0(f);
         N0(f2) = 100 * N0(f2);
         N0(f3) = -N0(f3);
         eqs = addeq(eqs, N0, epsil);
         waitbar(i / ntrial, wb);
      end;
      close(wb);
   end;
   if ~isempty(eqs)
      [xx, index] = sort(min(eqs, [], 2) + sum(eqs, 2) / 1000);
      eqs = eqs(flipud(index), :); %sort in an order that the trivial / negative equilibria are last
   end;
   eqlst = cell(size(eqs,1),1);
   for i = 1:size(eqs, 1)
      N0 = eqs(i, :)';
      [Jacobian, eigenval] = i_eigen(isempty(g_grind.Jacobian),g_grind.solver.iters,N0);
      [stable, issaddle, isspiral] =  i_stability(eigenval, g_grind.solver.isdiffer);
      if issaddle && ~isspiral
         s1 = 'saddle';
      else
         if stable
            s1 = 'stable ';
         else
            s1 = 'unstable ';
         end;
         if isspiral
            s1 = sprintf('%sspiral',s1);
         else
            s1 = sprintf('%snode',s1);
         end;
      end;
      if max(abs(N0)) < 1E-3
         s = '0';
      elseif abs(min(N0)) < 1E-3
         s = '+/0';
      elseif min(N0) < 0
         s = '--';
      else
         s = '++';
      end
      s2 = sprintf('%g ', round(N0 * 1000) / 1000);
      eqlst{i} = sprintf('%s (%s): %s', s1,s, s2);
   end;
end;
if ntrial>0
   g_grind.findeqs.eqlist = eqlst;
   g_grind.findeqs.settings = settings;
end;
i_keep(oldN0);
catch err
%   err=lasterror;
   i_keep(oldN0);
   rethrow(err);
end;
if nargout == 1
   eqlist = eqlst;
else
   % fprintf('Found %d equilibria (might not be all!):\n', size(eqlst, 1));
   [Ans, ok] = listdlg( 'ListString', eqlst , ...
      'SelectionMode' ,'multiple' ,'ListSize' , [400 300],...
      'Name','findeqs',       'PromptString'  ,...
      sprintf('Found %d equilibria (might not be all!), select one:', length(eqlst)));
   if ok
      for i=1:length(Ans)
          s = eqlst{Ans(i)};
          N1 = str2num(s(strfind(s,':') + 2:end)); %#ok
         i_keep(N1);
         findeq;
      end;
   end;
   
end;

function eqs = addeq(eqs, N0, epsil)
[N1, found] = findeq(0, N0);
if found && (isempty(eqs) ||  min(sum((eqs - repmat(N1(:)', size(eqs, 1), 1)).^2, 2)) > epsil)
   eqs(end + 1, :) = N1(:)';
end;

function settings = getparsettings(ntrial, scal)
global g_grind;
if g_grind.statevars.vector
   maxn = 500;
else
   maxn = length(g_grind.pars) + 1;
end;
settings = zeros(1, maxn);
j = 0;
for i = 1:length(g_grind.pars)
   p = evalin('base', char(g_grind.pars{i}));
   jj = numel(p);
   j2 = j + jj;
   while j2 + 1 > maxn
      addn = j2 + 1;
      maxn = maxn + addn;
      settings = [settings, zeros(1,addn)]; %#ok
   end;
   settings(j + 1:j2) = p(1:jj);
   j = j2;
end;
settings(j + 1) = ntrial;
settings(j + 2) = scal;
settings = settings(1:j + 2);

