%OUTFUN   Make matrix with a variable or an equation
%   Make a matrix from the last run with any auxiliary/state variables from the model.
%   Note that there is no rerun of the model if parameters have been changed.
%
%
%   Usage:
%   A=OUTFUN('FUN') - makes a matrix A with the output FUN.
%   A=OUTFUN('FUN','-p') - makes a matrix A with the output FUN using the last 
%      paranal result
%   A=OUTFUN('FUN','-p2') - makes a matrix A with the output FUN using the last 
%      paranal2d result
%   A=OUTFUN('FUN','-c') - makes a matrix A with the output FUN using the last 
%      conteq result
%   A=OUTFUN('FUN','-m') - makes a matrix A with the output FUN using the last 
%      mcarlo results.
%
%
%   Examples:
%   A=OUTFUN('x*y/par') x, y = state variable par = parameter
%   A=OUTFUN('cons') cons = permanent variable
%
%
%   See also out, time, paranal, conteq

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function [F,ndxs] = outfun(s, opt, acatfun)
global g_grind g_paranal g_paranal2d g_cont g_mcarlo;
i_parcheck;
if nargin == 0
   if ~isfield(g_grind, 'outfun')
      g_grind.outfun.fun = '';
      g_grind.outfun.opt = 1;
   end;
   res = i_outfundlg('init', g_grind.outfun);
   if isempty(res)
      return;
   end;
   s = res.fun;
   opt = ['-', res.optstr];
   g_grind.outfun.fun = res.fun;
   g_grind.outfun.opt = res.opt;
elseif nargin == 1
   opt = '-n';
end;
if (nargin < 3)
   acatfun = [];
end;
if isempty(s)
   F = [];
   return;
end;
%ps=[];
s = outf('changeshortcut', s);
if isstruct(opt)
   F = analysestruc(opt, s);
   return;
end;
usemcarlo = strncmpi('-m', opt, 2);
useparanal = strncmpi('-p', opt, 2);
if useparanal && ~isempty(strfind(opt, '2'))
   useparanal = 0;
   useparanal2d = 1;
else
   useparanal2d = 0;
end;
usecont = strncmpi('-c', opt, 2);
if useparanal2d
   if ~isfield(g_paranal2d, 'parname1')
      error('GRIND:outfun:paranal2d','Error: for outfun(''-paranal2d'') it is needed that paranal2d has been run');
   end;
   F = analysestruc(g_paranal2d, s);
   if ~isempty(acatfun)
      [F,ndxs] = catfun(acatfun, F, g_paranal2d.p);
   end;
elseif useparanal
   if ~isfield(g_paranal, 'parname')
      error('GRIND:outfun:paranal','Error: for outfun(''-paranal'') it is needed that paranal has been run');
   end;
   F = analysestruc(g_paranal, s);
   if ~isempty(acatfun)
      [F,ndxs] = catfun(acatfun, F, g_paranal.p);
   end;
elseif usecont
   if ~isfield(g_cont, 'parname')
      error('GRIND:outfun:conteq','Error: for outfun(''-conteq'') it is needed that conteq has been run');
   end;
   F = analysestruc(g_cont, s);
   if ~isempty(acatfun) %not very usefull, but each parameter has only one value, so that is not usefull either
      [F,ndxs] = catfun(acatfun, F);
   end;
elseif usemcarlo
   if ~isempty(g_mcarlo) &&(~isfield(g_mcarlo, 'pars')||~isfield(g_mcarlo, 'Ys'))
      error('GRIND:outfun:mcarlo','Error: for outfun(''-mcarlo'') it is needed that mcarlo has been run');
   end;
   F = analysemcarlo(s);
else
   %   THIS GIVES PROBLEMS WITH PARANAL
   %   N0 = i_initvar;
   %   if i_settingschanged(N0, g_grind.ndays)
   %  disp('running');
   %      i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
   %   end;
   F = i_getoutfun(s);
   if ~isempty(acatfun)
      [F,ndxs] = catfun(acatfun, F);
   end;
end;

function F = analysemcarlo(s)
global g_Y g_t g_mcarlo g_permanent;
oldgY = g_Y;
oldgt = g_t;
oldperm = g_permanent;
[n, npar] = size(g_mcarlo.p);
nt = size(g_mcarlo.Ys, 1);
%pars = cell(npar, 1);
%shortcut for speed
if strcmp(s,'t')
   F = repmat(g_mcarlo.t, g_mcarlo.n, 1);
   return;
end;
if  ~isempty(g_mcarlo.paranalpar)
   par = g_mcarlo.allpars{g_mcarlo.paranalpar}.name;
   if strcmp(s, par) || strcmp(s, '<param1>')
      F = repmat(g_mcarlo.paranal.p, g_mcarlo.n, 1);
      return;
   elseif strcmp(s, '#')
      F = g_mcarlo.paranal.p;
      return;
   else
      %replace parameter with oufun('#','-m') (returns g_mcarlo.paranal.p
      symbs=regexp(s,'\w+','match');
      ipars = strcmp(symbs, par);
      if any(ipars)
         isymbs=regexp(s,'\w+','start');
         isymbs = isymbs(ipars);
         for i = length(isymbs):-1:1
            if length(s) < isymbs(i) + length(par) + 1
               s=[s(1:isymbs(i)-1) 'outfun(''#'',''-m'')'];
            elseif  isymbs(i) < 1
               s=['outfun(''#'',''-m'')' s(isymbs(i)+length(par):end)];
            else
               s=[s(1:isymbs(i)-1) 'outfun(''#'',''-m'')' s(isymbs(i)+length(par):end)];
            end;
         end;
      end;
   end;
end;
try
   F = zeros(nt * n, 1);
   g_t = g_mcarlo.t;
   for i = 1:n
      for k = 1:npar
         multassignin('base', g_mcarlo.pars{k}.name, g_mcarlo.p(i,k));
      end;
      g_Y = g_mcarlo.Ys(:, :, i); %#ok
      if ~isempty(g_mcarlo.perm)
         g_permanent.Y = g_mcarlo.perm(:, :, i);
      end;
      g_permanent.t = g_t;
      FF = i_getoutfun(s);
      if size(FF, 2) > size(F, 2)
         F = repmat(F, 1, size(FF, 2));
      end;
      F((i - 1) * nt + 1:i * nt, :) = FF;
   end;
   for i = 1:length(g_mcarlo.pars)
      if g_mcarlo.pars{i}.selected
         multassignin('base', g_mcarlo.pars{i}.name, g_mcarlo.pars{i}.value);
      end;
   end;
   g_Y = oldgY;
   g_t = oldgt;
   g_permanent  = oldperm;
   clear oldgY oldgt FF oldperm;
catch err
   %   err=lasterror;
   for i = 1:length(g_mcarlo.pars)
      if g_mcarlo.pars{i}.selected
         multassignin('base', g_mcarlo.pars{i}.name, g_mcarlo.pars{i}.value);
      end;
   end;
   g_Y = oldgY; %#ok
   g_t = oldgt; %#ok
   rethrow(err);
end;

function F = analysestruc(strc, s)
global g_Y g_t g_permanent g_grind;
if isfield(strc, 'parname1')
   par1 = char(strc.parname1);
else
   par1 = char(strc.parname);
end;
oldpar1 = evalin('base', par1);
hast = isfield(strc, 't');
haspar2 = isfield(strc, 'parname2');
hasperm = isfield(strc, 'perm');
if hasperm %newest MATLAB version problem
   hasperm =  ~isempty(strc.perm);
end;
if haspar2
   oldpar2 = evalin('base', char(strc.parname2));
end;
%shortcuts for speed
if strcmp(s, 't') && hast
   F = strc.t;
   return;
elseif strcmp(s, '<param1>')
   F = strc.p;
   return;
elseif strcmp(s, '<param2>')
   F = strc.p2;
   return;
end;

g_l_X = i_getno(s);
if g_l_X.ispar
   if strcmp(s, par1)
      F = strc.p;
   elseif haspar2 && strcmp(s, strc.parname2)
      F = strc.p2;
   else
      F = evalin('base', s);
      F = repmat(F(:)',size(strc.p, 1), 1);
   end;
   return;
elseif g_l_X.isvar
   if ~isempty(g_l_X.vecno) % &  isempty(strfind(g_afun, '('))
      F = strc.Y(:, g_grind.statevars.dims{g_l_X.vecno}.from:g_grind.statevars.dims{g_l_X.vecno}.to);
   else
      F = strc.Y(:, g_l_X.no);
   end;
   return;
elseif g_l_X.isperm
   p =  g_grind.permanent{g_l_X.no};
   F = strc.perm(:, p.from:p.to);
   return;
end
oldgY = g_Y;
oldgt = g_t;
try
   if haspar2
      pars = {par1, strc.parname2};
   else
      pars = {par1};
   end;
   if isempty(intersect(symvar(s), pars))&&isempty(strfind(s, '''domeigen'''))...
         &&isempty(strfind(s, '''realeigen'''))&&isempty(strfind(s, '''imageigen'''))
      %about 100x as fast as the loop below
      g_Y = strc.Y; %#ok
      if hast
         g_t = strc.t;
      else
         g_t = zeros(size(strc.p));
      end;
      if hasperm
         g_permanent.Y = strc.perm;
         g_permanent.t = g_t;
      end;
      F = i_getoutfun(s);
      if size(strc.p, 2) > size(F, 2)
         F = repmat(F, 1, size(strc.p, 2));
      end;
   else
      %Slow loop is needed if control parameters are in the loop
      p = strc.p(1);
      if haspar2
         p1 = strc.p1(1);
      end;
      F = zeros(size(strc.p));
      j = 1;
      for i = 1:size(strc.Y, 1)
         if (i==size(strc.Y, 1)) || (p ~= strc.p(i+1)) || (haspar2 && (p1 ~= strc.p1(i+1)))
            multassignin('base', par1, strc.p(i));
            if haspar2
               multassignin('base', char(strc.parname2), strc.p1(i));
            end;
            g_Y = strc.Y(j:i, :); %#ok
            if hast
               g_t = strc.t(j:i, :); %#ok
            else
               g_t = zeros(i - j + 1, 1); %#ok
            end;
            if hasperm
               g_permanent.Y = strc.perm(j:i, :);
               g_permanent.t = strc.t(j:i, :);
            end;
            FF = i_getoutfun(s);
            
            if size(FF, 2) > size(F, 2)
               F = repmat(F, 1, size(FF, 2));
            end;
            F(j:i, :) = FF;
            if i + 1 < length(strc.p)
               j = i + 1;
               p = strc.p(i + 1);
               if haspar2
                  p1 = strc.p1(i + 1);
               end;
            end;
         end;
      end;
      multassignin('base', par1, oldpar1);
      if haspar2
         multassignin('base', char(strc.parname2), oldpar2);
      end;
   end;
   g_Y = oldgY;
   g_t = oldgt;
catch err
   %   err=lasterror;
   multassignin('base', par1, oldpar1);
   if haspar2
      multassignin('base', char(strc.parname2), oldpar2);
   end;
   g_Y = oldgY; %#ok
   g_t = oldgt; %#ok
   rethrow(err);
end;
