%OPTIMPARS   Fit data by optimizing parameters
%   Optimize parameters with the local simplex method or a global search (Shuffled Complex 
%   Evolution (SCE-UA)). First, a data set must be entered using setdata (a dialog box 
%   for entering data appears if this has not yet been done). Criterion is the sum 
%   of squares of normalized state variables.
%
%   Usage:
%   OPTIMPARS - Enter parameters and options in dialog boxes (Simplex method).
%   OPTIMPARS -g = Enter parameters and options in dialog boxes (global method).
%   OPTIMPARS P1 P2 P3 .. - Optimize the list of parameters P1, P2, P3 using simplex optimization.
%   OPTIMPARS -local P1 P2 P3 .. - use simplex method (default).
%   OPTIMPARS -global P1 [L1 U1] P2 [L3 U3] P3 [L3 U3] .. - using global optimization and ranges for 
%   each parameter [L1=Lower P1, U1=Upper P1].
%   OPTIMPARS -reset - Reset the parameters as before the previous run.
%
%
%   See also setdata, fglobalmin, fminsearch

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function optimpars(varargin)
global g_data g_grind g_Y;
i_parcheck;
if nargin == 1  && isstruct(varargin{1})
   g_data = varargin{1};
   if strncmp(g_data.optmethod, 'Nelder', 6)
      optmethod = 'fminsearch';
   else
      optmethod = 'fglobalmin';
   end;
else
   parmins = [];
   parmaxs = [];
   optmethod = 'fminsearch';
   if (nargin == 1) && strncmpi(varargin{1}, '-r', 2)
      %reset previous parameter settings
      if ~isempty(g_data) && isfield(g_data, 'X0')
         assignpars(g_data.X0);
         g_Y = [];
         disp('Parameters resetted to previous values');
         return;
      else
         error('GRIND:optimpars:CannotReset','Cannot reset parameters, no previous values available');
      end;
   end;
   if (nargin >= 1) && strncmpi(varargin{1}, '-g', 2)
      optmethod = 'fglobalmin';
      if nargin == 1
         setdata;
         if isfield(g_data, 'pars') && ~isempty(g_data.pars)
            if ~isfield(g_data, 'parmins')  || isempty(g_data.parmins)
               g_data.parmins = zeros(length(g_data.pars), 1);
               g_data.parmaxs = g_data.parmins + 1;
               if isfield(g_data, 'X0')
                  g_data.parmins = g_data.X0 * 0.5;
                  g_data.parmaxs = g_data.X0 * 1.5;
               else
                  g_data.parmins = zeros(length(g_data.pars), 1);
                  g_data.parmaxs = g_data.parmins + 1;
               end;
               
            end;
            inputs = cell(1,length(g_data.pars)*2);
            for i = 1:length(g_data.pars)
               inputs{(i-1)*2+1} = g_data.pars{i};
               inputs{(i-1)*2+2} = [g_data.parmins(i) g_data.parmaxs(i)];
             %   inputs = {inputs{:} g_data.pars{i} [g_data.parmins(i) g_data.parmaxs(i)]};
           end;
            s1=strrep(i_cell2str(inputs),'''','');
            s1 = s1(2:length(s1) - 1);
         else
            s1 = '';
         end;
         answer=inputdlg('Specify parameters and ranges in a list: P1 [L1 U1] P2 [L2 U2] ...', 'optimpars', 1, {s1});
         if isempty(answer{1})
            errordlg('You need to specify the parameters in a list e.g. "optimpars r [0.1 0.5] K  [1 10]"');
            error('GRIND:optimpars:NoParList','You need to specify the parameters in a list e.g. "optimpars r [0.1 0.5] K  [1 10]"');
         end;
         s=['''' strtrim(answer{1})];
         s=strrep(s,'  ',' ');
         s= strrep(s, ' [', ''' [');
         s= strrep(s, '] ', '] ''');
         inputs=eval(['{' s '}']);
      else
         inputs = varargin(2:end);
      end;
      varargin = cell(1,round(length(inputs)/2));
      parmins=zeros(1,round(length(inputs)/2));
      parmaxs=parmins;
      i = 1;
      par = 1;
      while i < length(inputs);
         varargin{par} = inputs{i};
         i = i + 1;
         minmax = i_checkstr(inputs{i});
         i = i + 1;
         parmins(par) = min(minmax);
         parmaxs(par) = max(minmax);
         par = par + 1;
      end;
   end;
   if (nargin >= 1) && strncmpi(varargin{1}, '-l', 2)
      varargin = varargin(2:end);
   end;
   if ~isfield(g_data, 'options') || isempty(g_data.options)
      opts = optimset('fminsearch');
      g_data.options.Display = opts.Display;
      g_data.options.TolX = 1E-6;
      g_data.options.TolFun = 1E-6;
      g_data.options.MaxFunEvals1 = opts.MaxFunEvals;
      g_data.options.MaxIter1 = opts.MaxIter;
      g_data.options.PosOnly = 0;
      g_data.options.ResetEachStep  = 0;
   end;
   if ~isfield(g_data, 'options2') || isempty(g_data.options2)
      opts =  optimset('fglobalmin');
      g_data.options2 = opts;
      g_data.options2.PosOnly = 0;
      g_data.options2.ResetEachStep = 0;
   end;
   
   if isempty(varargin)
      setdata;
  %    prompt={'Table with data:','List of variables in columns of the above table:'};
      if isfield(g_data, 'pars') && ~isempty(g_data.pars)
         s1=strrep(i_cell2str(g_data.pars),''' ''','  ');
         s1 = s1(3:length(s1) - 2);
      else
         s1 = '';
      end;
      answer=inputdlg('Specify parameters in a list delimited with spaces', 'optimpars', 1, {s1});
      if isempty(answer{1})
         errordlg('You need to specify the parameters in a list e.g. "optimpars r K"');
         error('GRIND:optimpars:NoParList','You need to specify the parameters in a list e.g. "optimpars r [0.1 0.5] K  [1 10]"');
      end;
      s=['''' strtrim(answer{1}) ''''];
      s=strrep(s,'  ',' ');
      s= strrep(s, ' ', ''' ''');
      varargin=i_checkstr(['{' s '}']);
   elseif isempty(g_data) || ~isfield(g_data, 'obs') || isempty(g_data.obs)
      setdata;
   end;
   if iscell(varargin{1})
      varargin = varargin{1};
   end;
   
   if strcmp(optmethod, 'fminsearch')
      prompt={'Termination tolerance on X (TolX)','Termination tolerance on the function value (TolFun)',...
            'Maximum number of function evaluations allowed (MaxFunEvals)','Maximum number of iterations allowed (MaxIter)',...
            'Penalty function (SS, normalized SS or user function)',...
         'Only positive pars? (1/0)','Reset each step? (1/0)'};
      answer = {num2str(g_data.options.TolX), num2str(g_data.options.TolFun), g_data.options.MaxFunEvals1, ...
         g_data.options.MaxIter1,'normalized SS', num2str(g_data.options.PosOnly),num2str(g_data.options.ResetEachStep)};
      answer = inputdlg(prompt, 'Change options', 1, answer);
      if ~isempty(answer)
         g_data.options.TolX = str2double(answer{1});
         g_data.options.TolFun = str2double(answer{2});
         g_data.options.MaxFunEvals = str2double(answer{3});
         g_data.options.MaxIter = str2double(answer{4});
         g_data.options.penaltyfun = answer{5};
         g_data.options.PosOnly = str2double(answer{6});
         g_data.options.ResetEachStep = str2double(answer{7});
      end;
   else
      prompt={'Maximum number of function evaluations','Convergence if normalized parameter space <',...
         'Maximum number of evolution loops before convergency', ...
         'The percentage change allowed in MaxIter loops before convergency', ...
         'Number of complexes (sub-populations)', ...
         'Penalty function (SS, normalized SS or user function)',...
         'Reset each step? (1/0)'};
      answer = {sstr(g_data.options2.MaxFunEvals), sstr(g_data.options2.TolX), sstr(g_data.options2.MaxIter), ...
         sstr(g_data.options2.TolFun), sstr(g_data.options2.MaxPCGIter), 'normalized SS',  ...
         sstr(g_data.options2.ResetEachStep)};
      answer = inputdlg(prompt, 'Change options', 1, answer);
      if ~isempty(answer)
         g_data.options2.MaxFunEvals = str2double(answer{1});
         g_data.options2.TolX = str2double(answer{2});
         g_data.options2.MaxIter = str2double(answer{3});
         g_data.options2.TolFun = str2double(answer{4});
         g_data.options2.MaxPCGIter = str2double(answer{5});
         g_data.options2.PosOnly = 0;
         g_data.options.penaltyfun = answer{6};
        g_data.options2.ResetEachStep = str2double(answer{7});
      end;
   end;
   g_data.pars = {};
   switch optmethod
    case 'fminsearch'
      g_data.optmethod = 'Nelder-Mead simplex (direct search) method';
    case 'fglobalmin'
      g_data.optmethod = 'Shuffled Complex Evolution (SCE-UA) Duan, et al. 1992';
      g_data.parmins = parmins;
      g_data.parmaxs = parmaxs;
   end;
   n = 0;
   for i = 1:size(varargin, 2)
      try
         kk = evalin('base', char(varargin{i}));
         if numel(kk) == 1
            g_data.pars = [g_data.pars varargin(i)];
            n = n + 1;
         else
            for k1 = 1 :size(kk, 1)
               for k2 = 1 :size(kk, 2)
                  g_data.pars = [g_data.pars {sprintf('%s(%d,%d)',varargin{i},k1,k2)}];
                  n = n + 1;
               end;
            end;
         end;
      catch err
      %   err=lasterror;
         errordlg(['Error reading parameter ' char(varargin{i})]);
         rethrow(err);
      end;
   end;
end;
g_grind.solver.reset_at_data=g_data.options.ResetEachStep;
g_data.X = [];
g_data.fval = [];
g_data.iter = 0;
g_data.stopped = 0;
g_data.minobs = min(g_data.obs);
g_data.maxobs = max(g_data.obs);

par0 = zeros(1, size(g_data.pars, 2));
for i = 1:size(g_data.pars, 2)
   par0(i) = evalin('base', char(g_data.pars{i}));
end;
oldtstep = g_grind.tstep;
g_data.X0 = par0;
g_grind.tstep = NaN;
try
   g_grind.solver.opt.OutputFcn = [];
   numberofvariables = length(g_data.pars); %#ok needed
 %  numberOfVariables = length(g_data.pars); %version <7.1
   g_data.options.MaxFunEvals  = eval(lower(g_data.options.MaxFunEvals1));
   g_data.options.MaxIter  = eval(lower(g_data.options.MaxIter1));
   H0 = i_optimpardlg;
   if strcmp(optmethod, 'fglobalmin')
      [X, fval, found] = fglobalmin('i_optimpars', g_data.parmins,g_data.parmaxs, g_data.options2);
   else
      [X, fval, found] = fminsearch('i_optimpars', par0, g_data.options);
   end;
   assignpars(X);
   if ishandle(H0)
      close(H0);
   end;
   g_grind.solver.opt.OutputFcn = str2func('i_odespeed');
   g_grind.tstep = oldtstep;
catch err
 %  err=lasterror;
   %if 0
   g_grind.solver.opt.OutputFcn = str2func('i_odespeed');
   g_grind.tstep = oldtstep;
   fval =    g_data.fval;
   if ishandle(H0)
      close(H0);
   end;
   found=g_data.stopped == 1;
   if found
      assignpars(g_data.X);
      disp('Optimization stopped by user, using best solution so far');
   elseif g_data.stopped == 2
      assignpars(g_data.X0);
      g_Y = []; %#ok
      rethrow(err);
   else
      rethrow(err);
   end;
end;
if found
   fprintf('Method            : %s',g_data.optmethod );
   fprintf('Penalty function  : %s\n',g_data.options.penaltyfun);
   adjr2(length(g_data.pars));
   disp(['Sum of squares    : ' num2str(fval)]);
   disp(['No. of iterations : ' num2str(g_data.iter)]);
   disp(' ');
end;
i = i_figno('time');
if ishandle(i)
   set(i, 'visible', 'off');
end;

function multassignin(ws, name, V, fbrack)
%multassignin supports assignments to elements of arrays (e.g. name='A(1,2)')
%fbrack should be strfind(name, '(') (added for efficiency)
if isempty(fbrack)
   assignin(ws, name, V);
else
   temp = evalin(ws, name(1, fbrack(1) - 1));
   s = ['temp' name(fbrack(1):length(name)) ' = ' num2str(V) ';' ];
   eval(s);
   assignin(ws, name(1, fbrack(1) - 1), temp);
end;
function assignpars(x)
global g_data;
g_data.pred = [];
for i = 1:size(g_data.pars, 2)
   fbrack = strfind(g_data.pars{i}, '(');
   multassignin('base', g_data.pars{i}, x(i),fbrack);
end;

function r = sstr(n)
if ischar(n)
   r = n;
else
   r = num2str(n);
end;

