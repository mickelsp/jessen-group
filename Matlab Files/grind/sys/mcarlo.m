%MCARLO   Monte Carlo uncertainty/sensitivity analysis
%   This function iterates a certain run, while changing the parameters and or initial conditions. 
%   Each parameter/initial condition can be selected separately and for each a distribution
%   for drawing each parameter independently can be selecte. You can select either a Uniform distribution 
%   (=rand), a normal distribution (=randn), a truncated norma distribution (negative number are discarded) 
%   or a log-normal distribution. You can also use this function to run a stochastic model several 
%   times without changing parameters.
%   The data (in g_mcarlo) are saved to a file, and after the run the data can be analysed. If you select
%   a uncertainty analysis time plots are generated with ranges (5%, 25% 50% 75% and 95% percentiles).
%   The sensitivity analysis performs a multivariate analysis of the parameter sensitivity (Klepper, 1989).
%   This analysis does not only find the overall effects of each parameter or initial condition, but it
%   also performs a cluster analysis, finding clusters of parameters with a similar effect on the results.
%   For this cluster analysis the statists toolbox is required (not for the sensitivity matrix).
%   
%   
%   Klepper, O. (1989) A model of carbon flows in relation to macrobenthic food supply in the Oosterschelde 
%   estuary (S.W. Netherlands), LUW Wageningen, The Netherlands. PhD thesis.
%   
%   
%   Usage:
%   MCARLO - all imput in dialog boxes.
%   MCARLO UNCERTAIN - uncertainty analysis, all other imput by user interface.
%   MCARLO g_mcarlo - All input is in this structure. Note that if g_mcarlo.i<g_mcarlo.n the previous run
%   is continued, else the previous data are used.
%   MCARLO SENS -load=mcarlo - sensitivity analysis while loading the data from mcarlo.mat
%
%   See also paranal, out, time

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function mcarlo(atype, opt)
i_parcheck;
evalin('base','global g_mcarlo;');
global g_mcarlo g_grind g_t g_Y t g_permanent g_paranal;
if (nargin == 1) && (isstruct(atype) || (ischar(atype)&&strcmp('g_mcarlo', atype)))
   if ~ischar(atype)
      g_mcarlo = atype;
   end;
   atype = '';
   opt = ' ';
elseif nargin  == 2
   if ischar(opt)&&(strncmpi(opt,'u',1)||strncmpi(opt,'s',1))
      h = opt;
      opt = atype;
      atype = h;
   end;
elseif nargin == 1
   if ischar(atype) && ~isempty(atype) && (atype(1)=='-')
      opt = atype;
      atype = '';
   else
      opt = '';
   end;
else
   opt = '';
   atype = '';
end;
atype = lower(atype);
if ischar(opt)
   if strncmpi(opt, '-l', 2)
      %load from file and continue if necessary
      f=strfind(opt, '=');
      if ~isempty(f)
         load(opt(f(1) + 1:end));
         if g_mcarlo.n == g_mcarlo.i
            disp('previous results loaded');
            if ~isempty(g_grind.permanent)
               time(1, '-s');
            end;
         end;
      end;
      if ~strcmpi(g_mcarlo.inifile, g_grind.inifile)
         warning('GRIND:mcarlo:inifile','Data generated with another inifile (%s) than currently opened',g_mcarlo.inifile);
      end;
   end;
end;
ndays = g_grind.ndays;
oldnsteps =  g_grind.tstep;
hasperm = ~isempty(g_grind.permanent);
if isempty(opt)
   i_mcarlodlg;
   uiwait;
   if ~isempty(g_mcarlo.paranalpar)
      Ans=inputdlg({'Number of iterations (paranal):', 'File for results'}, 'mcarlo',1 ,{'1000','mcarlo.mat'});
      if ~isempty(Ans)
         p = g_mcarlo.allpars{g_mcarlo.paranalpar};
         g_mcarlo.paranal.par = {p.name  ''};
         g_mcarlo.paranal.start = [p.minmax(1) 0];
         g_mcarlo.paranal.nend = [p.minmax(2) 10];
         g_mcarlo.paranal.steps = [p.steps 50];
         g_mcarlo.paranal.stabil = p.nstabilwrite(1);
         g_mcarlo.paranal.writing = p.nstabilwrite(2);
         g_mcarlo.paranal.hysteresis = strcmpi(p.distr, 'hysteresis');
         if g_mcarlo.paranal.writing(1)==0
            n = g_mcarlo.paranal.steps(1) * (g_mcarlo.paranal.writing(1) + 1);
         else
            n = (g_mcarlo.paranal.steps(1) + 1) * (g_mcarlo.paranal.writing(1) + 1)-1;
         end;
         Ans = [{'0'; num2str(ndays); num2str(n)}; Ans];
      end;
   else
      Ans=inputdlg({'Number of time steps to stabilize before each run','Number of time steps to simulate (each iteration):', 'Number of values of outcomes:',...
         'Number of iterations:', 'File for results'}, 'mcarlo',1 ,{'0',num2str(ndays),'10','1000','mcarlo.mat'});
   end;
   if isempty(Ans)
      disp('Cancel pressed')
      return;
   end;
   g_mcarlo.stabil=str2double(Ans{1});
   ndays = str2double(Ans{2});
   g_grind.tstep = str2double(Ans{3});
   if isempty(g_grind.tstep)
      error('GRIND:mcarlo:notfixed','None of the outcomes should be fixed');
   end;
   n = str2double(Ans{4});
   outfile = Ans{5};
   npar = 0;
   for j = 1:length(g_mcarlo.pars)
      if g_mcarlo.pars{j}.selected
         npar = npar + 1;
      end;
   end;
   g_mcarlo.n = n;
   g_mcarlo.i = 0;
   g_mcarlo.outfile = outfile;
   g_mcarlo.inifile = g_grind.inifile;
   g_mcarlo.p = zeros(n, npar);
   if g_grind.tstep<=1
      g_mcarlo.Ys = zeros(1, g_grind.statevars.dim, n);
   else
       g_mcarlo.Ys = zeros(g_grind.tstep + 1, g_grind.statevars.dim, n);
   end;
   if hasperm
      time(1, '-s')
      g_mcarlo.perm = zeros(size(g_mcarlo.Ys,1), g_grind.permanent{end}.to,n);
   else
      g_mcarlo.perm = [];
   end
end;
g_grind.solver.opt.OutputFcn = [];
try
   runned = 0;
   if ~isempty(g_mcarlo.paranalpar)
      g_grind.tstep =  max(1, g_mcarlo.paranal.writing);
   end;
   hwait = waitbar(0, 'Running mcarlo...');
   for i = g_mcarlo.i + 1:g_mcarlo.n
      runned = 1;
      drawnow;
      % setting the parameter values
      for j = 1:length(g_mcarlo.pars)
         switch g_mcarlo.pars{j}.distr
          case 'Uniform'
            ppar = g_mcarlo.pars{j}.min + rand(size(g_mcarlo.pars{j}.value)) * (g_mcarlo.pars{j}.max - g_mcarlo.pars{j}.min);
          case 'Normal'
            ppar =  drawrandn(size(g_mcarlo.pars{j}.value), g_mcarlo.pars{j}.value, g_mcarlo.pars{j}.sd);
          case 'TruncNormal'
            ppar = -1;
            while ppar < 0
               ppar =  drawrandn(size(g_mcarlo.pars{j}.value), g_mcarlo.pars{j}.value, g_mcarlo.pars{j}.sd);
            end;
          case 'LogNormal'
            ppar =  drawrandlogn(size(g_mcarlo.pars{j}.value), g_mcarlo.pars{j}.value, g_mcarlo.pars{j}.sd);
          otherwise
            error('GRIND:mcarlo:UnknownDistr','Distribution %s not supported', g_mcarlo.pars{j}.distr);
         end;
         multassignin('base', char(g_mcarlo.pars{j}.name), ppar);
         g_mcarlo.p(i, j) = ppar;
      end;
      %can also be state variables be changed
      N0 = i_initvar;
      if ~isempty(g_mcarlo.paranalpar)  
         i_keep(N0);
         %paranal (only one parameter allowed)
         if ~g_mcarlo.paranal.hysteresis  %paranal
            era;
            g_grind.paranal.silent=1;
            paranal(g_mcarlo.paranal)
            if i == 1
               g_mcarlo.t = g_paranal.t;
               g_mcarlo.paranal.p = g_paranal.p;
            end;
            g_mcarlo.Ys(:, :, i) = g_paranal.Y;
            if hasperm
               g_mcarlo.perm(:, :, i) =  g_paranal.perm;
            end;
         else  %hysteresis - back and forth run of paranal
            era;
            paranal(g_mcarlo.paranal);
            paranal('-1')
            if i == 1
               g_mcarlo.t = [g_paranal.prevt; g_paranal.t];
               g_mcarlo.paranal.p = [g_paranal.prevp; g_paranal.p];
            end;
            g_mcarlo.Ys(:, :, i) = [g_paranal.prevY; g_paranal.Y];
            if hasperm
               g_mcarlo.perm(:, :, i) = [g_paranal.prevperm; g_paranal.perm];
            end;
         end;
      else
         %running a normal model run
         if g_mcarlo.stabil>0 %optionally stabilize
             i_keep(N0);
             stabil(g_mcarlo.stabil,'-s');
             N0=i_initvar;
         end;
         if g_grind.tstep<=1
             stabil(ndays,'-s');
             g_Y=g_Y(end,:);
             g_t=g_t(end);
         else
             i_ru(g_grind.odefile,t, ndays, N0, 0);
         end;
         if i == 1
            g_mcarlo.t = g_t;
         end;
         g_mcarlo.Ys(:, :, i) = g_Y;
         if hasperm
            pperm = defpermanent('-g', []);
            g_mcarlo.perm(:, :, i) =  pperm;
         end;
      end;
      g_mcarlo.i = i;
      if mod(i, 10) == 0
         waitbar(i / g_mcarlo.n, hwait);
         save(g_mcarlo.outfile, 'g_mcarlo');
      end;
   end;
   if runned %save only if new data are generated
      save(g_mcarlo.outfile, 'g_mcarlo');
      fprintf('Saved data to %s\n', g_mcarlo.outfile);
   end;
   close(hwait);
   %set the parameters/statevars back to their original values
   g_grind.tstep = oldnsteps;
   for j = 1:length(g_mcarlo.pars)
      if g_mcarlo.pars{j}.selected
         multassignin('base', char(g_mcarlo.pars{j}.name), g_mcarlo.pars{j}.value);
      end;
   end;
catch err
   %  err=lasterror;
   if ishandle(hwait)
       close(hwait)
   end;
   g_grind.tstep = oldnsteps;
   save(g_mcarlo.outfile, 'g_mcarlo');
   for j = 1:length(g_mcarlo.pars)
      if g_mcarlo.pars{j}.selected
         multassignin('base', char(g_mcarlo.pars{j}.name), g_mcarlo.pars{j}.value);
      end;
   end;
   rethrow(err);
end;

%analyse the results either in an uncertainty analysis or sensitivity
%analysis
if isempty(atype)
   atype=questdlg('Select the analysis', 'Monte Carlo'...
      ,'Uncertainty analysis of the results','Sensitivity analysis of parameters','Sensitivity analysis of parameters');
end
if ~isempty(atype) && (atype(1) == 'U') %uncertainty analysis
   % the uncertainty analysis plots the timeplots but with the 5% 25% 50% 75% 95% percentiles
   for No = 1:size(g_grind.timevars, 2)
      if ~isempty(g_grind.timevars{No})
         i_makefig('time', No - 1);
         hold on;
         plotedit('off');
         apen = g_grind.pen;
         apen.i = 1; %  do not cycle colors between runs
         apen.nextpen = 1;
         
         if strcmp(g_grind.outt{No}, 't')
            tt = g_mcarlo.t; %efficiency
         else
            %repeated calls to outfun are not efficient, maybe needs to
            %be reimplemented
            tt = mean(reshape(outfun(g_grind.outt{No}, '-m'), length(g_mcarlo.t), g_mcarlo.n), 2);
         end;
         %use the variables defined in 'out'
         for i = 1:length(g_grind.timevars{No})
            s = g_grind.timevars{No}{i};
            if ~strncmp(s, 'Observed ', 9)
               xx = reshape(outfun(s, '-m'), length(g_mcarlo.t), g_mcarlo.n);
               xx =  i_makepercentiles(xx, [0.05, 0.25, 0.5, 0.75, 0.95]);
               plot(tt, xx(:,[1 5]), apen.pen, 'Color', apen.color,'Linewidth',0.5, 'LineStyle',':');
               plot(tt, xx(:,[2 4]), apen.pen, 'Color', apen.color,'Linewidth',1);
               plot(tt, xx(:,3), apen.pen, 'Color', apen.color,'Linewidth',2);
               apen = nextpen(apen);
            end;
         end;
         xlabel( g_grind.outt{No});
         s = i_disptext(g_grind.timevars{No}{1});
         ylabel(s);
      end;
   end;
elseif ~isempty(atype)
   %sensitivity analysis
   %disp('sensitivity analysis not yet implemented')
   % here code to create (1) a sensitivity matrix
   % 1 matrix per kolom (voor alle SI verschillend)
   
   % - define parameters
   Nsteps = length(g_mcarlo.t);
   Niter = g_mcarlo.i;
   Npars = size(g_mcarlo.p, 2);
   
   %2 - get data
   %dit is alleen van belang voor vector variabelen,
   %bvzou ook relevante willen selecteren, bijv met out
   
   if isfield(g_mcarlo, 'out')
      s = strtrim(sprintf('%s\n', g_mcarlo.out{:}));
   elseif g_grind.statevars.vector
      s = strtrim(sprintf('%s\n', g_grind.statevars.vectnames{:}));
   else
      s = strtrim(sprintf('%s\n', g_grind.statevars.names{:}));
   end
   answer=inputdlg({sprintf('Enter variables or functions for sensitivity analysis\n(one per line)')},'mcarlo',[10,60],{s});
   if isempty(answer)
      disp('Cancel pressed');
      return;
   end;
   hwait = waitbar(0, 'Analysing results...');
   drawnow;
   vectout = str2cell(answer{1});
   vectout = vectout(~strcmp(vectout, ''));
   %    s = strtrim(sprintf('%s ',s{:}));
   %    s=strrep(s,'  ',' '); %remove double spaces
   %    s=strrep(s,'  ',' ');
   %    f = [0 strfind(s, ' ') length(s) + 1];
   %    vectout = cell(length(f) - 1, 1);
   %    for i = 2:length(f)
   %       vectout{i - 1} = s(f(i - 1) + 1:f(i) - 1);
   %    end
   Nstate = 0;
   g_Y = g_mcarlo.Ys(:, :, 1); %#ok
   g_t = g_mcarlo.t;
   if hasperm
      g_permanent.Y = g_mcarlo.perm(:, :, 1);
      g_permanent.t = g_t;
   end;
   for i = 1:length(vectout);
      ss = outfun(vectout{i}); % the problem is that this may have a variable dimension
      Nstate = Nstate + size(ss, 2);
   end;
   g_Y = [];
   g_t = [];
   mc_output = NaN + zeros(Niter, Nsteps * Nstate);
   g_mcarlo.out = cell(Nstate, 1);
   k = 1;
   l = 1;
   for i = 1:length(vectout)
      ss = outfun(vectout{i}, '-m');
      waitbar(i / (length(vectout) + 1), hwait);
      nvar = size(ss, 2);
      for j = 1:nvar
         if nvar == 1
            g_mcarlo.out{l} = vectout{i};
         else
            g_mcarlo.out{l} = sprintf('%s(%d)', vectout{i}, j);
         end;
         mc_output(:, k:k + Nsteps - 1) = reshape(ss(:,j), Nsteps, Niter)';
         %I didn't manage to do this in one step
         k = k + Nsteps;
         l = l + 1;
      end;
   end;
   plabels = cell(Npars, 1);
   for j = 1:Npars
      plabels{j} = g_mcarlo.pars{j}.name;
   end
   vlabels = cell(Nsteps * Nstate, 2);
   for i = 1:Nstate
      for j = 1:Nsteps;
         vlabels{(i - 1) * Nsteps + j, 1} = g_mcarlo.out{i};
         vlabels{(i - 1) * Nsteps + j,2} = sprintf('%g', g_mcarlo.t(j));
      end;
   end;
   close(hwait);
   %do the multivariate sensitivity analysis
   if ~isempty(g_mcarlo.outfile)
      save(g_mcarlo.outfile, 'g_mcarlo');
   end;
   i_mcsensanal(mc_output, g_mcarlo.p, vlabels, plabels)
end

function res = drawrandlogn(siz, mean, sd)
% draw logNormal distribution with expected mean of mean and expected variance of sd^2
s = ln(1 + (sd / mean)^2)^0.5; %this is the sd needed to get a standard deviation of sd
mu = ln(mean) - 0.5 * s; %mu is needed to get an expected value of mean
res = exp(randn(siz) .* s + mu);
function res = drawrandn(siz, mean, sd)
res = randn(siz) .* sd + mean;



