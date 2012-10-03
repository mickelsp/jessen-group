%PARANAL   1D parameter analyser
%   Change one parameter step-by-step and show the attractor by simulation. 
%   A 2D or 3D figure is created with the results (on X axis the parameter,
%   on Y axis the parameter/function of the x-axis of phase plane and on 
%   the Z-axis the parameter/function of the y-axis of the phase plane (if any).
%
%   Usage:
%   PARANAL - user is prompted for information
%   PARANAL 1 - replot the previous analysis
%   PARANAL -1 - reanalyse the model in the opposite direction
%   PARANAL -out (-o) - change the default output in a dialog box.
%   PARANAL -out plotno [ 'funy1 funy2' funz1] [minx maxx] [miny maxy] [minz maxz] 
%   - sets the output in a command line: plotno  = number of plot,  is first parameter
%   funy1 = 1st variable or function for y axis [minx maxx] range for xaxis, etc.
%   PARANAL -default (-d) - reset the default output.
%   PARANAL -list (-l) - list the outputs for paranal.
%   PARANAL -save=filename (-s) - save the results of the last or current paranal to a mat 
%   file with name "filename.mat".
%   PARANAL -load=filename (-lo) - load the results of a previous paranal session 
%   (opens if needed the same inifile).
%
%   See also ax, conteq, paranal2d

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function paranal(varargin)
global g_paranal g_grind;
if (nargin==0)|| ~ischar(varargin{1}) || ~strncmpi(varargin{1}, '-lo', 3)
   i_parcheck;
end;
outfile = '';
if isfield(g_grind,'pars')&&isempty(g_grind.pars)
   error('GRIND:paranal:NoPars','No parameters to analyse');
end;
answer = [];
plotprev = 0;
if nargin > 0
   if isfield(g_grind,'paranal')&&isfield(g_grind.paranal,'dlg');
      answer=g_grind.paranal.dlg;
   end;
   for ii = 1:nargin
      arg = varargin{ii};
      if isstruct(arg)
         answer = i_paranaldialog('initstruct', arg);
      elseif ischar(arg) && strncmpi(arg, '-s', 2)
         outfile = 'paranal.mat';
         f=strfind(arg, '=');
         if ~isempty(f)
            outfile = arg(f(1) + 1:end);
         else
            outfile= fullfile(grindpath,outfile);
         end;
         if ~isempty(g_paranal) && (nargin==1)
            saveparanal(outfile);
            return;
         end;
      elseif ischar(arg) && strncmpi(arg, '-lo', 3)
         infile = 'paranal.mat';
         f=strfind(arg, '=');
         if ~isempty(f)
            infile = arg(f(1) + 1:end);
         else
            infile= fullfile(grindpath,infile);
         end;
         f = strfind(infile, '.');
         if isempty(f)
            infile = [infile '.mat']; %#ok (faster than sprintf)
         end;
         if ~exist(infile, 'file')
            [infile,path]=uigetfile('*.mat','File for paranal results');
            if isempty(infile)
               disp('File not found');
               return;
            end;
            cd(path);
         end
         loadparanal(infile);
          if nargin == 1
            paranal(1);
            return;
         end;
      elseif ischar(arg) && strncmpi(arg, '-d', 2)
         g_grind.paranal.plots = {};
         g_grind.paranal.plots{1}.xaxis = {'<param1>'};
         g_grind.paranal.plots{1}.xlim = [0 10];
         g_grind.paranal.plots{1}.yaxis = {g_grind.xaxis.var};
         g_grind.paranal.plots{1}.ylim = g_grind.xaxis.lim;
         g_grind.paranal.plots{1}.zaxis = {g_grind.yaxis.var};
         g_grind.paranal.plots{1}.zlim = g_grind.yaxis.lim;
         g_grind.paranal.defaultplots = 1;
         g_grind.paranal.currno = 1;
         return;
      elseif ischar(arg) && strncmpi(arg, '-l', 2)
         if ~isfield(g_grind, 'paranal')
            paranal('-default');
         end;
         for i = 1:length(g_grind.paranal.plots)
            fprintf('paranal -out %d [''%s'' ''%s'' ''%s''] [%g %g] [%g %g] [%g %g]\n',i, strtrim(sprintf('%s ',g_grind.paranal.plots{i}.xaxis{:})),...
               strtrim(sprintf('%s ',g_grind.paranal.plots{i}.yaxis{:})),strtrim(sprintf('%s ',g_grind.paranal.plots{i}.zaxis{:})),...
               g_grind.paranal.plots{i}.xlim, g_grind.paranal.plots{i}.ylim, g_grind.paranal.plots{i}.zlim);
         end;
         return;
      elseif ischar(arg) && strncmpi(arg, '-o', 2)
         if nargin == ii
            i_paranaloutdlg('init');
         else
            p.xaxis = {''};
            p.yaxis = {''};
            p.zaxis = {''};
            ip = 1;
            j = 0;
            for i = ii + 1:nargin
               v = str2num(varargin{i}); %#ok
               if isempty(v)
                  %analyse axes
                  s = varargin{i};
                  if s(1) == '['
                     s = s(2:end - 1);
                  end;
                  f=strfind(s,'''');
                  if isempty(f)
                     f = strfind(s, ' ');
                     ff = repmat(f,2,1); %this is a trick to get double values
                     f = [0; ff(:); length(s) + 1];
                  end;
                  if length(f) > 1
                     p.xaxis = mystr2cell(s(f(1) + 1:f(2) - 1));
                  end
                  if length(f) > 3
                     p.yaxis = mystr2cell(s(f(3) + 1:f(4) - 1));
                  end
                  if length(f) > 5
                     if f(5) < f(6) - 1
                        p.zaxis = mystr2cell(s(f(5) + 1:f(6) - 1));
                     end;
                  end
               elseif length(v) == 1
                  ip = v;
               elseif length(v) == 2
                  switch j
                   case 0
                     p.xlim = v;
                   case 1
                     p.ylim = v;
                   case 2
                     p.zlim = v;
                  end;
                  j = j + 1;
               end
            end;
            g_grind.paranal.plots{ip} = p;
         end
         return;
      else
         plotprev = i_checkstr(arg);
         if plotprev == 2
            paranal(1);
            plotprev = 1;
         end;
         if (plotprev ==1) && isempty(answer)
            p.par = {g_paranal.parname};
            p.start = min(g_paranal.p);
            p.nend = max(g_paranal.p);
            answer = i_paranaldialog('initstruct', p);
         end;
      end;
   end;
end;
if isempty(answer)
   if isfield(g_grind, 'paranal')&&isfield(g_grind.paranal, 'dlg')
      answer = g_grind.paranal.dlg;
   end
   if plotprev == 0
      a1 = i_paranaldialog(answer);
      if isempty(a1)
         return
      elseif ~isempty(answer) && isfield(g_paranal, 'Y') && ~isempty(g_paranal.Y) && ~isempty(g_paranal.p) && ...
            prod(double(a1.start==answer.start)) && strcmp(a1.par{1}, answer.par{1}) ...
            && strcmp(a1.par{2}, answer.par{2}) && prod(double(a1.steps==answer.steps)) ...
            && (a1.stabil==answer.stabil) && (a1.writing==answer.writing) && ...
            prod(double(a1.nend==answer.nend)) && ...
            strcmp(questdlg('Do you want to use data of the previous paranal run?','paranal',...
            'Yes','No','No'),'Yes')
         plotprev = 1;
      end;
      answer = a1;
   end;
end;
if ~isempty(answer)
   if ~isfield(g_grind, 'paranal')
      g_grind.paranal.defaultplots = 1;
   end;
   answer.par{1} = checkmat(answer.par{1});
   answer.par{2} = checkmat(answer.par{2});
   if plotprev == -1
      s = answer.start(1);
      answer.start(1) = answer.nend(1);
      answer.nend(1) = s;
   end;
   g_grind.paranal.dlg = answer;
   if ~isempty(answer.par{2})
      paranal2d(answer);
      clear answer;
      return;
   end;
   if answer.lines == 2
      i_warningdlg('GRIND:paranal:contour1D','Contour plot not supported in 1D paranal, making scatterplot instead');
   end;
   if answer.lines == 3
      i_warningdlg('GRIND:paranal:surface1D','Surface plot not supported in 1D paranal, making scatterplot instead');
   end;
   if plotprev == -1
      g_paranal.Yprev = g_paranal.Y;
      g_paranal.pprev = g_paranal.p;
      g_paranal.tprev = g_paranal.t;
   end;
   %   function i_paranal(parname, start, nend, nsteps, nstabilizing, ndays, lines, plotprev, outputtype, disturb);
   i_paranal(plotprev);
   if (plotprev == 1) && isfield(g_paranal, 'Yprev') &&~isempty(g_paranal.Yprev)
      Y = g_paranal.Y;
      g_paranal.Y = g_paranal.Yprev;
      g_paranal.Yprev = Y;
      p = g_paranal.p;
      g_paranal.p = g_paranal.pprev;
      g_paranal.pprev = p;
      t1 = g_paranal.t;
      g_paranal.t = g_paranal.tprev;
      g_paranal.tprev = t1;
   end;
   if ~isempty(outfile)
      saveparanal(outfile);
   end;
   clear answer;
end;
function saveparanal(outfile)
global g_paranal g_grind;  %#ok
path = pwd;                %#ok
inifile = g_grind.inifile; %#ok
dlg = g_grind.paranal.dlg;   %#ok
save(outfile,'g_paranal','path','inifile','dlg');
fprintf('Saved paranal results to "%s"\n', outfile)
function loadparanal(infile)
global g_paranal g_grind; %#ok
path = '';
inifile = '';
dlg=[];
load(infile);
if isempty(path)
   errordlg('Not a valid file for paranal');
   error('GRIND:paranal:NoFile','Paranal: not a valid file');
end;
if ~strcmp(path, pwd);
   cd(path);
end;
if isempty(g_grind)||(~isfield('inifile',g_grind)&&~strcmp(inifile, g_grind.inifile));
   use(inifile);
   load(infile);
end;
g_grind.paranal.dlg=dlg;
function p = checkmat(par)
if ~isempty(par) && isempty(strfind(par, '('))
   s=evalin('base',sprintf('size(%s)',par));
   answer{1} = par;
   if (max(s) > 1) && (min(s) == 1)
      answer=inputdlg('Enter which element of parameter you want to use',sprintf('Parameter is %dx1 vector',max(s)),1,{[par '(1)']});
   end;
   if (max(s) > 1) && (min(s) > 1)
      answer=inputdlg('Enter which element of parameter you want to use',sprintf('Parameter is %dx%d matrix',s),1,{[par '(1,1)']});
   end;
   p = answer{1};
else
   p = par;
end;
function A = mystr2cell(s)
s = outf('changeshortcut', s);
A=str2cell(strrep(s,' ',sprintf('\n')));
