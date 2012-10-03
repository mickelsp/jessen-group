% internal function of parameter analyser
function i_paranal(plotprev)
%parname, start, nend, nsteps, nstabilizing, ndays, lines, plotprev, outputtype, disturb)
%   i_paranal(answer.par{1}, answer.start(1), answer.nend(1), answer.steps(1), answer.stabil, ...
%      answer.writing,answer.lines,plotprev, answer.outputtype, answer.disturbance);

global t g_grind g_Y g_t g_paranal;
parname = g_grind.paranal.dlg.par{1};
start = g_grind.paranal.dlg.start(1);
nend = g_grind.paranal.dlg.nend(1);
nsteps = g_grind.paranal.dlg.steps(1);
nstabilizing = g_grind.paranal.dlg.stabil;
ndays = g_grind.paranal.dlg.writing;
lines = g_grind.paranal.dlg.lines;
outputtype = g_grind.paranal.dlg.outputtype;
silent=0;
if isfield(g_grind.paranal,'silent') &&( g_grind.paranal.silent==1);
   g_grind.paranal.silent=0;
   silent=1;
end;
%disturb=g_grind.paranal.dlg.disturbance;

if ~isfield(g_grind, 'paranal') || ~isfield(g_grind.paranal, 'defaultplots') || g_grind.paranal.defaultplots
   paranal('-defaults');
end;
parlim = [start, nend];
if start > nend
   parlim = [nend, start];
end;
for i = 1:length(g_grind.paranal.plots)
   if strcmp(g_grind.paranal.plots{i}.xaxis{1}, '<param1>')
      g_grind.paranal.plots{i}.xlim = parlim;
   end;
   if strcmp(g_grind.paranal.plots{i}.yaxis{1}, '<param1>')
      g_grind.paranal.plots{i}.ylim = parlim;
   end;
   if strcmp(g_grind.paranal.plots{i}.zaxis{1}, '<param1>')
      g_grind.paranal.plots{i}.zlim = parlim;
   end;
end;
if lines > 1
   lines = 0;
end;

outputlist = i_paranaldialog('outputlist');
outputtype = lower(outputlist{outputtype});

% disterror = 1E-4;
% if ~isempty(disturb)
%    if iscell(disturb)
%       disturb = disturb{1};
%    end;
%    if ~isempty(disturb) & isempty(strfind(disturb, ';'))
%       disturb = [disturb ';'];
%    end;
% end;
g_grind.solver.opt.OutputFcn = [];
% iY = i_varno(g_grind.paranal.plots{1}.yaxis{1});
% iZ = i_varno(g_grind.paranal.plots{1}.zaxis{1});
if silent
    h=[];
else
    parfig = i_figno('paranal');
[h, isnew] = i_makefig('paranal');
if  ~isempty(g_grind.paranal.plots{1}.zaxis{1}) && isnew
   set(gca, 'View', [8.5 10])
end;
i_odespeed(0, 0, 'init');
set(h, 'Name', 'Parameter analysis');
%oldhold = ishold;
hold on;
plotedit off;
nextpen;
if lines
   ppen = g_grind.pen.pen;
else
   ppen = '.';
end
end;
oldpar = evalin('base', char(parname));
oldtstep = g_grind.tstep;
% funx = [];
% funy = [];
% ps = [];
%try
hasperm = ~isempty(g_grind.permanent);
if plotprev ~= 1
   g_paranal.Y = [];
   g_paranal.t = [];
   g_paranal.p = [];
   %   pperm=[];
   if hasperm
      g_paranal.perm = [];
   end;
   g_paranal.parname = parname;
   %  wb=waitbar(0,'Calculating...');
   %try
   %   oldi = 1;
   %  newi = 0;
   t1 = t;
   N0 = i_initvar;
   OldN0 = N0;
   nstat = size(N0, 1);
   ppar = start;
   if isempty(g_grind.paranal.plots)
      xaxis1 = '';
      yaxis1 = '';
      zaxis1 = '';
   else
      xaxis1 = strrep(g_grind.paranal.plots{1}.xaxis{1}, '<param1>', parname);
      if ~isempty(g_grind.paranal.plots{1}.yaxis{1})
         yaxis1 = strrep(g_grind.paranal.plots{1}.yaxis{1}, '<param1>', parname);
      else
         yaxis1 = '';
      end;
      if ~isempty(g_grind.paranal.plots{1}.zaxis{1})
         zaxis1 = strrep(g_grind.paranal.plots{1}.zaxis{1}, '<param1>', parname);
      else
         zaxis1 = '';
      end;
   end;
   for i = 1:nsteps + 1
      %  waitbar(i/nsteps,wb);
      g_grind.solver.opt.OutputFcn = [];
      multassignin('base', char(parname), ppar);
      if ~g_grind.solver.isdiffer
         g_grind.tstep = 2;
      else
         g_grind.tstep = NaN;
      end;
      if nstabilizing > 0
         if i == 1 %extra initial stabilizing
            i_ru(g_grind.odefile,t1, 2 * nstabilizing, N0, 0);
         else
            i_ru(g_grind.odefile,t1, nstabilizing, N0, 0);
         end
         N0 = transpose(g_Y(size(g_Y, 1), :));
         drawnow;
      end
      g_grind.tstep = oldtstep;
      if ndays == 0
         g_t = i - 1;
         g_Y = transpose(N0);
      else
         if ~isempty(g_t)
            t1 = g_t(size(g_t, 1));
         else
            t1 = 0;
         end
         %          if ~isempty(disturb)
         %             ke;
         %             disterror = funeval(disturb);
         %             N0 = i_initvar;
         %          end;
         i_ru(g_grind.odefile, t1, ndays, N0, 0);
         drawnow;
         t1 = g_t(size(g_t, 1));
      end;
      p = repmat(ppar,size(g_t, 1), 1);
      N0 = transpose(g_Y(size(g_Y, 1), :));
      %little invasion to avoid hanging in a trivial equilibrium (tricky
      %if the model can have negative numbers
      epsil = 0.001;
      for j = 1:nstat
         if (N0(j) < epsil) && (N0(j) > 0)
            N0(j) = epsil;
         elseif (N0(j) > - epsil)  && (N0(j) < 0)
            N0(j) =  -epsil;
         elseif N0(j) == 0
            N0(j) = sign(OldN0(j)) * epsil;
         end
      end
      ppar = ppar + (nend - start) / nsteps;
      g_paranal.Y = [g_paranal.Y; g_Y];
      g_paranal.t = [g_paranal.t; g_t];
      g_paranal.p = [g_paranal.p; p];
      if hasperm
         pperm = defpermanent('-g', []);
         g_paranal.perm = [g_paranal.perm; pperm];
      end;
      %[xdata, ydata, zdata] = i_paranalfun(1,iY, iZ, outputtype, parname, g_t, g_Y, p, pperm, disterror);
      xdata = outfun(xaxis1);
      [ydata,ndx] = outfun(yaxis1, '-n',outputtype);
      if ~isempty(ndx) && (length(ndx)~=length(xdata))
         xdata = xdata(ndx);
      end;
      zdata = outfun(zaxis1, '-n',outputtype);
      if length(ydata) ~= length(zdata)
         zdata = outfun(zaxis1);
         zdata = zdata(ndx, :);
      end;
      if ~silent && g_grind.drawnow && (gcf == parfig)
         ud = get(parfig, 'userdata');
         if ~isempty(ud) && isfield(ud, 'stop') && ud.stop
            break;
         end;
         if i == 1
            ax = gca;
            set(ax, 'Xlim', g_grind.paranal.plots{1}.xlim);
            set(ax, 'Ylim', g_grind.paranal.plots{1}.ylim);
            set(ax, 'Zlim', g_grind.paranal.plots{1}.zlim);
            h = makeplot(1,ppen, xdata,ydata,zdata);
            set(h, 'Color', g_grind.pen.drawcolor);
            xlabel(mydisptext(g_grind.paranal.plots{1}.xaxis{1}, parname));
            ylabel(mydisptext(g_grind.paranal.plots{1}.yaxis{1}, parname));
            if ~isempty(g_grind.paranal.plots{1}.zaxis{1})
               zlabel(mydisptext(g_grind.paranal.plots{1}.zaxis{1}, parname));
            end;
         else
            if isempty(zdata)
               set(h, 'xdata', xdata(:,1), 'Ydata', ydata(:,1),  'Color', g_grind.pen.drawcolor);
            else
               set(h, 'xdata', xdata(:,1), 'Ydata', ydata(:,1), 'zdata', zdata(:,1), 'Color', g_grind.pen.drawcolor);
            end;
            drawnow;
         end;
      end;
   end;
   multassignin('base', char(parname), oldpar);
   if ishandle(h)
      delete(h);
   end;
end

%gather all data
%[ps, funx,funy, iY] = i_paranalfun(1,iY, iZ, outputtype, parname, g_paranal.t, g_paranal.Y, g_paranal.p, disterror);
%h = i_makefig('paranal');
if ~silent
i_odespeed(1, 1, 'done');
try
   for i = 1:length(g_grind.paranal.plots)
      nx = length(g_grind.paranal.plots{i}.xaxis);
      ny = length(g_grind.paranal.plots{i}.yaxis);
      nz = length(g_grind.paranal.plots{i}.zaxis);
      for k = 1:max([nx, ny, nz])
         ps = outfun(g_grind.paranal.plots{i}.xaxis{min(nx, k)}, '-p');
         [funx,ndx] = outfun(g_grind.paranal.plots{i}.yaxis{min(ny,k)},'-p', outputtype);
         if ~isempty(ndx) && (length(ndx)~=length(ps))
            ps = ps(ndx);
         end;
         funy = outfun(g_grind.paranal.plots{i}.zaxis{min(nz,k)},'-p', outputtype);
         if length(funy) ~= length(funx)
            funy = outfun(g_grind.paranal.plots{i}.zaxis{min(nz, k)}, '-p');
            funy = funy(ndx, :);
         end;
         %      iY = i_varno(g_grind.paranal.plots{i}.yaxis{1});
         %      iZ = i_varno(g_grind.paranal.plots{i}.zaxis{1});
         %      [ps, funx,funy, iY] = i_paranalfun(i, iY, iZ, outputtype, parname, g_paranal.t, g_paranal.Y, g_paranal.p, iif(hasperm,g_paranal.perm,[]), disterror);
         [h] = i_makefig('paranal',i - 1);
         
         set(gca, 'YLimMode', 'auto');
         set(gca, 'ZLimMode', 'auto');
         set(h, 'userdata', []);
         makeplot(i, ppen, ps, funx, funy);
      end
      set(h, 'Name', sprintf('Parameter analysis - %d',i));
      if ~isempty(g_grind.paranal.plots{i}.zaxis{1})
         set(gca, 'View', [8.5 10])
      end
      xlim(g_grind.paranal.plots{i}.xlim);
      xlabel(mydisptext(g_grind.paranal.plots{i}.xaxis{1}, parname));
      if ~strcmp(outputtype, 'unchanged')
         ylabel(mydisptext(sprintf('%s_{%s}', g_grind.paranal.plots{i}.yaxis{1}, outputtype),parname));
         zlabel(mydisptext(sprintf('%s_{%s}',g_grind.paranal.plots{i}.zaxis{1}, outputtype),parname));
      else
         ylabel(mydisptext(g_grind.paranal.plots{i}.yaxis{1}, parname));
         zlabel(mydisptext(g_grind.paranal.plots{i}.zaxis{1}, parname));
      end;
      set(findobj(gcf,'type','line'),'linewidth',g_grind.pen.linewidth);
      i_plotdefaults(gcf);
      hold on;
      plotedit on;
      %   if ~oldhold
      %      hold off;
      %   end;
   end;
   g_Y = [];
   g_t = [];
catch err
%    err=lasterror;
   %in case of error restore the old parameter (if functions
   %are used, the parameter has become a matrix)
   g_grind.tstep = oldtstep;
   multassignin('base', char(parname), oldpar);
   rethrow(err);
end;
end;
g_grind.tstep = oldtstep;
multassignin('base', char(parname), oldpar);


% function [xdata, ydata, zdata, iY] = i_paranalfun(iplot,iY, iZ, outputtype, parname, ts, Ys, ps, perms, disterror)
% global g_Y g_t g_grind;
% if isempty(iY) & ~isempty(strfind(outputtype, 'minima+maxima'))
%    m = max(g_Y);
%    iY=find(max(m) == m);
% end;
% if ~strcmp(outputtype, 'returntime')
%    [Ys, ps] = catfun(outputtype, Ys, ps, iY);
% end;
% if max(ps) > min(ps)
%    xdata = [];
%    ydata = [];
%    zdata = [];
%    %   istart = 1;
%    iend = 1;
%    while iend <= size(ps, 1)
%       p = ps(iend);
%       istart = iend;
%       while (iend <= size(ps, 1)) & (ps(iend, 1) == p)
%          iend = iend + 1;
%       end;
%       multassignin('base', char(parname), p);
%       g_Y = Ys(istart:iend - 1, :);
%       g_t = ts(istart:iend - 1, 1);
%       if ~isempty(perms)
%          g_permanent.t = g_t;
%          g_permanent.Y = perms;
%       end;
%       [ydat, zdat] = i_paranalsimpfun(iplot, iY, iZ, outputtype, parname, disterror);
%       %      if ~strcmp(g_grind.paranal.plots{iplot}.xaxis{1}, parname)&~strcmp(g_grind.paranal.plots{iplot}.xaxis{1}, '<param1>')
%       xdat =   getoutfun(g_grind.paranal.plots{iplot}.xaxis{1}, parname);
%       %      else
%       %         xdat =  ones(size(ydat, 1), 1) * p;
%       %     end;
%       ydata = [ydata; ydat];
%       zdata = [zdata; zdat];
%       xdata = [xdata; xdat];
%    end;
% else
%    g_Y = Ys;
%    g_t = ts;
%    [ydata, zdata] = i_paranalsimpfun(iplot, iY, iZ, outputtype, parname, disterror);
%    xdata = ps(1) * ones(size(ydata, 1), 1);
% end;

function s = mydisptext(s, parname)
if strcmp(s, '<param1>')
   s = parname;
end;
s = i_disptext(s);

% function [ydata, zdata] = i_paranalsimpfun(iplot, iY, iZ, outputtype, parname,   disterror)
% global g_Y g_grind;
% if strcmp(outputtype, 'returntime')
%    ydata = returntime(disterror, NaN, '-equil');
%    zdata = ydata;
% else
%    ydata = [];
%    if ~isempty(g_grind.paranal.plots{iplot}.yaxis{1})
%       for i = 1:length(g_grind.paranal.plots{iplot}.yaxis)
%          if (i > 1) | isempty(iY)
%             ydata = [ydata, getoutfun(g_grind.paranal.plots{iplot}.yaxis{i}, parname)];
%          else
%             ydata = g_Y(:, iY);
%          end;
%       end;
%    end;
%    zdata = [];
%    if ~isempty(g_grind.paranal.plots{iplot}.zaxis{1})
%       for i = 1:length(g_grind.paranal.plots{iplot}.zaxis)
%          if (i > 1) | isempty(iZ)
%             zdata = [zdata, getoutfun(g_grind.paranal.plots{iplot}.zaxis{i}, parname)];
%          else
%             zdata = g_Y(:, iZ);
%          end;
%       end;
%    end;
% end;

% function Adat = getoutfun(s, parname)
% Adat =   i_getoutfun(strrep(s, '<param1>', parname));
%
%
% function ERR = funeval(l_fun)
% ERR = 1E-4;
% eval(i_globalstr(who('global')));
% eval(l_fun);

function h = makeplot(iplot,ppen, ps, funx, funy)
global g_grind
maxn = max([size(ps, 2), size(funx, 2), size(funy, 2)]);
oldhold = ishold;
if ~isempty(funy)
   for i = 1:maxn
      h = plot3(ps(:,min(size(ps,2),i)), funx(:,min(size(funx,2),i)), funy(:,min(size(funy,2),i)), ppen, ...
         'MarkerSize', 3, 'EraseMode', 'none', 'Color', g_grind.pen.color2);
      set(gca, 'DrawMode', 'fast');
      hold on;
      nextpen;
      box on;
   end;
else
   for i = 1:maxn
      h = plot(ps(:,min(size(ps,2),i)), funx(:,min(size(funx,2),i)), ppen, ...
         'MarkerSize', 3, 'EraseMode', 'none', 'Color', g_grind.pen.color2);
      hold on;
      nextpen;
   end;
end;
if maxn > 1
   legend(g_grind.paranal.plots{iplot}.yaxis);
end;
if ~oldhold
   hold off;
end;
