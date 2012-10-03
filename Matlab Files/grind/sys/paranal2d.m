%PARANAL2D   2D parameter analyser
%   Change two parameters step-by-step and show the attractor by simulation. 
%   Contour plots are created with mean value of the attractor.
%
%   Usage:
%   PARANAL2 - user is prompted for information
%   
%   See also paranal

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function [Res, X, Y] = paranal2d(answer)
evalin('base','global g_paranal2d;');
global g_paranal g_grind g_paranal2d;
if isfield(g_grind,'pars')&&isempty(g_grind.pars)
   error('GRIND:paranal2d:NoPars','No parameters to analyse');
end;
i_parcheck;
plotprev = 0;
if (nargin == 1) && ischar(answer) && (answer == '1')
   parname1 = g_grind.paranal.dlg.par{1};
   parname2 = g_grind.paranal.dlg.par{2};
   lines = g_grind.paranal.dlg.lines;
   outputtype = g_grind.paranal.dlg.outputtype;
else
   if nargin == 1
      if ischar(answer) && (strcmp(answer,'-1'))
         g_paranal2d.prevp = g_paranal2d.p;
         g_paranal2d.prevp1 = g_paranal2d.p1;
         g_paranal2d.prevY = g_paranal2d.Y;
         answer = g_grind.paranal.dlg;
         h = answer.start(1);
         answer.start(1) = answer.nend(1);
         answer.nend(1) = h;
         g_grind.paranal.dlg = answer;
      end;
      %    answer = i_paranaldialog('initstruct', answer);
   else
      if isfield(g_grind, 'paranal')
         answer = g_grind.paranal.dlg;
      else
         answer.paranal2d = 1;
      end
      if ~isfield(answer, 'paranal')
         answer.paranal2d = 1;
      end;
      a1 = i_paranaldialog(answer);
      if isempty(a1)
         return
      elseif isfield(g_grind, 'paranal2d')&&~isempty(answer)&&isfield(g_paranal, 'Y') && ~isempty(g_paranal.Y) && ~isempty(g_paranal.p) && ...
            prod(a1.start == answer.start) && strcmp(a1.par{1}, answer.par{1}) ...
            && strcmp(a1.par{2}, answer.par{2}) && prod(a1.steps == answer.steps) ...
            && (a1.stabil == answer.stabil) && (a1.writing == answer.writing) && ...
            prod(a1.nend == answer.nend) && ...
            strcmp(questdlg('Do you want to use data of the previous paranal2d run?','paranal2d',...
            'Yes','No','No'),'Yes')
         plotprev = 1;
      end;
      answer = a1;
      g_grind.paranal.dlg = answer;
      g_grind.paranal.defaultplots=1;
   end;
   ipar1 = 1;
   ipar2 = 2;
   if isempty(answer.par{ipar2}) || strcmpi(answer.par{ipar2}, 'none')
      ipar2 = 1;
      ipar1 = 2;
   end;
   parname1 = answer.par{ipar1};
   nsteps1 = answer.steps(ipar1);
   writing1 = answer.writing(ipar1);
   parname2 = answer.par{ipar2};
   outputtype = answer.outputtype;
   if iscell(parname2)
      parname2 = parname2{1};
   end;
   start2 = answer.start(ipar2);
   nend2 =  answer.nend(ipar2);
   nsteps2 =  answer.steps(ipar2);
   nopar2 = isempty(parname2) | strcmpi(parname2, 'none');
   lines = answer.lines;
   g_paranal2d.parname1=parname1;
   g_paranal2d.parname2=parname2;
   %try
   if plotprev ~= 1
      if ~nopar2
         oldpar2 = evalin('base', char(parname2));
         fbrack = strfind(char(parname2), '(');
      else
         parname2 = 'iteration no.';
         start2 = 1;
         nend2 = nsteps1;
      end;
%     ppar = start2;
      if isnan(g_grind.tstep)
        np=nsteps1*writing1+1; %assume each step one outcome
      else
        np=nsteps1*g_grind.tstep+1;
      end;          
      X1 = zeros(np*nsteps2,1);
      k=1;
      Res1 =  zeros(np*nsteps2,g_grind.statevars.dim);
      Y1 =  zeros(np*nsteps2,1);
      t1 =X1;
      for i = 1:nsteps2
         ppar = start2 + (i-1) * (nend2 - start2) / (nsteps2-1);
         if ~nopar2
            multassignin('base', char(parname2), ppar, fbrack);
         end;
         i_paranal(plotprev)
         close(gcf);
         X1(k:k+size(g_paranal.p,1),:) = [NaN; g_paranal.p];
         t1(k:k+size(g_paranal.p,1),:) = [NaN; g_paranal.t];
         Y1(k:k+size(g_paranal.p,1),:) = [NaN; ones(size(g_paranal.p)) * ppar];
         Res1(k:k+size(g_paranal.p,1),:) = [ones(1, size(g_paranal.Y, 2)) + NaN; g_paranal.Y];
         k=k+size(g_paranal.p,1)+1;
     end;
     if size(X1,1)>k
        X1=X1(1:k-1,:);
        t1=t1(1:k-1,:);
        Y1=Y1(1:k-1,:);
        Res1=Res1(1:k-1,:);
     end;
      if ~nopar2
         multassignin('base', char(parname2), oldpar2, fbrack);
      end;
      g_paranal2d.p = X1;
      g_paranal2d.p1 = Y1;
      g_paranal2d.Y = Res1;
      g_paranal2d.t = t1;
    end;
end;
%catch
%   if ~nopar1
%      multassignin('base', char(parname1), oldpar1, fbrack);
%   end;
%end;
if outputtype~=1
   outputlist=i_paranaldialog('outputlist');
   outputtype=lower(outputlist{outputtype});
   [pY,p]=catfun(outputtype,g_paranal2d.Y,[g_paranal2d.p,g_paranal2d.p1]);
   p1=p(:,2);
   p=p(:,1);
else
   p=g_paranal2d.p;
   p1=g_paranal2d.p1;
   pY=g_paranal2d.Y;
end;
if lines >=2
   [X,Y,Z]=makegrids(p,p1,pY);
end;
for i = 1:g_grind.statevars.dim
   H = i_makefig('paranal2d', i);
   if lines ==2
      contourf(X,Y,Z{i});
      colorbar;
   elseif lines==3
      surf(X,Y,Z{i});
      colorbar;
   else
   if lines == 0
      plot3(p, p1, pY(:, i), '.k');
      set(gca,'DrawMode', 'fast');
   elseif lines == 1
      plot3(p, p1, pY(:, i), '-k');
      set(gca,'DrawMode', 'fast');
   end;
   box on;
   end;
   hold on;
   set(H, 'Name', ['paranal2d of ', i_statevars_names(i)]);
   xlabel(parname1);
   ylabel(parname2);
   zlabel(i_statevars_names(i));
   shading flat;
end;
if nargout > 0
   Res = Res1;
   X = X1;
   Y = Y1;
end;

function [X,Y,Z]=makegrids(p,p1,Y1)
pp=sort(p);
pp(isnan(pp))=Inf;
x=pp(diff(pp)>0);
pp=sort(p1);
pp(isnan(pp))=Inf;
y=pp(diff(pp)>0);
[X,Y]=meshgrid(x,y);
x1=p(~isnan(p));
y1=p1(~isnan(p1));
z1=Y1(~isnan(Y1(:,1)),:);
Z=cell(1,size(z1,2));
for Col=1:size(z1,2)
   Z{Col}=griddata(x1,y1,z1(:,Col),X,Y);
end;

function multassignin(ws, name, V, fbrack)
%multassignin supports assignments to elements of arrays (e.g. name='A(1,2)')
%fbrack should be strfind(name,'(') (added for efficiency)
if isempty(fbrack)
   assignin(ws, name, V);
else
   temp = evalin(ws, name(1, fbrack(1) - 1));
   s = ['temp' name(fbrack(1):length(name)) ' = ' num2str(V) ';' ];
   eval(s);
   assignin(ws, name(1, fbrack(1) - 1), temp);
end;
