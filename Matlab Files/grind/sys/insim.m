%INSIM   Phenology and populatIoN SIMulator (Mols)
%  INSIM is a specific model system. The command INSIM reads a input text file and
%  creates a GRIND inifile from this text file. See further INSIM manual (in prep.).
%  Use the GRIND command 'time' to run the model.
%
%  Usage:
%  INSIM - user is prompted for a input text file. It creates an ini file with 
%  the same name.
%  INSIM FILE.TXT - enter filename on the command line
%  INSIM -sizes - show minimal sizes of boxcartrains for different deltat's
%  INSIM -check - check current sizes of boxcartrains (is done automatically 
%  after opening file)
%  INSIM -deltat - increase deltat to the maximum possible value
%
%  See also boxcartrain, use

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function [outp] = insim(filename)
global g_grind;
if nargin < 1
   [filename, pathname] = uigetfile('*.txt', 'Get input file for INSIM');
   if filename == 0
      disp('Cancel pressed');
      return;
   else
      filename=[pathname filename];
   end;
else
   if strncmpi(filename, '-s', 2) %-sizes option
      deltats = [g_grind.solver.opt.MaxStep 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 1];
      outp1 = zeros(size(g_grind.boxcar.names));
      if nargout == 0
         fprintf('delta-t:\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n', deltats(:));
      end;
      for j = 1:length(g_grind.boxcar.names)
         A = evalin('base', sprintf('t_%s', g_grind.boxcar.names{j}));
         devrate  = A(:, 2);
         cv_devrate = A(:, 3);
         [a, b] = i_statevarnos(g_grind.boxcar.names{j});
         s = sprintf('%s(%d)', g_grind.boxcar.names{j}, b - a + 1);
         deltat = deltats(1);
         outp1(j) = floor(min(1 ./ (1E-30 + deltat .* devrate + cv_devrate.^2)));
         for i = 1:length(deltats)
            deltat = deltats(i);
            s = (sprintf('%s\t%d', s, floor(min(1 ./ (1E-30 + deltat .* devrate + cv_devrate.^2)))));
         end;
         if nargout == 0
            fprintf('%s\n', s);
         end;
      end;
      if nargout == 1
         outp = outp1;
      end;
      return;
   elseif strncmpi(filename, '-c', 2) %-check option
      sizes = insim('-sizes');
      allok = 1;
      for j = 1:length(g_grind.boxcar.names)
         [a, b] = i_statevarnos(g_grind.boxcar.names{j});
         if sizes(j) < b - a + 1
            ButtonName = questdlg(sprintf('INSIM: boxcartrain of %s is probably too large (if all temperatures occur)\n\nDo you want to shrink it to %d?\n\n' ...
               ,g_grind.boxcar.names{j},sizes(j)),'INSIM','Yes','No','Cancel','Yes');
            switch ButtonName
             case 'Yes'
               setdimension(g_grind.boxcar.names{j}, sizes(j));
             case 'Cancel'
               return;
             otherwise
               allok = 0;
            end;
         end;
      end;
      if allok
         fprintf('The sizes of the boxcartrains are ok (time step=%g)\n', g_grind.solver.opt.MaxStep);
      end;
      return;
   elseif strncmpi(filename, '-d', 2) %-deltat option
      deltats=zeros(length(g_grind.boxcar.names),1);
      for j = 1:length(g_grind.boxcar.names)
         A = evalin('base', sprintf('t_%s', g_grind.boxcar.names{j}));
         [a, b] = i_statevarnos(g_grind.boxcar.names{j});
         N = b - a + 1;
         gamm = 1 / N;
         devrate  = A(:, 2);
         cv_devrate = A(:, 3);
%         s = [g_grind.boxcar.names{j}];
         deltats(j) = min((1 - (N .* cv_devrate.^2)) .* gamm ./ (devrate + 1e-30));
      end;
      mindeltat = min(deltats);
      if mindeltat-1E-10 < g_grind.solver.opt.MaxStep
         disp('The time step cannot be increased with the current boxcar sizes (see INSIM -S'); 
         return;
      end;
      ButtonName=questdlg(sprintf('Time step can be set to %g, ok to change?',mindeltat),'INSIM','OK','Other value', 'Cancel','OK');
      switch ButtonName
       case 'OK'
         g_grind.solver.opt.MaxStep = mindeltat;
         disp('integration time step changed')
       case 'Other value'
         answ=inputdlg('Value for delta-t','INSIM',1,{ num2str(mindeltat) });
         if ~isempty(answ)
            g_grind.solver.opt.MaxStep = str2double(answ{1}); 
            disp('integration time step changed')
         end;
      end;
      
      return;
   end;
end;
[pathname, outname,ext] = fileparts(filename);
if isempty(ext)
   filename = [filename '.txt'];
end;
[p, outname] = fileparts(filename);
outname = fullfile(pathname, [outname '.ini']);
%outname=inputdlg('Name of the output file:','INSIM',1,{outname});
A = tabread(filename,  -1);
if isempty(A)
   error('GRIND:insim:EmptyFile','INSIM: file %s is empty', filename);
end;
[NRows, NCols] = size(A);
fout = fopen(outname, 'w');
%try
fprintf(fout, '%%model\n');
i = 1;
while ~isempty(strfind(A{i, 1},'%'))
   fprintf(fout, '%s\n',A{i, 1});
   i = i + 1;
end;
i = findline(1, A, NRows, NCols, '[overview]');
if isempty(i)
   error('GRIND:insim:overview','INSIM: [overview] not found');
end;
NA = [];
while (i <= NRows) && isempty(NA)
   i = i + 1;
   if ischar(A{i, 4}) %Classes column should be numeric
      NA = str2num(A{i, 4}); %#ok
   else
      NA = [];
   end;
end;
j = 0;
Symbols = cell(20,1); Input = cell(20,1); Classes = zeros(20,1);
while (i < NRows) && ~isempty(A{i, 2})
   j = j + 1;
   Symbols{j} = A{i, 2};
   Input{j} = A{i, 3};
   if ~isempty(Input{j})&&(strncmp('**',Input{j},2)||strcmp('-',Input{j}))
      Input{j} = 'NaN';
   end;
   Classes(j) = str2double(A{i, 4});
   i = i + 1;
end;
if j==0
   error('GRIND:insim:NoStatevars','INSIM: no state variables defined');
elseif j<20
  Symbols = Symbols(1:j); Input = Input(1:j); Classes = Classes(1:j);
end;
for i = 1:length(Input)
   s = symvar(Input{i});
   for n = 1:length(s)
      if isempty(getspecindex(s{n}, Symbols))
         error('GRIND:insim:UnknownSymb','INSIM: input of %s (%s) has unknown symbol', Input{i}, Symbols{i})
      end;
   end;
end;
for i = 1:NRows
   s =  lower(A{i, 1});
   insimcomms = '[overview]\n\t[reproduction]\n\t[temperature]\n\t[table]\n\t[initial conditions]\n\t[stepsize]';
   if ~isempty(s) && (s(1) == '[') && isempty(strfind(insimcomms, s))
      error('GRIND:insim:UnknownComm','INSIM: command ''%s'' unknown\nValid commands include:\n\t%s',s,sprintf(insimcomms));
   end;
end;
i = findline(1, A, NRows, NCols, '[reproduction]');
Repr = cell(20, 3);
j = 0;
while ~isempty(i)
   k = 2;
   while (k < NCols) && ~isempty(A{i, k})
      R = A{i, k};
      j = j + 1;
      Repr{j, 1} = R;
      f = strfind(R, '->');
      if ~isempty(f)
         Repr{j, 2} = getspecindex(R(1:f - 1),Symbols);
         Repr{j, 3} = getspecindex(R(f + 2:end),Symbols);
         if isempty(Repr{j, 2})
            error('GRIND:insim:stageUnknown','INSIM Reproduction: Stage (%s) not found', R(1:f - 1));
         end;
         if isempty(Repr{j, 3})
            error('GRIND:insim:stageUnknown','INSIM Reproduction: Stage (%s) not found', R(f + 2:end));
         end;
      else
         error('GRIND:insim:NoArrow','INSIM: No "->" sign in reproduction field');
      end;
      k = k + 1;
   end;
   i = findline(i + 1, A, NRows, NCols, '[reproduction]');
end;
Repr = Repr(1:j, :);
i = findline(1, A, NRows, NCols, '[temperature]');
fprintf(fout, 'definepars Tcrit;\n');
if ~isempty(i) && strcmpi(A(i, 2), 'min/max')
   Tminmax = 1;
   fprintf(fout, 'defextern Tmax 15 ''-c'' ''-f'';\n');
   fprintf(fout, 'defextern Tmin 10 ''-c'' ''-f'';\n');
   fprintf(fout, 'T=Tdev+((Tmin+Tmax)-(Tmax-Tmin)*cos(2*pi*t))*0.5;\n');
else
   Tminmax = 0;
   fprintf(fout, 'defextern Tmean 15 ''-c'' ''-f'';\nT=Tdev+Tmean;\n');
end;
for j = 1:length(Symbols)
   i = findline(1,A, NRows, NCols, '[table]', Symbols{j});
   if ~isempty(i)
      fprintf(fout,'p_%s=parlookup(t_%s,T,1);\n', Symbols{j}, Symbols{j});
   else
      btn=questdlg(sprintf('INSIM: warning: no lookup table found for %s, continue?', Symbols{j}),'INSIM','No');
      if strcmp(btn, 'No')
         return;
      end;
   end;
end;
for j = 1:size(Repr, 1)
   fprintf(fout,'r_%s_%s=parlookup(tr_%s_%s,T,1);\n', Symbols{Repr{j,2}},Symbols{Repr{j,3}},Symbols{Repr{j,2}},Symbols{Repr{j,3}});
end;
for j = 1:length(Symbols)
   if Classes(j) > 1
      fprintf(fout,'b_%s=boxcartrain(%s,p_%s(1),p_%s(2));\n', Symbols{j}, Symbols{j}, Symbols{j}, Symbols{j});
   end;
end;
for j = 1:length(Symbols)
   if Classes(j) > 1
      if strcmpi(Input{j}, 'NaN')
         fprintf(fout,'%s(1:%d)''=b_%s.flow -p_%s(3).*%s', Symbols{j}, Classes(j),...
            Symbols{j}, Symbols{j}, Symbols{j});
      else
         fprintf(fout,'%s(1:%d)''=b_%s.flow + boxcarinflow(%s,%s)-p_%s(3).*%s', Symbols{j}, Classes(j),...
            Symbols{j}, Symbols{j}, symreplace(Input{j},'b_%s.outflow',Symbols), Symbols{j}, Symbols{j});
      end;
   else
      fprintf(fout,'%s(1:%d)''=%s -p_%s(3).*%s', Symbols{j}, Classes(j),...
         symreplace(Input{j},'b_%s.outflow',Symbols), Symbols{j}, Symbols{j});
   end;
   for n =  1:size(Repr, 1)
      if Repr{n, 3} == j
         fprintf(fout,'+ boxcarinflow(%s,r_%s_%s*%s)',Symbols{Repr{n,3}},Symbols{Repr{n,2}}, ...
            Symbols{Repr{n, 3}}, Symbols{Repr{n, 2}});
      end;
   end;
   fprintf(fout, ';\n');
end;
fprintf(fout, '%%commands\nsimtime 0 366 183\nt=1;\nTdev=0;\nTcrit=10;\n');
for j = 1:length(Symbols)
   fprintf(fout,'%s=zeros(size(%s));\n', Symbols{j}, Symbols{j});
end;
i = findline(1,A, NRows, NCols, '[initial conditions]');
if ~isempty(i)
   i = i + 1;
   while (i <= NRows) && ~isempty(A{i, 1})
      fprintf(fout,'%s=%s;\n',A{i, 1},A{i, 2});
      i = i + 1;
   end;
else
   %default initial conditions
   fprintf(fout, '%s(1)=100;\n', Symbols{1});
end;
for j = 1:length(Symbols)
   i = findline(1, A, NRows, NCols, '[table]',Symbols{j});
   NA = [];
   while (i <= NRows) && isempty(NA)
      i = i + 1;
      if ischar(A{i, 1})
         NA = str2num(A{i, 1}); %#ok
      else
         NA = [];
      end;
   end;
   T = zeros(20,4);
   %i = i + 2;
   k = 1;
   maxdevr = 0;
   while (i <= NRows) && ~isempty(A{i, 1})
      T(k, 1) = str2double(A{i, 1});
      if ~isempty(A{i, 2})
         devr = str2double(A{i, 2});
      else
         devr = [];
      end;
      if devr > maxdevr
         maxdevr = devr;
      end;
      if ~isempty(A{i, 3})
         sd = str2double(A{i, 3});
      else
         sd = [];
      end;
      if isempty(devr) || isnan(devr) || (devr == 0)
         T(k, 2) = NaN;
         devr = NaN;
      else
         T(k, 2) = 1 / devr;
      end;
      if isempty(sd) || isnan(sd) || (devr == 0)
         T(k, 3) = NaN;
      else
         T(k, 3) = sd / devr; % coefficient of variation EXPRESSED IN TIME
      end;
      mort = str2double(A{i, 4});
      if isempty(mort) || isnan(mort) 
         T(k, 4) = NaN;
      else
         T(k, 4) = mort; % coefficient of variation EXPRESSED IN TIME
      end;

      i = i + 1;
      k = k + 1;
   end;
   T=T(1:k-1,:);
   if maxdevr < 1 % if devr < 1 we assume that all things are expressed as rates
      T(:, 2) = 1 / T(:, 2);
      T(:, 3) = 1 / T(:, 3); % Should be expressed in TIME!!!
   end;
   s='\nt_%s=[';
   if ~isempty(T)
      T = sortrows(T, 1);
      for k = 1:size(T, 1)
         s = sprintf('%s%g,%g,%g,%g;\n',s,T(k,1),T(k,2),T(k,3),T(k,4));
      end;
      f = strfind(s, ';');
      if ~isempty(f)
         s(f(end)) = ']';
      else
         s = sprintf('%s]',s);
      end;
      fprintf(fout, s, Symbols{j});
   end;
end;
for j = 1:size(Repr, 1)
   i = findline(1, A, NRows, NCols, '[reproduction]',Repr{j,1});
   NA = [];
   while (i <= NRows) && isempty(NA)
      i = i + 1;
      if ischar(A{i, 1})
         NA = str2num(A{i, 1}); %#ok
      else
         NA = [];
      end;
   end;
   T = zeros(20, NCols);
   h = i  - 1;
   % i = i + 2;
   k = 0;
   while (i <= NRows) && ~isempty(A{i, 1})
      n = 1;
      k = k + 1;
      while (n <= NCols) && ~isempty(A{h, n})
         T(k, n) = str2double(A{i, n});
         n = n + 1;
      end;
      i = i + 1;
   end;
   T=T(1:k, :);
   s='\ntr_%s_%s=[';
   if size(T, 2) - 1 ~= Classes(Repr{j, 2})
      error('GRIND:insim:ReproError','Reproduction error: number of classes in reproduction %s unequal (%d) to defined classes (%d)' ...
         , Symbols{Repr{j, 2}}, size(T, 2) - 1, Classes(Repr{j, 2}));
   end;
   T = sortrows(T, 1);
   for k = 1:size(T, 1)
      s = sprintf('%s %g', s, T(k, 1));
      for n = 2:size(T, 2)
         s = sprintf('%s,%g',s,T(k,n));
      end;
      s = sprintf('%s;\\n',s);
   end;
   s(length(s) - 2) = ']';
   fprintf(fout, s, Symbols{Repr{j,2}},Symbols{Repr{j,3}});
end;
i = findline(1, A, NRows, NCols, '[stepsize]');
if ~isempty(i)
   fprintf(fout,'solver(''euler'',%s);\n',A{i,2});
else
   fprintf(fout,'solver(''euler'',0.1);\n');
end;

s = 'out -1 ';
for j = 1:length(Symbols)
   if Classes(j) > 1
      s=sprintf('%soutf(''sum'',''%s'') ',s,Symbols{j});
   else
      s =sprintf('%s',s,Symbols{j});
   end;
end;
i = findline(1, A, NRows, NCols, '[temperature]');
if ~isempty(i)
   i = i + 1;
   while (i <= NRows) && ~isempty(A{i, 2})
      if ~strcmpi(A{i, 2}, 'nan')
         p = fileparts(A{i, 2});
         if isempty(p)
            loadname=fullfile(pathname,A{i, 2});
         else
            loadname=A{i,2};
         end;
         fprintf(fout, 'loaddata ''%s''; ', loadname);
      end;
      i = i + 1;
   end;
end;
if Tminmax
   fprintf(fout, '\n%s\nout -2 Tmin Tmax\nout -3 tempsum(''Tmax'',Tcrit,365)\n',s);
else
   fprintf(fout, '\n%s\nout -2 Tmean\nout -3 tempsum(''Tmean'',Tcrit,365)',s);
end;


fclose(fout);
use(outname)
disp('***********************************');
disp('INSIM: GRIND model ready for use');
disp('time = make a time plot');
disp('vectplot -b = show the contents of the boxcartrains');
disp('commands = help on other GRIND commands');
insim('-check');
%catch
%   fclose(fout);
%end;
function res = getspecindex(s, Symbols)
for j = 1:length(Symbols)
   if strcmp(s, Symbols{j})
      res = j;
      return;
   end;
end;
res = [];

function res = symreplace(s, repstr, Symbols)
res = s;
for j = 1:length(Symbols)
   L = length(Symbols{j});
   f = strfind(res, Symbols{j});
   for i = length(f):-1:1
      if f(i) + L > length(res)
         res = [res(1:f(i) - 1) sprintf(repstr, Symbols{j})];
      elseif isempty(strfind('abcdefghijklmnopqrstuvwABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_', res(f(i) + L)))
         res = [res(1:f(i) - 1) sprintf(repstr, Symbols{j}) res(f(i) + L:length(res))];
      end;
   end;
end;

function i = findline(istart,A, NRows, NCols, m, m2)
if nargin < 6
   m2 = [];
end;
i = istart;
found = 0;
while ~found
   while (i <= NRows) && ~strcmpi(A{i, 1}, m)
      i = i + 1;
   end;
   if i > NRows
      i = [];
      found = 1;
   elseif ~isempty(m2)
      j = 2;
      while (j <= NCols) && ~strcmpi(A{i, j}, m2)
         j = j + 1;
      end;
      if  j <= NCols
         found = 1;
      else
         i = i + 1;
      end;
   else
      found = 1;
   end;
end;

