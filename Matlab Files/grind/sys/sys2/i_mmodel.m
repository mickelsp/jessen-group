%i_mmodel internal GRIND function to create an odefile
% g_grind.model, g_grind.commands  should contain resp the model and the funcs
function i_mmodel(mode1, comman1, inifil1, schem)
if (nargin>=3) || (nargin==1)
   evalin('base','initgrind');
end;
global g_grind;
if nargin >= 3
   i_startup;
   initgrind;
   g_grind.model = mode1;
   g_grind.commands = comman1;
   g_grind.inifile = inifil1;
   if nargin == 4
      g_grind.scheme = schem;
   elseif ~isfield('scheme', g_grind)
      g_grind.scheme = {};
   end;
   clear('mode1','comman1','inifil1');
elseif nargin == 1
   i_startup;
   g_grind.model = mode1.model;
   g_grind.commands = mode1.commands;
   g_grind.inifile = mode1.inifile;
   if isfield(mode1, 'scheme')
      g_grind.scheme = mode1.scheme;
   end;
end;
for g_l_i = 1:length(g_grind.model)
   g_l_s = strtrim(g_grind.model{g_l_i});
   if ~isempty(g_l_s) && (g_l_s(1) ~= '%')
      g_l_f2 = strfind(g_l_s, '/dt');
      if ~isempty(g_l_f2)
         g_l_f1=strfind(g_l_s, '=');
         if ~isempty(g_l_f1) && (g_l_f1(1) > g_l_f2(1))
            g_grind.model{g_l_i}=[g_l_s(2:g_l_f2(1)-1) '''' g_l_s(g_l_f2(1)+3:length(g_l_s))];
         end;
      end;
   end;
end;
g_grind.model = removepoints(g_grind.model, 0);
if ~isempty(g_grind.commands)
   g_grind.commands = removepoints(g_grind.commands, 1);
   for g_l_h = 1:size(g_grind.commands, 2)
      if isempty(strfind(g_grind.commands{g_l_h}, '%'))
         g_l_c = strtrim(g_grind.commands{g_l_h});
         g_l_l = length(g_l_c);
         while (g_l_l > 1) && (g_l_c(g_l_l) == ' ')
            g_l_l = g_l_l - 1;
         end;
         if (g_l_l > 0) && (g_l_c(g_l_l) ~= ';')
            g_grind.commands{g_l_h} = [g_grind.commands{g_l_h} ';'];
         end;
      end;
   end;
end;
%global g_grind;
if ~isempty(g_grind.model) && iscell(g_grind.model{1})
   g_grind.model = g_grind.model{1};
end;
if isempty(g_grind.pars)
   i_fillpars;
end
if ~isempty(g_grind.pars)
   evalin('base', i_globalstr(g_grind.pars,[],'clear')); % clear the parameters first to avoid warnings
   evalin('base', i_globalstr(g_grind.pars));
end;
for g_l_i = length(g_grind.model):-1:1
   if isempty(strtrim(g_grind.model{g_l_i}))
      g_grind.model = {g_grind.model{1:g_l_i - 1} g_grind.model{g_l_i + 1:size(g_grind.model, 2)}};
   end;
end;
clear(g_grind.odefile);
gmex('-d');
oldpath = pwd;
cd(grindpath);
aodefile = [g_grind.odefile '.m'];
g_grind.statevars.names = cell(1, size(g_grind.model, 2));
g_rhs = cell(1, size(g_grind.model, 2));
funcs = cell(1, size(g_grind.model, 2));
localfuncs = cell(1, size(g_grind.model, 2));
ilocfun=0;
irhs = 0;
ifun = 0;
g_grind.solver.haslag = 0;
infunction = 0;
LF = sprintf('\n');
dimtot = -1;
diffchar = [];
dim2 = [];
g_grind.funcnames.names = {};
for i1 = 1:size(g_grind.model, 2)
   s = strtrim(g_grind.model{i1});
   f = strfind(s, 'boxcartrain(');
   if ~isempty(f)
      % with boxcartrain the default solver is Euler
      g_grind.solver.name = 'euler';
      g_grind.solver.opt.MaxStep = 0.3;
      f2 = strfind(s, ',');
      f3 = strfind(s, '=');
      s=['[' s(1:f3(1)-1) ',' s(f(1)+12:f2(1)-1) ']' s(f3(1):f(1)+11) '''' s(f(1)+12:f2(1)-1) ''',' s(f(1)+12:length(s))];
   end;
   f = strfind(s, 'dwiener(');
   if ~isempty(f)
      % with Ito stochastic differential equations the default solver is Euler
      g_grind.solver.name = 'euler';
      g_grind.solver.opt.MaxStep = 0.1;
      f2 = strfind(s(f(1) + 8:end), '(');
      f3 = strfind(s(f(1) + 8:end), ')');
      if isempty(f2) || (f2(1) > f3(1))
         s = [s(1:f3(1) + f(1) + 7 - 1),',t',s(f3(1) + f(1) + 7:end)]; %add ,t to the equation
      else
         %TO DO analyse STACKED ()
      end;
   end;
   f = [strfind(s, 'rand(') strfind(s,'randn(') strfind(s, 'drawuniform(')];
   if ~isempty(f)
      warning('GRIND:model:euler','Stochastic functions need a fixed time solver\n%s',...
          'Solver changed to euler with step of 0.1')
      g_grind.solver.name = 'euler';
      g_grind.solver.opt.MaxStep = 0.1;
   end;
   
   if strncmpi(s, 'setevent(', 9)
      error('GRIND:model:setevent','"setevent" should be in the lower parameters panel');
   end;
   if infunction
      ilocfun=ilocfun+1;
      localfuncs{ilocfun} = s;
%      fwrite(ffunc, [makearrayop(s) LF]);
      if strncmp(s, 'return', 6)
         infunction = 0;
%         fclose(ffunc);
      end;
   elseif strncmp(s, 'function ', 9)
%       f=strfind(s, '=');
%       if isempty(f)
%          f = strfind(s, ' ');
%       end;
%       if ~isempty(f)
%          f2 = strfind(s, '(');
%          if isempty(f2)
%             f2 = length(s);
%          end;
%          fname = strtrim(s(f(1) + 1:f2 - 1));
%          if strfind(fname, ';')
%             fname = strtrim(fname(1:strfind(fname, ';') - 1));
%          end;
%          [ffunc,message]=fopen([fname '.m'],'w');
%          if ffunc < 0
%             errordlg(message);
%             error('GRIND:model:FileError',message);
%          end;
%       end;
%       fwrite(ffunc, [s LF]);
       ilocfun=ilocfun+1;
      localfuncs{ilocfun} = s;
      infunction = 1;
   elseif strncmp(s,'if ',3)||strncmp(s,'for ',4)||strncmp(s,'while ',6)||strncmp(s,...
         'case ',5)||strncmp(s,'otherwise',9)||strncmp(s,'switch ',6)
      ifun = ifun + 1;
      funcs{ifun} = s;
   elseif ~isempty(s) && (s(1) ~= '%')
      s=strrep(s,'( ','(');
      feq=strfind(s, '=');
      if ~isempty(feq) && any(feq(1)+1==feq) %~isempty(find(feq(1)+1==feq))
         feq = feq(3:end);
      end;
      feq = min(feq);
      fdiff=strfind(s, '''');
      if isempty(fdiff) || (mod(length(fdiff), 2) == 0)
         %     flogical=~(isempty(strfind(s, '=='))&& isempty(strfind(s,'<')) && isempty(strfind(s, '>')));
         %     if ~flogical
         flogical=min([strfind(s, '==') strfind(s, '<') strfind(s, '>') strfind(s, '~=')]);
         fdiffer = strfind(s, '(');
         if ~isempty(flogical) && (flogical < feq)
            fdiffer = [];
         end;
         if  ~isempty(fdiffer) && ~isempty(feq) && (fdiffer(1) < feq)
            i = 1;
            while (i <= length(fdiffer)) && (fdiffer(i) < feq)
               if ~ismember(s(fdiffer(i)+1), '01234567890-+: ')
                  diffchar = s(fdiffer(i):fdiffer(i) + 1);
                  fdiffer = strfind(s, diffchar);
                  i = 1;
                  break;
               end;
               i = i + 1;
            end;
            if (i > 1) || (fdiffer(1) > feq)
               fdiffer = [];
            end;
            if ~isempty(fdiffer) && (fdiffer(1) < feq) && ~isempty(diffchar)&&(diffchar(2) ~= 't') && ~isempty(g_grind.pars)
               l_k = 1;
               %   i = i - 1;
               while (l_k <= length(g_grind.pars)) && ~strcmp(g_grind.pars{l_k}, diffchar(2))
                  l_k = l_k + 1;
               end;
               if l_k <= length(g_grind.pars)
                  if length(g_grind.pars) <= 1
                     g_grind.pars = {};
                  else
                     g_grind.pars = {g_grind.pars{1:l_k - 1} g_grind.pars{l_k + 1:length(g_grind.pars)}};
                  end;
               end;
            end;
            %   end;
         end;
         if ~isempty(fdiffer) && ~isempty(feq) && (fdiffer(1) > feq)
            fdiffer = [];
         end;
      else
         fdiffer = [];
      end;
      fcolon = strfind(s, ':');
      if ~isempty(fcolon) && (~isempty(fdiffer) || ~isempty(fdiff)) && (feq(1) > fcolon(1))
         if dimtot == -1
            if irhs > 0
               error('GRIND:model:MixScalarMatrix','Mix of vector/matrix notation and scalar notation not supported, use e.g. X(1:1)'' instead of X''');
            end;
            dimtot = 0;
            dim = [];
            dim2 = [];
            start = [];
         end;
         fbrack1 = strfind(s, '(');
         fbrack2 = strfind(s, ')');
         if (length(fbrack1) > 1) && ~isempty(fdiffer) && (fbrack1(1) == fdiffer(1))
            fbrack2 = min(fbrack2(2:length(fbrack2)));
            fbrack1 = min(fbrack1(2:length(fbrack1)));
         else
            fbrack2 = min(fbrack2);
            fbrack1 = min(fbrack1);
         end;
         if ~isempty(fbrack1)
            if length(fcolon) == 1
               d1 = str2double(s(fbrack1 + 1:fcolon(1) - 1));
               d2 = str2double(s(fcolon(1) + 1:fbrack2 - 1));
               d = d2 - d1 + 1;
               dim = [dim d];
               start = [start d1 - 1];
               dimtot = dimtot + d;
               s = [s(1:fbrack1 - 1) s(fbrack2 + 1:length(s))]; 
            else
               s1 = s(fbrack1 + 1:fbrack2 - 1);
               fcomma = strfind(s1, ',');
               kcolon = strfind(s1, ':');
               d1 = str2double(s1(1:kcolon(1) - 1));
               d2 = str2double(s1(kcolon(1) + 1:fcomma(1) - 1));
               d3 = str2double(s1(fcomma(1) + 1:kcolon(2) - 1));
               d4 = str2double(s1(kcolon(2) + 1:length(s1)));
               dimen1 = d2 - d1 + 1;
               dimen2 = d4 - d3 + 1;
               dim = [dim dimen1];
               dim2 =  [dim2 dimen2];
               start = [start d1 - 1];
               dimtot = dimtot + dimen1 * dimen2;
               s = [s(1:fbrack1 - 1) s(fbrack2 + 1:length(s))]; 
            end;
            if ~isempty(fdiffer)
               fdiffer = strfind(s, diffchar);
            end;
         end;
      else
         if dimtot ~= -1
            error('GRIND:model:MixScalarMatrix','Mix of vector/matrix notation and scalar notation not supported, use e.g. X(1:1)'' instead of X''');
         end;
      end;
      
      %    k1=strfind(s, '=');
      %   k=strfind(s, '''');
      if ~isempty(fdiff) && ~isempty(feq) && (feq > fdiff(1))
         j1 = 1;
         if (length(fdiff) > 1) && (feq > fdiff(2))
            error('GRIND:model:SyntaxError',['Syntax error in: ' s(1:fdiff) LF s ]);
         end;
         while (s(j1)~='''')&&(s(j1)~=' ')&&(s(j1)~='=')
            j1 = j1 + 1;
         end;
         irhs = irhs + 1;
         g_grind.statevars.names{irhs} = s(1:j1 - 1);
         while (s(j1) ~= '=')
            j1 = j1 + 1;
         end;
         g_rhs{irhs} = s(j1 + 1:size(s, 2));
      else
         if ~isempty(fdiffer) && ~isempty(feq) && (fdiffer(1) < feq)
            %difference equation
            k1 = fdiffer(1) + 2;
            k2 = k1;
            while s(k2) ~= ')'
               k2 = k2 + 1;
            end;
            if k2 == k1
               difford = 0;
            else
               difford = str2double(s(k1:k2 - 1));
            end;
            g_grind.solver.name = 'i_differ';
            g_grind.solver.isdiffer = 1;
            g_grind.ndays = 100;
            j1 = 1;
            while (s(j1)~='(')&&(s(j1)~=' ')&&(s(j1)~='=')
               j1 = j1 + 1;
            end;
            irhs = irhs + 1;
            g_grind.statevars.names{irhs} = s(1:j1 - 1);
            while (s(j1) ~= '=')
               j1 = j1 + 1;
            end;
            g_rhs{irhs} = s(j1 + 1:size(s, 2));
            %        g_rhs{irhs}=[s '-' char(g_grind.statevars.names{irhs})];
            if ~isempty(strfind(g_rhs{irhs}, 'lag('))
               g_grind.solver.haslag = 1;
               error('GRIND:model:LagInDiffer','time lag not supported in difference equations');
            end;
         else
            ifun = ifun + 1;
            funcs{ifun} = s;
            f=strfind(s, '=');
            if (length(f)==1) && isempty(strfind(s, '['))
               funcname = strtrim(s(1:f(1) - 1));
               g_grind.funcnames.names = [g_grind.funcnames.names {funcname}];
            end;
         end;
      end;
   end;
end;
%functions=functions(1:jfunc);
if irhs == 0
   error('GRIND:model:NoEquation','There is not any difference/differential equation entered');
end;
g_rhs = g_rhs(1:irhs);
for i = 1:irhs
   g_rhs{i} = extractlags(g_rhs{i});
end;
funcs = funcs(1:ifun);
%for i = 1:ifun
%   funcs{i} = extractlags(funcs{i});
%end;
if g_grind.solver.isdiffer
   [g_rhs,difford2] = clearlist(g_rhs, irhs,diffchar);
   [funcs,dord2] = clearlist(funcs, irhs,diffchar);
   if isempty(difford2)
      difford2 = dord2;
   elseif ~isempty(dord2) && (dord2 ~= difford2)
      error('GRIND:model:MixOfOrders','Error in difference equation: mix of different orders not allowed');
   end;
   if ~isempty(difford2)
      if (difford2 >= difford)
         error('GRIND:model:OrderTooSmall','Difference equation cannot be of an order <1')
      else
         g_grind.solver.iters = difford - difford2;
      end;
   end;
end;
if (dimtot > 0) && (irhs > 0)
   statvars = g_grind.statevars.names(1:irhs);
   g_grind.statevars.vector = 1;
   %rmfield(g_grind.statevars,'names')
   g_grind.statevars.names = {}; %cell(1, dimtot);
   g_grind.statevars.dims = cell(1, irhs);
   g_grind.statevars.vectnames = statvars;
   d = 0;
   for j = 1:length(statvars)
      %vectnames{j}=statvars{j};
      clear(statvars{j});
      v.dim1 = dim(j);
      if isempty(dim2)
         v.dim2 = 1;
         %  for i = 1:dim(j)
         %     g_grind.statevars.names{d + i} = sprintf('%s(%d)', statvars{j}, start(j) + i);
         %  end;
         d = d + dim(j);
      else
         v.dim2 = dim2(j);
         for l_k = 1:dim2(j)
            %nbase = d + (l_k - 1) * dim(j);
            %    for i = 1:dim(j)
            %       g_grind.statevars.names{nbase + i} = sprintf('%s(%d,%d)',statvars{j},start(j) + i,start(j) + l_k);
            %   end;
         end;
         d = d + dim(j) * dim2(j);
      end;
      g_grind.statevars.dims{j} = v;
      %g_grind.statevars = g_grind.statevars(1:dim);
   end;
   if sum(start) > 0
      %no matrix notation
      dimtot = 0;
   end;
   g_grind.statevars.dim = dimtot;
else
   dimtot = 0;
   dim = ones(irhs, 1);
   g_grind.statevars.names = g_grind.statevars.names(1:irhs);
   g_grind.statevars.dim = irhs;
end;
s = '';
g_grind.funcs = '';
g_grind.permanent = {};
for i1 = 1:ifun
   if (strncmp(funcs{i1}, 'defextern ', 10)||strncmp(funcs{i1}, 'defextern(', 10))
      s1 = defextern(stripcomments(char(funcs{i1})));
      g_grind.funcs = [g_grind.funcs s1 LF];
   elseif (strncmp(funcs{i1}, 'defpermanent ', 13)||strncmp(funcs{i1}, 'defpermanent(', 13))
      s1 = stripcomments(char(funcs{i1}));
      defpermanent('-i_mmodel', s1(14:end));
   elseif ~strncmp(funcs{i1}, 'global ', 7)&&~(strncmp(funcs{i1}, 'definepars ', 11)||strncmp(funcs{i1}, 'definepars(', 11))
      if strncmp(funcs{i1},'if ',3)||strncmp(funcs{i1},'end',3)||strncmp(funcs{i1},'else',4)||...
            strncmp(funcs{i1},'for ',4)||strncmp(funcs{i1},'while ',6)||strncmp(funcs{i1},...
            'case ',5)||strncmp(funcs{i1},'otherwise',9)||strncmp(funcs{i1},'switch ',6)
         g_grind.funcs = [g_grind.funcs stripcomments(char(funcs{i1})) LF];
         s = '  ';
         if strncmp(funcs{i1}, 'end', 3)
            s = '';
         end;
      else
         g_grind.funcs = [g_grind.funcs s stripsemicolon(stripcomments(char(funcs{i1}))) ';' LF];
      end;
   end;
end;
g_grind.funcs = sortfuncs(g_grind.funcs); % analysis of dependence of auxiliary variables!
funcs = g_grind.funcs;
funcs = extractlags(funcs);
for i1 = 1:length(g_rhs)
   s = g_rhs{i1};
   
   f = strfind(s, 'implicitdisperse(');
   if isfield(g_grind, 'implicdisp')
      ivar = length(g_grind.implicdisp);
   else
      ivar = 0;
   end;
   for ii = length(f):-1:1
      % with implicitdisperse the default solver is Euler
      g_grind.solver.name = 'Euler';
      g_grind.solver.opt.MaxStep = 0.3;
      f1 = f(ii) + 17;
      f2 = f1;
      stack = 1;
      while (f2 < length(s)) && (stack > 0)
         if (s(f2) == '(')
            stack = stack + 1;
         elseif (s(f2) == ')')
            stack = stack - 1;
         end;
         if (stack==1) && (s(f2)==',')
            comm = f2;
         end;
         f2 = f2 + 1;
      end;
      f2 = f2 - 1;
      avar = s(f1:comm - 1);
      dispers = s(comm + 1:f2 -  1);
      ivar = ivar + 1;
      implicitdisperse(ivar, avar, dispers); %initiation
      s = [s(1:f(ii) - 1) sprintf('implicitdisperse(%d,%s)',ivar,avar) s(f2 + 1:end)];
   end;
   g_rhs{i1}  = s;
end;


[fid] = fopen(aodefile, 'wt');
if fid < 0
   cd(oldpath);
   [fid,message] = fopen(aodefile, 'wt');
   if fid > 0
      warning('GRIND:model:writetosys','Cannot write to GRIND/SYS, using current directory instead');
   else
      errordlg(message);
      error('GRIND:model:FileError',message);
   end;
end;
fwrite(fid, ['%function created by GRIND' LF]);
if g_grind.solver.haslag
   fwrite(fid,sprintf('function g_X2=%s(t,g_X1,Lags)\n',g_grind.odefile));
else
   fwrite(fid,sprintf('function g_X2=%s(t,g_X1)\n', g_grind.odefile));
end;
if ~isempty(g_grind.pars)
   fwrite(fid, [i_globalstr(g_grind.pars) LF]);
end;
gl = 0;
if ~isempty(g_grind.permanent) || g_grind.solver.haslag
   fprintf(fid, 'global g_grind;\n');
   gl = 1;
end;
if g_grind.solver.haslag
   fwrite(fid,['if nargin==2,Lags=repmat(g_X1,length(g_grind.lags));end;' LF]);
end
if ~isempty(g_grind.permanent)
   fprintf(fid, 'i_updatepermanent(t);\n');
end;
if (dimtot > 0) || g_grind.statevars.vector
   if dimtot == 0 %not in matrix notation
      % v.name = g_grind.statevars.names{1}(1:strfind(g_grind.statevars.names{1}, '(') - 1);
      dimtot = irhs;
      v.dim1 = dimtot;
      v.dim2 = 1;
      g_grind.statevars.vector = 1;
      g_grind.statevars.dims = {v};
      g_grind.statevars.dim = dimtot;
      g_grind.statevars.vectnames = {g_grind.statevars.names{1}(1:strfind(g_grind.statevars.names{1}, '(') - 1)};
      statvars = g_grind.statevars.names;
   end;
   if ~gl
      fprintf(fid, 'global g_grind;\n');
   end;
   fprintf(fid,'g_X2=zeros(g_grind.statevars.dim,1);\n');
   d = 1;
   for j = 1:length(dim)
      if dimtot == 0
         g_grind.statevars.dims{j}.from = d;
         g_grind.statevars.dims{j}.to = d;
         fprintf(fid, '%s = g_X1(g_grind.statevars.dims{%d}.from);\n',statvars{j},d);
         d = d + 1;
      elseif isempty(dim2)
         g_grind.statevars.dims{j}.from = d;
         g_grind.statevars.dims{j}.to = d + dim(j) - 1;
         fprintf(fid, '%s = g_X1(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to);\n',statvars{j},j,j);
         d = d + dim(j);
      else
         dd = dim(j) * dim2(j);
         g_grind.statevars.dims{j}.from = d;
         g_grind.statevars.dims{j}.to = d + dd - 1;
         fprintf(fid, '%s = reshape(g_X1(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to),g_grind.statevars.dims{%d}.dim1,g_grind.statevars.dims{%d}.dim2);\n',statvars{j},j,j,j,j);
         d = d + dd;
      end;
   end;
else
   if irhs > 3
      fprintf(fid,'g_X2=zeros(%d,1);\n',irhs);
   end;
   for d = 1:irhs
      inp = sprintf('g_X1(%d)', d);
      for i = 1:irhs
         g_rhs{i} = i_changevar(g_rhs{i}, g_grind.statevars.names{d}, inp);
      end;
      %fprintf(fid, '%s = g_X1(%d);\n',g_grind.statevars.names{d},d);
   end;
end;
nperm = zeros(size(g_grind.permanent));
for i = 1:irhs
   for d = 1:length(g_grind.permanent)
      [g_rhs{i}] = i_changevar(g_rhs{i}, g_grind.permanent{d}.name, sprintf('g_grind.permanent{%d}.currvalue',d));
      nperm(d) = nperm(d) + 1;
   end;
end;
for d = 1:length(g_grind.permanent)
   [funcs] = i_changevar(funcs, g_grind.permanent{d}.name, sprintf('g_grind.permanent{%d}.currvalue',d));
   nperm(d) = nperm(d) + 1;
   if nperm(d) <= 1 %if it occurs only once it needs not initial value (can be more precise)
      g_grind.permanent{d}.initiate = 0;
   end;
end;
if (dimtot == 0)
   if ~isempty(g_grind.externvars);
      fprintf(fid, 'global g_grind;\n');
   end;
   Fun = funcs;
   for d = 1:irhs
      inp = sprintf('g_X1(%d)', d);
      Fun = i_changevar(Fun, g_grind.statevars.names{d}, inp);
   end;
   fwrite(fid, Fun);
else
   fwrite(fid, funcs);
end;
if dimtot > 0
   d = 1;
   for j = 1:length(dim)
      if dim(j) == 1
         fprintf(fid,'g_X2(g_grind.statevars.dims{%d}.from,1) = %s;\n',j, stripsemicolon(stripcomments(char(g_rhs{j}))) );
         d = d + dim(j);
      elseif isempty(dim2)
         fprintf(fid,'g_X2(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to) = %s;\n',j,j, stripsemicolon(stripcomments(char(g_rhs{j}))));
         d = d + dim(j);
      else
         g_rhs{j} = strtrim(g_rhs{j});
         if g_rhs{j}(length(g_rhs{j})) == ';'
            g_rhs{j} = g_rhs{j}(1:length(g_rhs{j} - 1));
         end;
         fprintf(fid,'g_X2(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to) = %s;\n',j,j, stripsemicolon(stripcomments(char(g_rhs{j}))));
         d = d + dim(j) * dim2(j);
      end;
   end;
else
   for i1 = 1:irhs
      fprintf(fid,'g_X2(%d,1) = %s;\n',i1, stripsemicolon(stripcomments(char(g_rhs{i1}))));
   end;
end;
if ~isempty(g_grind.permanent)
   fprintf(fid, 'i_updatepermanent;\n');
end;
%if ~isempty(functions)
%   for i=1:length(functions)
%      if ~isempty(functions{i})
%         fwrite(fid,functions{i});
%      end;
%   end;
%end;
    for i=1:ilocfun
      if ~isempty(localfuncs{i})
         fprintf(fid,'%s\n',localfuncs{i});
      end;
    end;
fclose(fid);
cd(oldpath);
if g_grind.solver.haslag && ~g_grind.solver.isdiffer
   g_grind.solver.name = 'ddesol';
end;
i_modelinit;
if g_grind.solver.isdiffer && g_grind.solver.nonautonomous
   g_grind.solver.nonautonomous = 0;
end;
if isempty(g_grind.statevars.names) && isempty(g_grind.statevars.vectnames)
   error('GRIND:model:NoStatevars','Error: no state variables detected');
else
   if g_grind.statevars.vector
      evalin('base', i_globalstr(g_grind.statevars.vectnames,[],'clear')); %avoid warnings
      evalin('base', i_globalstr(g_grind.statevars.vectnames)); 
   else
      evalin('base', i_globalstr(g_grind.statevars.names,[],'clear')); %avoid warnings
      evalin('base', i_globalstr(g_grind.statevars.names));
   end;
   i_initvar(1); %initiates the state variables
end;
if ~isempty(g_grind.commands) && iscell(g_grind.commands{1})
   g_grind.commands = g_grind.commands{1};
end;
try
   g_grind.starting=1;
   resetpars(1);
    if isfield(g_grind,'starting')
       g_grind=rmfield(g_grind,'starting');
    end;
catch err
    if isfield(g_grind,'starting')
        g_grind=rmfield(g_grind,'starting'); %#ok;
    end;
    rethrow(err);
end;
%randomize; neccesary
%************************************************************************************************************

function [res,difford2] = clearlist(alist, irhs, diffchar)
global g_grind;
difford2 = [];
if ~isempty(alist)
   n = size(alist, 2);
   res = cell(1, n);
   for i = 1:n
      s = char(alist{i});
      for j = 1:irhs
         statevr = [char(g_grind.statevars.names{j}) diffchar]; 
         f = strfind(s, statevr);
         while ~isempty(f)
            if (f(1) == 1) || isempty(strfind( ...
                  'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890', s(f(1) - 1)))
               f2 = strfind(s(f(1):size(s, 2)), ')') + f(1) - 1;
               f = f(1) + strfind(s(f(1):size(s, 2)), '(') - 1; 
               if f(1) + 2 == f2(1)
                  dord2 = 0;
               else
                  dord2 = str2double(s(f(1) + 2:f2(1) - 1));
               end;
               if isempty(difford2)
                  difford2 = dord2;
               elseif dord2 ~= difford2
                  error('GRIND:model:MixofOrders','Error in difference equation: mix of different orders not allowed');
               end;
               s = [s(1:f(1) - 1) s(f2(1) + 1:size(s, 2))];
               if length(statevr) > length(s)
                  f = [];
               else
                  f = strfind(s, statevr);
               end;
            else
               f = f(1) + strfind(s(f(1) + 1:size(s, 2)), statevr);
            end
         end;
      end
      res{i} = char(s);
   end
else
   res = alist;
end;


function s = makearrayop(s1)
s=strrep(s1,'.*','*');
s=strrep(s,'*','.*');
s=strrep(s,'.^','^');
s=strrep(s,'^','.^');
s=strrep(s,'./','/');
s=strrep(s,'/','./');

function s = stripcomments(s)
%global g_grind;
if iscell(s)
    for j=1:length(s)
        i = strfind(s{j}, '%'); %strip off comments
        if ~isempty(i)
            s{j} = s{j}(1:i(1) - 1);
        end;
    end;
else
    i = strfind(s, '%'); %strip off comments
    if ~isempty(i)
        s = s(1:i(1) - 1);
    end;
end;
% if g_grind.version.matlabrelease == 11
%    s=strrep(s,'&&','&');
%    s=strrep(s,'||','|');
% end;
function s = extractlags(s1)
global g_grind;
s = s1;
ilag = strfind(s1, 'lag(');
if ~isempty(ilag)
   if ~isfield(g_grind, 'lags')
      g_grind.lags = {};
   end;
   g_grind.solver.haslag = 1;
   i = 1;
   while ~isempty(ilag);
      j1 = ilag(1) + 4;
      if s(j1)=='''', j1=j1+1; end;
      j2 = j1;
      while s(j2) ~= ','
         j2 = j2 + 1;
      end;
      j2 = j2 - 1;
      if s(j2)=='''',j2=j2-1; end;
      avar = i_getno(s(j1:j2));
      avar = avar.no;
      while s(j2) ~= ','
         j2 = j2 + 1;
      end;
      j1 = j2 + 1;
      while s(j2) ~= ')'
         j2 = j2 + 1;
      end;
      j2 = j2 - 1;
      found = 0;
      for l_k = 1:length(g_grind.lags)
         if strcmp(g_grind.lags, s(j1:j2))
            found = 1;
         end;
      end;
      if ~found
         g_grind.lags = [g_grind.lags; {s(j1:j2)}];
         l_k = length(g_grind.lags);
      end;
      s = sprintf('%sLags(%d,%d)%s',s(1:ilag(i) - 1),avar,l_k,s(j2 + 2:length(s)));
      ilag = strfind(s, 'lag(');
   end;
   g_grind.solver.name = 'ddesol';
end;

function res = sortfuncs(funcs1)
global g_grind;
if isempty(funcs1)
   res = funcs1;
   return;
end;
funcs = removepoints(str2cell(funcs1), 0);
permnames=cell(1,length(g_grind.permanent));
for i = 1:length(g_grind.permanent)
   permnames{i} = g_grind.permanent{i}.name;
end;
if g_grind.statevars.vector
   filter = [g_grind.pars, g_grind.statevars.vectnames,permnames];%,'g_grind','t'};
else
   filter = [g_grind.pars, g_grind.statevars.names,permnames];%,'g_grind','t'};
end;

funnames =  cell(1,length(funcs));
funs = cell(1,length(funcs));
hasif = 0;
for i = 1:length(funcs)
   if strncmp(funcs{i},'if ',3)||strncmp(funcs{i},'end;',4)||strncmp(funcs{i},'for',3)
      hasif = 1;
   end;
   f=strfind(funcs{i}, '=');
   if ~isempty(f)
      s = funcs{i};
      fun.name = strtrim(s(1:f(1) - 1));
      % if ~isempty(setdiff(fun.name,filter)) %not assignment to permanent or statevar
      funnames{i} = fun.name;
      fun.comm = funcs{i};
      s = s(f(1) + 1:end);
      fun.references = i_symvar(s, filter);
      funs{i} = fun;
      %  end;
   else
      warning('GRIND:model:incompleteformula','"%s" incomplete formula?\n', funcs{i});
      fun.references = {};
      funs{i} = fun;
   end;
end;
if hasif
   res = funcs1;
   return;
end;
%remove references to other than functions:
for i = 1:length(funs)
   if ~isempty(funs{i})
      ref = cell(1,length(funs{i}.references));
      l_k = 0;
      for j = 1:length(funs{i}.references)
         found = 0;
         j1 = 1;
         while ~found && (j1<=length(funnames))
            if strcmp(funnames{j1}, funs{i}.references{j})
               found = 1;
            end;
            j1 = j1 + 1;
         end;
         if found
            l_k = l_k + 1;
            ref{l_k} = funs{i}.references{j};
         end
      end;
      funs{i}.references = ref(1:l_k);
   end;
end;

%Sort

sortedfuncs = cell(1,length(funs));
refs = cell(1,length(funs));
j = 0;
changedsome = 1; %safety for errorous funcs;
while (j < length(funs)) && changedsome
   changedsome = 0;
   for i = 1:length(funs)
      if ~isempty(funs{i}.name) && allunion(funs{i}.references, refs)
         changedsome = 1;
         j = j + 1;
         sortedfuncs{j} = funs{i}.comm;
         refs{j} = funs{i}.name;
         funs{i}.name = '';
      end;
   end;
end;
res = sprintf('%s\n', sortedfuncs{:});

function res = allunion(s1, s2)
res = 1;
for i = 1:length(s1)
   found = 0;
   j = 1;
   while ~found && (j<=length(s2))
      if strcmp(s1{i}, s2{j})
         found = 1;
      end;
      j = j + 1;
   end;
   if ~found
      res = 0;
      return;
   end;
end;

function s = stripsemicolon(s)
if ~isempty(s) && (s(end)==';')
   s = s(1:end - 1);
end;

function alist = removepoints(alist, docompound)
h = 1;
i = 1;
siz = length(alist);
s = stripcomments(alist);
sfull = alist;
while h <= length(alist)
   alist{i} = '';
   h2=stackedfindend(s,'[',']',h);
   if isempty(h2) || (h2==h)
      h2=stackedfindend(s,'{','}',h);
      if docompound && (isempty(h2) || (h2==h))
         h2=stackedfindend(s,{'for','if','while','switch','case'},'end',h,1);
      end;
   end;
   %merge h2 lines
   if ~isempty(h2)
      for j = h:h2 - 1
         alist{i} = [char(alist{i}) char(sfull{j}) sprintf('\n') ];
         h = h + 1;
      end;
   end;
   %merge ... lines
   while (h < siz) && (length(strfind(s{h}, '...')) == 1)
      alist{i} = [char(alist{i}) char(sfull{h}) sprintf('\n')];
      h = h + 1;
   end;
   alist{i} = [char(alist{i}) char(sfull{h})];
   h = h + 1;
   i = i + 1;
end;
alist = alist(1: i - 1);

function h = stackedfindend(list, sstart, send, i1,wholewords)
if nargin < 4
   i1 = 1;
end;
if nargin < 5
   wholewords = 0;
end;
h = i1;
stack = 0;
if isempty(findstrings(sstart, list{h}, wholewords))
   h = [];
else
   found = 0;
   while (h < length(list)) && ~found
      stack = stack + length(findstrings(sstart, list{h},wholewords)) - length(findstrings(send, list{h},wholewords));
      found=stack == 0;
      if ~found
         h = h + 1;
      end;
   end;
end;

function f = findstrings(strng, s, wholewords)
if ~iscell(strng)
   f = strfind(s, strng);
   if ~isempty(f) && wholewords
      f = checkwhole(strng, s, f);
   end;
else
   f = [];
   for i = 1:length(strng)
      f1 = strfind(s, strng{i});
      if ~isempty(f1) && wholewords
         f1 = checkwhole(strng{i}, s, f1);
      end;
      f = [f f1];
   end;
end;

function f = checkwhole(substr, s, f1)
len = length(substr);
f = [];
for i = 1:length(f1)
   if ((f1(i)==1) || ~isnumletter(s(f1(i) - 1))) &&  ((f1(i)+len - 1==length(s)) || ~isnumletter(s(f1(i)+len)))
      f = [f f1(i)];
   end;
end;
function res = isnumletter(ch)
res = isletter(ch);
if ~res
   res = ~isempty(strfind('1234567890_', ch));
end;

