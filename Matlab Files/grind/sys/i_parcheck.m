function [Res, S] = i_parcheck(nodlg)
global g_grind;
if nargin == 0
   nodlg = 0;
end;
res = 1;
s = '';
startinggrind=isfield(g_grind,'starting');
LF =  sprintf('\n');
if isfield(g_grind, 'pars') && ~startinggrind
   for i = 1:length(g_grind.pars)
      if (strcmp('eps', char(g_grind.pars{i})))
         i_warningdlg('GRIND:model:eps','"eps" is a reserved word in MATLAB, do not use for parameter name');
      end;
      x = evalin('base', char(g_grind.pars{i}));
      if isempty(x)
         res = 0;
         s = sprintf('%sParameter %s is not initialized\n',s,char(g_grind.pars{i}));
      end
   end
end;
if res
   if isempty(g_grind)
      s = 'No model defined, use <a href="matlab: model">model</a> or <a href="matlab:vismod">vismod</a> to define a model';
      res=0;
   elseif ~isfield(g_grind, 'statevars') || isempty(g_grind.statevars) ||...
         (isempty(g_grind.statevars.names)&&isempty(g_grind.statevars.vectnames))
      s = 'No state variables defined';
      res = 0;
   else
      if ~g_grind.solver.isdiffer
         g_grind.solver.iters = 1;
      end;
      g_grind.solver.opt.OutputFcn = str2func('i_odespeed');
      if ~isempty(g_grind.odefile) && ~startinggrind
         try
            N0 = i_initvar;
         catch %#ok
           if g_grind.statevars.vector
              evalin('base',i_globalstr(g_grind.statevars.vectnames));
           else
              evalin('base',i_globalstr(bg_grind.statevars.names));
           end;
           N0=i_initvar;
        end;
        if ~isempty(g_grind.permanent)
           defpermanent('-deactivate',[]);
        end;
         try
            if ~isempty(N0)
               feval(g_grind.odefile, 1, N0);
            else
               s = 'State variables not initialized';
               res = 0;
            end;
         catch err
 %           err = lasterror;
            s=err.message;
            f=strfind(s,'opentoline(');
            if ~isempty(f)
%opentoline('filename',line,col)
                 f2=strfind(s,')');
                 if ~isempty(f2);
                     comm=s(f(1):f2(1));
                     comm=strrep(comm,'opentoline','getline');
                     f1=strfind(s,'</a>');
                     if ~isempty(f1)
                        s1=s(f1(end)+5:end);
                     else
                        s1='';
                     end;
                     s2=[eval(comm) s1];
                 end
            else
            s2 = s;
            f=strfind(s, '==>');
            if length(f) == 2
               s2 = s(f(2) + 3:length(s));
               s = [s(1:f(1) - 1) LF s2];
               f2 = strfind(s, LF);
               s2 = removegX(s);
               f3 = strfind(s2, LF);
               if f2(1) - f3(1) > 0
                  s2 = [s2(1:f3(1) + 1) s2(f3(1) + 1 + f2(1) - f3(1):length(s2))];
               end;
            end;
            end;
            s = ['SYNTAX error in model:' LF s2];
            res = 0;
         end;
      end;
      for i = 1:length(g_grind.externvars)
         if isempty(g_grind.externvars{i}.data)
            warning('GRIND:externvars:nodata','Some external variables have no data, use <a href="matlab: setdata">setdata</a> to add data');
            break
         end;
      end;
   end;
end;
if ~res
   if ~nodlg
      s = ['****** ERROR in GRIND ******' LF s];
 %     errordlg(s);
      disp(' ');
      error('GRIND:parcheck:ModelNotCorrect',s);
   else
      disp('Model not correct:');
      disp(s);
   end;
end;
if nargout > 0
   Res = res;
   S = s;
end;
function s1=getline(file,line,col)
fid=fopen(file,'r');
for i=1:line
    s1=fgetl(fid);
end;
s1=[s1(1:col) '#' s1(col+1:end)];
s1=removegX(s1);
col=strfind(s1,'#');
s1=strrep(s1,'#','');
s1=sprintf('%s\n%s^\n',s1,repmat(' ',1,col));
fclose(fid);


function s1 = removegX(s)
global g_grind;
if ~isempty(g_grind) && isfield(g_grind, 'statevars') && ~isempty(g_grind.statevars)
   s1 = s;
   for d = 1:length(g_grind.statevars.names)
      inp = sprintf('g_X1(%d)', d);
      s1 = i_changevar(s1, inp, g_grind.statevars.names{d});
      inp = sprintf('g_X2(%d,1)',d);
      s1=i_changevar(s1,inp,[g_grind.statevars.names{d} '''']);
   end;
end;

