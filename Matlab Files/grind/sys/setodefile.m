%SETODEFILE   Use specific ODEFILE
%   use setodefile to use a specific ODEFILE (is stored in g_grind.odefile). The names of the
%   parameters are extracted from this file into g_grind.pars(should be declared as global variables);
%   This command is only used for special models. Parameters need to be declared separately.
%
%   Usage: SETODEFILE AFILE = AFile is set as the new ODEFILE (see also MATLAB's ODEFILE).
%   SETODEFILE AFILE {'statevar1','statevar2'} = AFile loaded and the list of statevariables is used.
%   SETODEFILE AFILE STATEVARLIST -difference = Treat the file as a difference equation.
%
%   See also model, use, g_grind.odefile, g_grind.pars 

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function setodefile(afile, varlist,opt)
global g_grind;
if isempty(g_grind)
   evalin('base','initgrind');
end;
if nargin==0
   help setodefile;
   error('GRIND:setodefile:NoFile','No file indicated');
end;
if nargin>=3
   if strncmpi(opt,'-d',2)
      g_grind.solver.isdiffer=1;
      g_grind.solver.name='i_differ';
   end;
end;
if nargin>=2
   if ischar(varlist)
      varlist=eval(varlist);
   end; 
   g_grind.statevars.names=varlist;
   g_grind.statevars.dim=length(varlist);
end;
[p,name] = FILEPARTS(afile);
g_grind.odefile = name;
g_grind.pars = cell(1, 50);
name=[name '.m'];
if ~exist(name, 'file')
   cd(findgrindfile(name));
end;
p = 1;
ID = fopen(name,'r');
while 1
   line = fgetl(ID);
   if ~ischar(line), break, end;
   k = strfind(line, 'global ');
   if ~isempty(k)
      isspace = 1;
      for i = 1:k - 1
         if line(i) ~= ' '
            isspace = 0;
         end
      end
      if isspace
         k = k + 7;
         while line(k) ~= ';'
            while line(k) == ' '
               k = k + 1;
            end;
            k0 = k;
            while (line(k)~=' ')&&(line(k)~=';')
               k = k + 1;
            end;
            g_grind.pars{p} = line(k0:k - 1);
            p = p + 1;
         end
      end
   end;
end
g_grind.pars = g_grind.pars(1:p - 1);
fclose(ID);
if g_grind.statevars.dim>0
   if g_grind.statevars.vector
      evalin('base',i_globalstr(g_grind.statevars.vectnames));
   else      
      evalin('base',i_globalstr(g_grind.statevars.names));
   end;   
end
evalin('base',i_globalstr(g_grind.pars));
g_grind.model={'%external odefile'};
i_modelinit;
if g_grind.statevars.dim==0
   disp(' ');
   disp('>>>>Please enter names of state variables as list');
end;
i_parcheck(1);
