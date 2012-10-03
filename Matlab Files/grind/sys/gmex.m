%GMEX   compile/create mex ODE file
%   Creates or compiles the current ODE file as mex. Use 
%   this command to optimize speed. The generated
%   C file has the same name as the current ini file. Often 
%   the file has to be adjusted manually.
%
%   Usage:
%   GMEX -g - generates a new c file, and overwrites existing!
%   GMEX - compiles the current c file 
%   GMEX -d - deletes the compiled c file. 
%   GMEX -d - deletes the compiled c file. 
%   GMEX -c FILE - compile FILE.c to the current ODE file. 
%
%   See also solver, euler, model, odefile

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function gmex(flag, file)
global g_grind;
if isempty (g_grind)
   error('GRIND:gmex:NoModel','No model selected, make first MATLAB version');
end;
if nargin == 0
   flag = '-c';
   [path, name] = fileparts(g_grind.inifile);
   file = [name, '.c'];
end;
if (nargin == 1)
   if (flag(1) ~= '-')
      file = flag;
      flag = '-c';
   else
      [path, name] = fileparts(g_grind.inifile);
      file = [name, '.c'];
   end;
end;
switch lower(flag)
 case '-c'
   if ~exist(file,'file')
      error('GRIND:gmex:FileUnknown','%s does not exist, perform first "gmex -g"',file);
   end;
   disp('compiling...');
   g = grindpath;
   if getrelease==11
      mex(file,'-O','-DR11',sprintf('-I%s',g),'-output',g_grind.odefile,'-outdir',g);
   else
      mex(file,'-O',sprintf('-I%s',g),'-output',g_grind.odefile,'-outdir',g);
   end;
   fprintf('%s compiled\n', file);
 case '-d'
   ouddir = cd;
   cd(grindpath);
   clear(g_grind.odefile);
   if exist([g_grind.odefile '.dll'],'file')
      delete([g_grind.odefile '.dll']);
   end;
   if exist([g_grind.odefile '.mexw32'],'file')
      delete([g_grind.odefile '.mexw32']);
   end;
   cd(ouddir);
case '-g'
   a=questdlg('Do you want to overwrite existing C file?','Warning: a C file exists','Yes','No','Yes');
   if strcmp(a,'No')
      error('GRIND:gmex:Cancelled','Cancelled');
   end;
   fid = fopen(file, 'w');
   fprintf(fid, '/*=================================================================\n');
   fprintf(fid, '*\n*  GRIND MEX ODE FILE\n');
   fprintf(fid, '*  include "grindmex.h" for all declarations\n');
   fprintf(fid, '*  define here the yprime function\n');
   fprintf(fid, '*\n*=================================================================*/\n');
   fprintf(fid, '#include "grindmex.h"\n\n');
   fprintf(fid,'static void yprime(double g_X2[], double *t, double g_X1[])\n{\n');
   fprintf(fid, '// define state variables\n');
   if ~g_grind.statevars.vector
      for i = 1:g_grind.statevars.dim;
         fprintf(fid, '#define %s g_X1[%d]\n', i_statevars_names(i), i - 1);
      end;
      fprintf(fid, '\n// define primes\n');
      for i = 1:g_grind.statevars.dim;
         fprintf(fid, '#define d%sdt g_X2[%d]\n', i_statevars_names(i), i - 1);
      end;
   else
      for i = 1:length(g_grind.statevars.vectnames);
         fprintf(fid, '#define %s g_X1[%d]\n', g_grind.statevars.vectnames{i}, i - 1);
      end;
      fprintf(fid, '\n// define primes\n');
      for i = 1:length(g_grind.statevars.vectnames);
         fprintf(fid, '#define d%sdt g_X2[%d]\n', g_grind.statevars.vectnames{i}, i - 1);
      end;
   end;
   fprintf(fid, '\n// define parameters \n');
   for i = 1:length(g_grind.pars)
      fprintf(fid, 'double %s = mexGetdouble("%s");\n', g_grind.pars{i}, g_grind.pars{i});
   end;
   
   fprintf(fid, '\n// define local variables\n');
   for i = 1:length(g_grind.funcnames.names)
      fprintf(fid, 'double %s;\n', g_grind.funcnames.names{i});
   end;
   
   fprintf(fid, '\n// differential equations\n');
   for i = 1:length(g_grind.model)
      s = g_grind.model{i};
      if s(1) ~= '%'
         fis=strfind(s, '=');
         if ~isempty(fis)
            facc=strfind(s,'''');
            if ~isempty(facc) && facc(1) < fis
               fprintf(fid, 'd%sdt %s;\n', s(1:facc(1) - 1), s(fis:end));
            else
               fprintf(fid, '%s;\n', s);
            end;
         else
            fprintf(fid, '%s;\n', s);
         end;
      else
         fprintf(fid, '//%s;\n', s);
      end;
   end;
   fprintf(fid, 'return;\n}\n\n');
   fclose(fid);
   fprintf('generated "%s"\n', file);
   disp('check the file and type "gmex" to compile');
   if g_grind.statevars.vector
      i_warningdlg('GRIND:gmex:vector','Vector statevars not fully supported, adjust the file manually');
   end;
   edit(file);
 otherwise
   error('GRIND:gmex:unknownflag','Enknown flag for gmex');
end;

% function [res,difford2] = clearlist(alist, irhs, diffchar,difford)
% 
% difford2 = [];
% if ~isempty(alist)
%    n = size(alist, 2);
%    res = cell(1, n);
%    for i = 1:n
%       s = char(alist{i});
%       for j = 1:irhs
%          statevr = [char(i_statevars_names(j)) diffchar];
%          f = strfind(s,statevr);
%          while ~isempty(f)
%             if (f(1) == 1) || isempty(strfind( ...
%                   'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890',s(f(1) - 1)))
%                f2 = strfind(s(f(1):size(s, 2)), ')') + f(1) - 1;
%                f = f(1) + strfind(s(f(1):size(s, 2)),'(') - 1;
%                if f(1) + 2 == f2(1)
%                   dord2 = 0;
%                else
%                   dord2 = str2num(s(f(1) + 2:f2(1) - 1)); %#ok
%                end;
%                if isempty(difford2)
%                   difford2 = dord2;
%                elseif dord2 ~= difford2
%                   error('Error in difference equation: mix of different orders not allowed');
%                end;
%                s = [s(1:f(1) - 1) s(f2(1) + 1:size(s, 2))];
%                f = strfind(s, statevr);
%             else
%                f = f(1) + strfind(s(f(1) + 1:size(s, 2)), statevr);
%             end
%          end;
%       end
%       res{i} = char(s);
%    end
% else
%    res = alist;
% end;



