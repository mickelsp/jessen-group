%werkt niet voor matrix notatie
function i_makepotfun(avar,potfun)
global g_grind;
oldpath = pwd;
cd(grindpath);
afile = potfun;
clear(afile); % make sure that a file is released from memory
afile = [afile '.m'];
fin = fopen([g_grind.odefile '.m'], 'r');
fout = fopen(afile, 'w');
line = fgetl(fin);
if ~strcmp(line, '%function created by GRIND');
   error('GRIND:potential:FileError','Can only process ODEFILES created by GRIND');
end;
LF = sprintf('\n');
fwrite(fout, [line LF]);
fwrite(fout, ['function g_X2=curr_pot(' avar ',g_X1);' LF]);
while ~strncmp(line, 'function ', 7)
   line = fgetl(fin);
end;
line = fgetl(fin);
if strncmp(line, 'global ', 7)
   fwrite(fout, [line LF]);
   line=fgetl(fin);
end
if g_grind.statevars.dim > 1
   if g_grind.statevars.vector
      error('GRIND:potentials:Vector','Vector statevars not supported');
   end
   s = i_globalstr(g_grind.statevars.names,avar);
   fwrite(fout, [s LF]);
end;
fwrite(fout, ['t=1;' LF]);
%while isempty(strfind(line, 'Input(')) && ~feof(fin)
%   line = fgetl(fin);
%end;
%if feof(fin)
%   error('Error creating potential function (does not support matrix notation)')
%end;
%while ~isempty(strfind(line, 'Input('))
%   line = fgetl(fin);
%end;
while isempty(strfind(line, 'g_X2')) && ~feof(fin)
   for j=1:g_grind.statevars.dim
      line=strrep(line,['g_X1(' num2str(j) ')'],g_grind.statevars.names{j});
   end;
   fprintf(fout, [line '\n']);
   line = fgetl(fin);
end;
if isempty(strfind(line, 'g_X2('))
   error('GRIND:potential:NoResult','Error creating potential function, cannot find Result');
end;
s=['g_X2(' num2str(i_varno(avar)) ',1)'];
found = 0;
while 1
   if ~isempty(strfind(line,s))
      k = length(s) + 3;
      if ~g_grind.solver.isdiffer
         s1 = ['g_X2=-(',trimsemicolon(strtrim(line(k:length(line)-1))),');\n'];
      else
         s1 = ['g_X2=-(',trimsemicolon(strtrim(line(k:length(line)-1))),'-' avar ');\n'];
      end;
      for j=1:g_grind.statevars.dim
         s1=strrep(s1,['g_X1(' num2str(j) ')'],g_grind.statevars.names{j});
      end;
      fprintf(fout, s1);
      found = 1;
   end;
   if feof(fin)
      break
   end;
   line = fgetl(fin);
end;
fclose(fin);
fclose(fout);
cd(oldpath);
if ~found
   error('GRIND:potential:Curr_pot','Creation of CURR_POT.M not successfull');
end;
function res = trimsemicolon(s)
if ~isempty(s)
   i2 = length(s);
   while (i2 > 1) && (s(i2) == ';')
      i2 = i2 - 1;
   end;
   res=s(1:i2);
else
    res = s;
end;