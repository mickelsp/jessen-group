function [varlist, amatrix] = i_loaddata(filename)
fid = fopen(filename, 'r');
if (fid == -1)
   error('GRIND:loaddata:NoFile','File not found or permission denied.');
end
global g_grind;
[g_grind.loaddata.path g_grind.loaddata.name g_grind.loaddata.ext] = fileparts(filename);
TAB = sprintf('\t');
firstline = fgetl(fid);
if strcmp(firstline, '%model') %inifile
    while ~feof(fid) && ~strcmp(firstline, '%[data]')
        firstline = fgetl(fid);
    end;
    if feof(fid)
        varlist = {};
        amatrix = [];
        warning('GRIND:loaddata:EmptyIniFile','The inifile "%s" has no data', filename);
        return;
    end;
    firstline = fgetl(fid);
end;
isTAB = 0;
if ~isempty(strfind(firstline, TAB))
   dlm = TAB;
   isTAB = 1;
elseif ~isempty(strfind(firstline, ','))
   dlm = ',';
elseif ~isempty(strfind(firstline, ' '))
   dlm = ' ';
else
   error('GRIND:loaddata:FileError','File format incorrect: no delimiter found in the file');
end;
if ~isTAB
   firstline = strrep(firstline, dlm, TAB);
end;
firstline= strrep(firstline,' ','_');
f = strfind(firstline, dlm);
if max(isletter(firstline)) %if firstline in A ...Z or a..z
   varlist = firstline;
   amatrix = {};
else
   varlist = sprintf('Column_%d\t',1:length(f)+1);
   amatrix = {firstline};
end;
changeddec = 0;
while 1
   line = fgetl(fid);
   if ~ischar(line)
      break
   else
      if ~isTAB
         line = strrep(line, dlm, TAB);
      elseif strfind(line, ',') %the user made a TAB delimited file with decimal commas
         line=strrep(line,',','.');
         changeddec = 1;
      end;
      f1 = strfind(line, TAB); %fill missing tabs at the end of lines (Excell problem)
      if length(f1)<length(f)
          line=[line ones(1,length(f)-length(f1)).*TAB]; %#ok
      end;
      line = addnans(line);
      amatrix = [amatrix, {line}]; %#ok
   end;
end
if changeddec
   warning('GRIND:loaddata:commas','Decimal commas are changed into dots');
end;
function line = addnans(line)
TAB = sprintf('\t');
i = length(line);
while i > 1
   if (line(i - 1)==TAB) && (line(i)==TAB)
      line = [line(1:i - 1), 'NaN', line(i:end)];
   end;
   i = i - 1;
end;
if line(1) == TAB
   line = ['NaN' line];
end;
if line(length(line)) == TAB
   line = [line 'NaN'];
end;

