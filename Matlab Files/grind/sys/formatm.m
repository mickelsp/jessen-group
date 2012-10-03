%FORMATM - set the inendation of a m file
%   This function can mess up the inendation, therefore
%    always check the results before approving
function formatm(afile)
if nargin == 0
   afile = uigetfile;
end;
if strcmp(afile, '1')
   d = what(pwd);
   d.m = sort(lower(d.m));
   for i = 1:length(d.m)
      formatm(char(d.m(i)));
   end;
   return;
end;
afile = which(afile);
nincrem = 3;
fid = fopen(afile, 'r');
incr = 0;
line = '';
nmax = 100;
k = 1;
lines = cell(nmax, 1);
iscomment = 0;
while ~feof(fid)
   if ~iscomment && ~isempty(strfind(line, '...')) && (strfind(line, '...')==length(line) - 2)
      tmpinc = nincrem;
   else
      tmpinc = 0;
   end;
   line = strtrim(fgetl(fid));
   %    if ~ischar(line), break, end;
   if ~isempty(line)
      iscomment = (line(1) == '%');
   else
      iscomment = 0;
   end;
   if ~iscomment
      if strcmp(line, 'end') || strcmp(line, 'end;')
         incr = incr  - nincrem;
      end;
      if strcmp(line, 'catch') || strncmp(line, 'catch ', 6) || strcmp(line, 'else') || strncmp(line, 'elseif ', 7) || strncmp(line, 'case ', 5)
         tmpinc = -nincrem;
      end;
      if strncmp(line, 'case ', 5) || strcmp(line, 'otherwise')
         tmpinc = -round(nincrem / 2);
      end;
      if incr < 0
         incr = 0;
      end;
      %repair bug
      line=strrep(line,'./','./');
      line=strrep(line,'.*','.*');
      
      %spacing (not done in lines with >1 strings)
      line=checkspace(line, ' ',{'-'},' ',{'= -','=-','E-','e-','* -','+ -','(-','( -',': -',':-'});
      line=checkspace(line, ' ', {'*'},' ',{'.*'});
      line=checkspace(line, ' ', {'/'},' ',{'./'});
      line=checkspace(line, ' ', {'\'},' ',{'.\'});
      line=checkspace(line, ' ', {'~=', '>=', '<=', '+', '.*', '==', '.\','|', '&', './'}, ' ',{'|','&'}); %No  convert
      line=checkspace(line, ' ', {'&&', '||'}, ' ', {'&','|'}); %No convert
      line=checkspace(line, ' ', {'&', '|'}, ' ', {'&&', '||'}); %No convert
      line=checkspace(line, ' ', {'='},' ',{'~=','>=','<=','=='});
      line=checkspace(line, ' ', {'>'},' ',{'>='});
      line=checkspace(line, ' ', {'<'},' ',{'<='});
      line = checkspace(line, '', {',', ';'}, ' ', {});
      line = checkspace(line, ' ', {'...'}, '', {});
   end;
   lines{k} = [char(ones(1, incr + tmpinc) * ' ') line];
   k = k + 1;
   if k > nmax
      nmax = nmax + 100;
      lines = [lines; cell(100, 1)]; %#ok
   end;
   if ~iscomment && isempty(strfind(line,'; end')) && isempty(strfind( line,', end')) && ...
      (strncmp(line, 'if ', 3) || strncmp(line, 'while ', 6) || ...
      strncmp(line, 'for ', 4)) || strncmp(line, 'switch ', 7)...
      || strncmp(line, 'switch(', 7) || strcmp(line, 'try')
      incr = incr + nincrem;
   end;
end;
fclose(fid);
for i = 1:k - 1
   disp(lines{i});
end;
disp(' ');
disp('Please check the above file');
overwrite = input('Press 1 if you want to overwrite the file > ');
if overwrite
   disp(['overwrite ' afile]);
   LF = sprintf('\n');
   fid = fopen(afile, 'w');
   for i = 1:k - 1
      fwrite(fid,[lines{i}, LF]);
   end;
   fclose(fid);
end;

function line2 = checkspace(line, spbefor, keywords, spafter, exclkw)
line2 = line;
fparen = strfind(line2,'''');
nparen = length(fparen);
if nparen < 3
   for i = 1:length(keywords)
      kw2 = [char(spbefor) keywords{i} char(spafter)];
      fkw = strfind(line2, keywords{i});
      if ~isempty(fkw)
         fkw2 = strfind(line2, kw2);
         k = 1;
         while isempty(fkw2) && (k <= length(exclkw))
            fkw2 = strfind(line2, exclkw{k});
            k = k + 1;
         end;
         if isempty(fkw2)
            if (nparen <= 1)
               line2 = strtrim(strrep(line2, keywords{i}, kw2));
            else
               ok = 1;
               for j = 1:length(fkw)
                  if (fkw(j) > fparen(1)) && (fkw(j) < fparen(2))
                     ok = 0;
                  end;
               end;
               if ok
                  line2 = strtrim(strrep(line2, keywords{i}, kw2));
               end;
            end;
         end;
      end;
   end
end;

