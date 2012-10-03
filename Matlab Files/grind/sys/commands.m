% COMMANDS - opens a html browser for the online help system
%   COMMANDS M - opens the MATLAB help system

%   Copyright 2011 WUR
%   Revision: 1.0.7 $ $Date: 02-Nov-2011 21:30:08 $
function commands(full)
%global c;
grindver='1.1.8';
s=grindpath;
f=strfind(s,filesep);
s=s(1:f(end)-1);
addpath(s);
if ~exist('i_use','file')
  addpath([grindpath filesep 'sys2']);
end;
if nargin == 1
   f = str2double(full);
   if ~isnan(f)
      full = f;
   elseif ischar(full) && (length(full)==1) && (upper(full(1))=='M')
      helpwin('contents');
      return;
   else
      disp('Opening html browser for the online help system');
      s=['hh "' grindpath '\grindhelp.chm::/' full '.html" &&'];
      dos(s);
      return;
   end;
end;
if (nargin == 0) || ~full
   disp('Opening html browser for the online help system');
   s = grindpath;
   s=['hh "' s '\grindhelp.chm" &&'];
   dos(s);
   %web([s(1:length(s)-4) filesep 'help' filesep 'reference.html']);
elseif (full == 2)
   gcd('grindhelp.html');
   afile=input('file name with help texts (empty=grindhelp.html)> ','s');
   if isempty(afile)
      afile = 'grindhelp.html';
   end;
   fid = fopen(afile, 'r');
   line = '';
   while isempty(strfind(line, '====='))
      line = fgetl(fid);
   end;
   while ~feof(fid)
       updatedat=datestr(now);
       copyright= sprintf('\n%%   Copyright %s WUR\n%%   Revision: %s $ $Date: %s $\n',...
           datestr(now,'yyyy'),grindver,updatedat);
      if ~isempty(strfind(line, '====='))
         f = strfind(line, ' ');
         if isempty(f)
            fclose(fid);
            break;
         end;
         if ~isempty(strfind(line, '(=)'))
            fname = '';
         else
            fname = which([line(f(1) + 1:f(2) - 1) '.m']);
         end;
         if isempty(fname)
            line = fgets(fid);
            while (~feof(fid)) && isempty(strfind(line, '====='))
               line = fgets(fid);
            end;
         else
            disp(['changing: ' fname]);
            fid2 = fopen(fname, 'r');
            line2 = '%';
            while line2(1) == '%'
               line2 = fgets(fid2);
            end;
            if isempty(line2) || (ischar(line2) && (length(strtrim(line2))<=1)) %copyright line
                line2 = fgets(fid2);
                while line2(1) == '%'
                   line2 = fgets(fid2);
                end;
            end;                
            F = fread(fid2);
            fclose(fid2);
            line = fgets(fid);
            fid3 = fopen(fname, 'w');
            while (~feof(fid)) && isempty(strfind(line, '====='))
               line = striphtml(line,1);
               f = strfind(line, '<a') ;
               f1 = strfind(line, '.html">') ;
               try
                   if ~isempty(f1)
               for k=length(f):-1:1
                   if exist(line(f(k)+9:f1(k)+4),'file')
                     line = [line(1:f(k)+8) 'matlab: commands ' line(f(k)+9:f1(k)-1) line(f1(k)+5:end)];
                   else
                     line = [line(1:f(k)+8) 'matlab: help ' line(f(k)+9:f1(k)-1) line(f1(k)+5:end)];
                   end
               end;
                   end;
               catch err
                   fprintf('error reading line: %s\n',line);
                   rethrow(err);
               end;
               if line(1) == ' '
                  fwrite(fid3, ['%' line(2:length(line))]);
               else
                  fwrite(fid3, ['%' line]);
               end;
               line = fgets(fid);
            end;
            fwrite(fid3, copyright);
            fwrite(fid3, line2);
            fwrite(fid3, F);
            fclose(fid3);
         end;
      end;
   end;
   ddir=pwd;
   model('-c');
   cd(grindpath);
   delete('*.ini')
   delete('*.~ini')
   delete('*.asv')
   delete('curr_ode*.m');
   cd(ddir);
   cd ..
   cd setup
   fid=fopen('grindversion.txt', 'w');
   fprintf(fid,'%s\n%s\n',grindver,updatedat);
   fclose(fid);
   s=grindpath;
   f=strfind(s,filesep);
   s=s(1:f(end)-1);
   fprintf('zipping grind to %s%sgrind.zip\n',filesep,s)
   zip('grind.zip',s)
   disp(copyright);
   cd(ddir);
   %   fclose(fid);
   %update date in Contents.m
%    cont = which('Contents.m');
%    fid = fopen(cont, 'r');
%    line1 = fgets(fid);
%    line2 = strtrim(fgets(fid));
%    d = datestr(date);
%    line2=[line2(1:length(line2)-length(d)-2) ' ' d sprintf('\n')];
%    F = fread(fid);
%    fclose(fid);
%    fid2 = fopen(cont, 'w');
%    fwrite(fid2, line1);
%    fwrite(fid2, line2);
%    fwrite(fid2, F);
%    fclose(fid2);
elseif (full == 3) || (full==4)
   %   gcd help;
   afile=input('file name with help texts (empty=grindhelp.html)> ','s');
   LF = sprintf('\n');
   if isempty(afile)
      afile = 'grindhelp.html';
   end;
   %  helppath=grindpath;
   %  helppath=[helppath(1:length(helppath)-4) filesep 'help' filesep];
   %  cd( helppath);
   fid = fopen(afile, 'r');
   if full == 4
      fid2 = fopen('manual.html','w');
      fwrite(fid2, ['<html><head><LINK REL="stylesheet" TYPE="text/css" HREF="style.css">' LF]);
      fwrite(fid2, ['<title>GRIND - manual</title></head>' LF]);
      fwrite(fid2, '<br><br><h1>Description and installing</h1>');
   end;
   line  = '';
   while isempty(strfind(line, '====='))
      line = fgetl(fid);
   end;
   line2='';
   while ~feof(fid)
      if ~isempty(strfind(line, '====='))
         f = strfind(line, ' ');
         if isempty(f)
            break;
         end;
         fname = line(f(1) + 1:f(2) - 1);
         disp(['changing: ' fname]);
         if full == 3
            fid2 = fopen([fname '.html'], 'w');
            fwrite(fid2, ['<html><head><LINK REL="stylesheet" TYPE="text/css" HREF="style.css">' LF]);
            fwrite(fid2, ['<title>' fname '- GRIND command</title></head>' LF]);
         else
            fwrite(fid2, '<br><br><br>');
         end;
         line = fgetl(fid);
         fwrite(fid2,['<h2>' line '</h2>' LF]);
         if full == 3
            fwrite(fid2, '<a href="reference.html">Reference</a> <a href="grind.html">GRIND</a><br><hr width="100%"><br><br>');
            fwrite(fid2, LF);
         end;
         seealso = 0;
         line = strtrim(fgetl(fid));
         while ~feof(fid) && isempty(strfind(line, '====='))
            bold = 0;
            if full == 4
               if  strcmp(line, '<!-- commands -->')
                  fwrite(fid2, '<br><br><h1>Commands</h1>');
               end;
               if strcmp(line, '<!-- vars -->')
                  fwrite(fid2, '<br><br><h1>Global variables</h1>');
               end;
            end;
            if seealso
               f = [-1, strfind(line,', '), length(line) + 1];
               l = [];
               for ii = 2:length(f)
                  name = line(f(ii - 1) + 2:f(ii) - 1);
                  if strcmp(name, 'ODE45')
                     seealso = 0;
                  end;
                  if seealso && (full == 3)
                     l= sprintf('%s<a href="%s.html">%s</a>, ',l,lower(name),name);
                  else
                     l = sprintf('%s%s, ',l ,name);
                  end;
               end;
               line = l(1:length(l) - 2);
            end;
            if strncmp(line,'See also ',9)
                seealso =1;
                bold=1;
                line2=line(10:end);
                line=line(1:8);
            end;
            switch line
             case 'See also:'
               seealso = 1;
               bold = 1;
             case 'Examples:'
               bold = 1;
             case 'Fields:'
               bold = 1;
             case 'Usage:'
               bold = 1;
            end;
            if bold
               fwrite(fid2,['<h3>' line '</h3>' LF]);
            else
               fwrite(fid2, [line LF]);
            end;
            if ~isempty(line2)
                line=strtrim(line2);
                line2='';
            else
               line = strtrim(fgetl(fid));
            end;
         end;
         if full == 3
            fwrite(fid2, ['</body></html>' LF]);
            fclose(fid2);
         end;
      end;
   end;
   if full == 4
      fclose(fid2);
   end;
   fclose(fid);
   !copy grindhelp.hhp grindhelp.~hhp
   fid1=fopen('grindhelp.~hhp','r');
   fid2=fopen('grindhelp.hhp','w');
   s = 'dff';
   while ~strcmpi(s, '[FILES]')
      s = fgetl(fid1);
      fprintf(fid2, '%s\n', s);
   end;
   s = fgetl(fid1);
   while isempty(s) || ((s(1)~='[') && (~feof(fid)))
      s = fgetl(fid1);
   end;
   d = dir('*.html');
   for i = 1:length(d)
      if ~(strcmp(d(i).name,'grindhelp.html') ||strncmp(d(i).name,'~temp',5))
         fprintf(fid2, '%s\n', d(i).name);
      end;
   end;
   fprintf(fid2, '\n%s\n', s);
   while ~feof(fid1)
      s = fgetl(fid1);
      fprintf(fid2, '%s\n', s);
   end;
   fclose(fid1);
   fclose(fid2);
   !grindhelp.bat
else
   d = what(grindpath);
   d.m = sort(lower(d.m));
   for i = 1:length(d.m)
      anam = d.m{i};
      if ~strncmp(anam,'i_',2) && ~strcmp(anam,'commands.m')%&&~strcmp(anam,'trim.m')
         anam = anam(1:strfind(anam, '.m') - 1);
         disp(['================= ' anam ' ====================']);
         h = help(anam);
         c = str2cell(h);
         m = length(c);
         if m > 0
            while (m > 0) && isempty(strtrim(c{m}))
               m = m - 1;
            end;
            for j = 1:m
               disp(c{j})
            end;
         end;
      end
   end;
end;

