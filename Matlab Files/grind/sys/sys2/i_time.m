function [outmat, leg] = i_time(varargin)
global g_Y g_t t g_grind g_data;
nocheck = 0;
OldY = [];
filename = '';
hol = 0;
adding = 0;
silent = 0;
i_parcheck;
%out1 = [];
ndays = [];
if nargin > 0
   i = 1;
   while i <= length(varargin)
     if ~ischar(varargin{i})
       ndays=varargin{i};
    elseif strncmpi(varargin{i}, '-', 1)
         if strncmpi(varargin{i}, '-n', 2)
            nocheck = 1;
         elseif strncmpi(varargin{i}, '-o', 2)
 %           out1 = g_grind.timevars;
            v = {};
            j = i + 1;
            while j <= length(varargin)
               if ~strncmpi(varargin{j}, '-o', 2)
             %     j = length(varargin) + 1;
             %  else
                  ii= str2double(varargin{j});
                  if ~isempty(ii)&&(ii>0)
                     varargin{j}=['-',varargin{j}];
                  end;
                  v = [v varargin(j)];
              %    i = j;
               end;
               j = j + 1;
            end;
            if isempty(v)
               uiwait(i_outdlg);
            else
               out(v{:});
               return;
            end;
         elseif strncmpi(varargin{i}, '-h', 2)
            hol = 1;
            hold on;
         elseif strncmpi(varargin{i}, '-r', 2)
            nocheck = -1;
         elseif strncmpi(varargin{i}, '-f', 2)
            f=strfind(varargin{i}, '=');
            if ~isempty(f)
               filename = varargin{i}(f(1) + 1:end);
               pathname = pwd;
            else
               [filename, pathname] = uiputfile('*.txt', 'Save As');
               if filename==0
                  filename=[];
               end;
            end;
         elseif strncmpi(varargin{i}, '-s', 2)
            silent = 1;
         elseif strncmpi(varargin{i}, '-a', 2)
            adding = 1;
            g_grind.solver.addmode = 1;
         elseif (strncmpi(varargin{i},'-p1',3)||strncmpi(varargin{i},'-pa',3))&&~strcmpi(varargin{i},'-paranal2')
            [OldY, Oldt] = setparanal(1);
            nocheck = 1;
         elseif strncmpi(varargin{i},'-p2',3)||strcmpi(varargin{i},'-paranal2')
            [OldY, Oldt] = setparanal(2);
            nocheck = 1;
         end;
      elseif isempty(ndays)
         ndays = i_checkstr(varargin{i});
      end;
      i = i + 1;
   end;
end;
if isempty(ndays)
   ndays = g_grind.ndays;
end;
N0 = i_initvar;
if (nocheck == -1) || (~nocheck && i_settingschanged(N0, ndays))
   %  disp('running');
   i_ru(g_grind.odefile, t, ndays, N0, 1);
end;
hastimevars = ~isempty(g_grind.timevars);
if hastimevars
   hastimevars = 0;
   for i = 1:length(g_grind.timevars)
      if ~isempty(g_grind.timevars{i})
         hastimevars = 1;
         break;
      end;
   end;
end;
if ~hastimevars
   i_warningdlg('GRIND:time:novars','No time variables defined, use <a href="matlab:out ">out</a> to define output');
end;
if (nargout > 0) || ~isempty(filename)
   outmat1 = cell(1,size(g_grind.timevars, 2));
   leg=cell(1,size(g_grind.timevars, 2));
   for No = 1:size(g_grind.timevars, 2)
      if ~isempty(g_grind.timevars{No})
         leg{No} = g_grind.outt(No);
         outdata = i_getoutfun(g_grind.outt{No});
         outdata2 = [];
         for i = 1:length(g_grind.timevars{No})
            leg{No}{i + 1} = i_disptext(char(g_grind.timevars{No}{i}));
            aa = i_getoutfun(g_grind.timevars{No}{i});
            if ~isempty(g_data) && (size(aa, 1) == size(g_data.t, 1))
               outdata2 = [outdata2, aa];
            else
               outdata = [outdata, aa];
            end;
         end;
         if ~isempty(outdata2)
            size1 = size(outdata);
            size2 = size(outdata2);
            outdata(size1(1) + 1:size1(1) + size2(1), 1:size1(2) + size2(2)) = ones(size2(1), size1(2) + size2(2)) + NaN;
            outdata(1:size1(1) + size2(1), size1(2) + 1:size1(2) + size2(2)) = ones(size1(1) + size2(1), size2(2)) + NaN;
            outdata(size1(1) + 1:size1(1) + size2(1), size1(2) + 1:size1(2) + size2(2)) = outdata2;
            outdata(size1(1) + 1:size1(1) + size2(1), 1) = g_data.t;
         end;
         outmat1{No} = outdata;
      end
   end;

   if ~isempty(filename)
      if ~isempty(pathname)
         if pathname(end)~=filesep
            pathname=[pathname, filesep];
         end;
         filename = [pathname, filename];
      end
      fid = fopen(filename ,'wt');
      for No = 1:size(g_grind.timevars, 2)
         fprintf(fid, '%s\t', leg{No}{:});
         fprintf(fid, '\n');
         for i = 1:size(outmat1{No}, 1)
            s=sprintf('%g\t', outmat1{No}(i,:));
            s=s(1:end-1);
            s=strrep(s,'NaN','');
            fprintf(fid, '%s\n',s);
         end;
         fprintf(fid, '\n');
         fprintf(fid, '\n');
      end;
      fclose(fid);
      fprintf('Data written to %s\n',filename);
   end;   
   if length(outmat1) == 1
      outmat1 = outmat1{1};
   end;
   if nargout>0
      outmat=outmat1;
   end;
end;
if ~silent
   i_timeplot;
end;
if adding
   g_grind.solver.addmode = 0;
end;

if hol
   hold off;
end;

%if ~isempty(out1)
%   g_grind.timevars = out1;
%end;
if ~isempty(OldY)
   g_Y = OldY;
   g_t = Oldt;
end;

function [OldY, Oldt] = setparanal(iparanal)
global g_Y g_t g_paranal;
if iparanal > 0
   OldY = g_Y;
   Oldt = g_t;
   g_Y = g_paranal.Y;
   g_t = g_paranal.t;
   if iparanal == 2
      g_Y = [g_paranal.Yprev; g_Y];
      g_t = [g_t; g_paranal.tprev + g_t(end)];
   end;
else
   OldY = [];
   Oldt = [];
end;
