% VARPASTE- Paste any matrix from clipboard
%   Paste any matrix from the Windows clipboard, which was stored 
%   as text tab delimited (for instance with Microsoft Excel). If there are text
%   strings in the clipboard a cell array is created else a matrix.
%
%   Usage:
%   A=VARPASTE - Paste the contents of the clipboard to the variable A (cell array or matrix).
%   A=VARPASTE(TYPE) - Paste the contents of the clipboard to the variable A, force A 
%   to be of a type TYPE. TYPE can be 'numeric','string','cell','logical' or 'dataset' (see
%   Statistics Toolbox for the latter).
%
%   See also varcopy

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function A = varpaste(forcetype)
if nargin == 0
   forcetype = '';
end;
lines = str2cell(clipboard('paste'));
if strncmpi(forcetype, 's', 1) %string
   A = lines;
   return;
end;
line = lines{1};
if (line(1)=='[') || (line(1)=='{')
   line=sprintf('%s',lines{:});
   try
     A = eval(line);
   catch %#ok
     A=stabread(lines);
   end;
else
   A = stabread(lines);
end
if strncmpi(forcetype, 'n', 1) && iscell(A) %numeric
   A = cell2num(A);
end;
if strncmpi(forcetype, 'c', 1) && isnumeric(A) %cell
   A = num2cell(A);
end;
if strncmpi(forcetype, 'l', 1) %logical
   M = false(size(A));
   hasnan = false;
   if iscell(A)
      for i = 1:size(A, 1)
         for j = 1:size(A, 2)
            if (isnumeric(A{i, j}) || islogical(A{i, j}))
               if isnan(A{i, j})
                  hasnan = true;
               else
                  M(i, j) = logical(A{i, j});
               end
            elseif any(strcmpi(A{i,j},{'on','yes'}))
               M(i, j) = true;
            end;
         end;
      end;
   end;
   if isnumeric(A)
      nans = isnan(A);
      if any(nans)
         M(nans) = 0;
         M(~nans) = logical(A(~nans));
         hasnan = 1;
      else
         M = logical(A);
      end;
   end;
   if hasnan
      warning('GRIND:varpaste:logicalNaN','Logical NaN are not defined: NaNs converted to false');
   end;
   A = M;
end;
if strncmpi(forcetype, 'd', 1) && iscell(A) %dataset
   ds = dataset;
   for i = 1:size(A, 2)
      AA = cell2num(A(2:end, i));
      if all(isnan(AA))
         AA = A(2:end, i);
      end;
      ds.(A{1, i}) = AA;
   end;
   A = ds;
end;



function A = stabread(lines)
hasstr = 0;
for j = 1:length(lines)
   line = lines{j};
   if ~isempty(line)
      [aa, n, err] = sscanf(line, '%f');
      if isempty(err)
         if isempty(aa)
             A=lines;  
             hasstr=1;
         else
         for i = 1:length(aa)
            A{j, i} = aa(i);
         end;
         end;
      else
         f = [0,strfind(line, sprintf('\t')),length(line) + 1];
         for i = 1:length(f) - 1
            s = line(f(i) + 1:f(i + 1) - 1);
            J = str2num(s); %#ok
            if  (~isempty(J) && ~isnan(J))
               res = J;
            else
               res = s;
               hasstr = 1;
            end;
            A{j, i} = res;
         end;
      end;
   end;
end;
if ~hasstr
   A = cell2num(A);
end;


