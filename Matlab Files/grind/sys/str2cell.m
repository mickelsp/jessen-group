%STR2CELL - Convert a string with end-of-line marks to a cell structure
%
%    Example:
%    c=str2cell(sprintf('line 1\n line 2\ntest'))
%    c =
%         'line 1'
%         ' line 2'
%         'test'
%    See also SPRINTF
function c = str2cell(s)
if size(s, 1) > 1
   c=cell(size(s,1),1);
   for i=1:size(s,1)
      c{i}=char(s(i,:));
   end;
   c=strtrim(c);
else
   f = strfind(s, sprintf('\n'));
   if isempty(f)
      c{1} = s;
   else
      c = cell(length(f), 1);
      p = 1;
      for i = 1:length(f)
         c{i} = s(p:f(i) - 1);
         p = f(i) + 1;
      end;
      if p < length(s)
         c{length(f) + 1} = s(p:length(s));
      end;
   end;
end;
