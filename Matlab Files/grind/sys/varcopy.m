%VARCOPY   Copy variable to clipboard
%   Copy a numeric or a cell array (containing only values or character strings) 
%   or a datset(statistics toolbox') to the Windows clipboard, 
%   as text tab delimited.
%
%   Usage:
%   VARCOPY(g_t) - Copy the variable g_t to the clipboard.
%   VARCOPY g_Y - Copy variable g_Y to clipboard.
%   VARCOPY(ds) - Copy dataset ds to clipboard.
%
%   See also varpaste

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function varcopy(A, precision)
if nargin == 0
   prompt = {'Enter variable:', 'Number of digits to copy:'};
   defaultanswer = {'', '15'};
   answer = inputdlg(prompt, 'varcopy - copy variable to clipboard', 1, defaultanswer);
   A = evalin('base', answer{1});
   precision = str2double(answer{2});
elseif nargin < 2
   precision = 15;
end;
if ~ischar(A) && ~isnumeric(A) && ~islogical(A) &&~isa(A, 'dataset') &&~iscell(A)
   error('GRIND:varcopy:typeNotSupp','Variables of type "%s" are not supported', class(A));
end;
if ischar(A)
   try
      M = evalin('base', A);
   catch %#ok
      M = [];
   end;
   if ~isempty(M)
      A = M;
   end;
end;
if ischar(A)
   clipboard('copy', A);
   return;
end
s = [];
pattrn = sprintf('%%.%dg\\t', precision);
for i = 1:size(A, 1)
   if isnumeric(A) || islogical(A)
      s1 = sprintf(pattrn, A(i, :));
   elseif iscell(A) || isa(A, 'dataset')
      s1 = [];
      for j = 1:size(A, 2)
         if ischar(A{i, j})
            s1 = sprintf('%s%s\t', s1, A{i, j});
         else
            s1 = sprintf(['%s' pattrn], s1, A{i, j});
         end
      end;
   end;
   s = sprintf('%s%s\n', s, s1(1:end - 1));
end;
if isa(A, 'dataset')
   vars = get(A, 'VarNames');
   s1 = sprintf('%s\t', vars{:});
   s = sprintf('%s\n%s', s1(1:end - 1), s);
end;
clipboard('copy', s);



