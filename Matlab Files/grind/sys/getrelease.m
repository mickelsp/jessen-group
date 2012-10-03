function [h, matlabsys] = getrelease
try
   matlabsys = program_name;
   matlabsys = matlabsys(1:6);
   h = 0;
catch %#ok
   matlabsys = 'matlab';
   v = version;
   v=(v(find(v=='('):find(v==')')));
   if strncmp(v, '(R11',4)
      h = 11;
   elseif strncmp(v, '(R12',4)
      h = 12;
   else
      v = version('-release');
      i = 1;
      while (i < length(v)) && (~isempty(find(v(i)=='1234567890.', 1)))
         i = i + 1;
      end;
      h = str2double(v(1:i - 1));
   end;
end;
