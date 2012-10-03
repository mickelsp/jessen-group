function anum=i_checkstr(astr)
if isempty(astr)
   anum=[];
   return;
end;
if ischar(astr)
  anum=str2num(astr); %#ok
  if isempty(anum)
     try
       anum=evalin('base',astr);
     catch
       anum=[];
     end;
  end;
else
   anum=astr;
end;



   