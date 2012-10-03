function [tr, datar] = i_concatdata(t1, data1, t2, data2)
if isempty(t1)
   tr = t2;
   datar = data2;
   return;
end;
if isempty(t2)
   tr = t1;
   datar = data1;
   return;
end;
ts = sort([t1; t2]);
l1 = size(data1, 2);
l2 = size(data2, 2);
datar = zeros(length(ts), l1 + l2) + NaN;
tr=zeros(length(ts),1)+NaN;
k = 0;
i = 1;
while i <= length(ts);
   tt = ts(i);
   k = k + 1;
   tr(k) = tt;
   i1=find(t1 == tt);
   if ~isempty(i1)
      datar(k, 1:l1) = data1(i1(1), :);
   end;
   i2=find(t2 == tt);
   if ~isempty(i2)
      datar(k, l1 + 1:l1 + l2) = data2(i2(1), :);
   end;
   while (i<=length(ts)) && (tt==ts(i))
      i = i + 1;
   end;
end;
tr = tr(1:k);
datar = datar(1:k, :);
