function x=poissonsignal(ts,tperiod,tduration)
x=ts*0;
i=1;
while i<=length(ts);
   t1=ts(i)+drawwaittime(tperiod);
   while (i<=length(ts))&&(ts(i)<t1)
      i=i+1;
   end;
   for k=1:tduration
      if i<=length(ts)
         x(i)=1;
         i=i+1;
      end;
   end;
end;
