function s=i_statevars_names(i)
global g_grind
if g_grind.statevars.vector
   j=1;
   while (j<=length(g_grind.statevars.dims))&&(g_grind.statevars.dims{j}.to<i)
      j=j+1;
   end;
   i=i-g_grind.statevars.dims{j}.from+1;
   [i1,i2]=ind2sub([g_grind.statevars.dims{j}.dim1,g_grind.statevars.dims{j}.dim2],i);
   if g_grind.statevars.dims{j}.dim2>1
      s=sprintf('%s(%d,%d)',g_grind.statevars.vectnames{j},i1,i2);
   else
      s=sprintf('%s(%d)',g_grind.statevars.vectnames{j},i1);
   end;
else
   if length(i)==1
      s=g_grind.statevars.names{i};
   else
      s=sprintf('%s;',g_grind.statevars.names{i});
   end
end;
