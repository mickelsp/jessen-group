function A1=cell2num(A)
   A1 = zeros(size(A));
   for i = 1:size(A, 2)
      for j = 1:size(A, 1)
         if islogical(A{j,i}) || isnumeric(A{j, i})
            aa=A{j, i};
            if isempty(aa)
               A1(j,i)=NaN;
            else
               A1(j,i)= aa(1);
            end;
         else
            A1(j, i) = str2double(A{j, i});
         end;
      end;
   end;
