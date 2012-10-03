function res=submatrix(A,row,col)
if isempty(col)
   res=A(row,:);
elseif isempty(row)
   res=A(:,col);
else
   res=A(row,col);
end;
