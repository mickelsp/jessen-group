function [res2,res,lyaps]=lyapunov1(ndays)
global g_Y g_t g_grind t;
if nargin<1
   ndays=g_grind.ndays;
else
   ndays=i_checkstr(ndays);
end;
N0 = i_initvar;
if i_settingschanged(N0,ndays)
   i_ru(g_grind.odefile, t, ndays, N0, 1);
end;
if ~isempty(g_Y)
   lyaps=g_Y;
   A=eye(g_grind.statevars.dim);
   lyaps2=zeros(length(g_t),size(A,1));
   for i=2:length(g_t)
      N0=g_Y(i,:)';
      [J,lya]=i_eigen(0,g_grind.solver.iters,N0);
      P=A+J*A.*(g_t(i)-g_t(i-1));
      [A,R]=qr(P);
      R=abs(diag(R));
      for j=1:size(A,1)
         if R(j)<1E-40
            lyaps(i,j)=0;
         else
            lyaps2(i,j)=sort(log(R(j)),'descend');
         end;
      end;
      lyaps(i,:)=real(lya)';
   end;
   res=mean(lyaps);
   res2=sum(lyaps2)/(g_t(length(g_t))-t);
end;
