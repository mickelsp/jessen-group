% CONNECTMAT - matrix with connections between neighbors
% This function returns a sparse matrix defining the connections between a matrix 
% the matrix gives the connections between the elements of the space matrix
% if A= the space matrix and M the connections then
% d*(M*A(:)) defines the diffusion between 4 neighbors
%
% Usage
% M=CONNECTMAT(SIZX,SIZY,ISROUND) SIZX is the size in X direction, SIZY is the size in the Y direction, ISROUND = 1 for 
% periodic 0 for bouncing boundaries;
function [M]=connectmat(s1,s2,isround)
if nargin<3
   if size(s1)==2
      if nargin==2
         isround=s2;
         s2=s1(2);
         s1=s1(1);
      elseif nargin==1
         s2=s1(2);
         s1=s1(1);
         isround=1;
      else
         error('GRIND:connectmat:UnknownSize','What is the size of the connections matrix?');
      end;       
   else
      if size(s1)==1
         s2=s1;
      end;
      isround=1;
   end;
end;
siz=[s1,s2];
ndxx=zeros(s1*s2*5,1);
ndxy=zeros(s1*s2*5,1);
vals=zeros(s1*s2*5,1);
k=1;
for i=1:s1
   for j=1:s2
      ii=sub2ind(siz,i,j);
      ndxx(k:k+4)=ii;
      ndxy(k)=ii;
      ndxy(k+1)=sub2ind(siz,i,getndx(j-1,s2,isround));
      ndxy(k+2)=sub2ind(siz,i,getndx(j+1,s2,isround));
      ndxy(k+3)=sub2ind(siz,getndx(i+1,s2,isround),j);
      ndxy(k+4)=sub2ind(siz,getndx(i-1,s2,isround),j);
      vals(k)=-4;
      vals(k+1:k+4)=1;
      %Less efficient similar code
      %total diffusion to M
     %(ii,ii)=-4;
      %to neighbors:
     %M(ii,sub2ind(siz,i,getndx(j-1,s2,isround)))=1;
     %M(ii,sub2ind(siz,i,getndx(j+1,s2,isround)))=1;
     %M(ii,sub2ind(siz,getndx(i+1,s1,isround),j))=1;
     %M(ii,sub2ind(siz,getndx(i-1,s1,isround),j))=1;
      k=k+5;
   end;
end;
M=sparse(ndxx,ndxy,vals);
function i=getndx(i,s1,isround)
if isround
  if i<=0
     i=i+s1;
  elseif i>s1
     i=i-s1;
  end;
else
  if i<=0
     i=1;
  elseif i>s1
     i=s1;
  end;
end;

  
