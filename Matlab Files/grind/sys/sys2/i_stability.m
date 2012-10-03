function [isstable, issaddle, isspiral] =  i_stability(eigenval, isdiffer)
if isdiffer
   isstable = max(abs(eigenval),[],1) < 1;
   issaddle = (min(abs(eigenval),[],1) < 1) & (max(abs(eigenval),[],1) > 1); %No convert
else
%    if size(eigenval, 1) == 1
%       isstable = real(eigenval) < 0;
%       issaddle = (real(eigenval) > 0) & (real(eigenval) < 0); %No convert
%    else
      isstable =  max(real(eigenval),[],1) < 0;
      issaddle =  (max(real(eigenval),[],1) > 0) & (min(real(eigenval),[],1) < 0); %No convert
%   end;
end;
if size(eigenval, 1) == 1
   isspiral = abs(imag(eigenval)) > 0;
else
   isspiral = max(imag(eigenval)) > 0;
end
