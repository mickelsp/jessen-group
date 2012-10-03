%ERA   Erase and close all figures
%   Close all open figures.
%
%   See also e2n, e3r, e2r, e2p

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function era(allfigures)
if nargin == 0
   allfigures = 0;
else
   if strcmpi(allfigures, 'all')
      allfigures = 1;
   end;
end;
if allfigures
   close all hidden;
else
   ch = findobj(get(0,'Children'), 'flat','type','figure');
   max = i_figno('maxno');
   i2 = i_figno('combfig');
   i3=i_figno('dialog');
   for i = 1:length(ch)
      if (ch(i) <= max) && (ch(i) ~= i2)&& (ch(i) ~= i3)
         figure(ch(i));
         close;
      end;
   end;
end;
