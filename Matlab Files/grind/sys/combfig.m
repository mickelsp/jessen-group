%COMBFIG   Combine several figures
%   Combine several existing figures in a nxn matrix of figures. Note: 
%   era (erase) does not clear this figure.
%
%   Usage:
%   COMBFIG - Combine all currently opened figures in one figure.
%   COMBFIG [FIGNO1 FIGNO2.. FIGNOn] - combine the listed figure numbers.
%   COMBFIG [FIGNO1 FIGNO2.. FIGNOn] [ROW COL] - combine the listed 
%   figures, but override the default number of rows/columns (max 4).
%   COMBFIG [FIGNO1 FIGNO2.. FIGNOn] [ROW COL] START - Define a 
%   starting position (default 1). This is usefull for combining several
%   versions of the same figure, for instance to combine 4 versions of
%   Figure 1: COMBFIG 1 [2 2] 1; COMBFIG 1 [2 2] 2;
%   COMBFIG 1 [2 2] 3; COMBFIG 1 [2 2] 4 [..];
%
%   See also copyfig, era

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function combfig(fignrs, Cells, starti)
if nargin  == 0
   fignrs = get(0, 'children');
else
   fignrs = i_checkstr(fignrs);
end;
if nargin  >= 2
   Cells = i_checkstr(Cells);
end;
if nargin  == 3
   starti = i_checkstr(starti);
else
   starti = 1;
end;
nfig = length(fignrs);
if nargin <= 1
   if nfig  == 0
      errordlg('Not enough figures to combine');
      error('GRIND:combfig:TooFewFigs','Not enough figures to combine');
   elseif nfig == 1
      Cells = [1 1];
   elseif nfig == 2
      Cells = [2 1];
   elseif nfig == 3
      Cells = [3 1];
   elseif nfig == 4
      Cells = [2 2];
   elseif nfig <= 6
      Cells = [3 2];
   elseif nfig <= 8
      Cells = [4 2];
   else
      Cells = [3 3];
   end
end;
scale = zeros(1, 2);
for i = 1:2
   if i == 1
      j = 2;
   else
      j = 1;
   end;
   if Cells(i) == 1
      scale(j) = 1;
   elseif Cells(i) == 2;
      scale(j) = 0.45;
   elseif Cells(i) == 3;
      scale(j) = 0.27;
   else
      scale(j) = 0.2;
   end;
end;
h = figure(i_figno('combfig'));
set(h, 'Name', 'Combined figure');
for i = 1:nfig
   s = subplot(Cells(1), Cells(2), i + starti - 1);
   pos = get(s, 'Position');
   delete(s);
   if ~ishandle(fignrs(i))
      error('GRIND:combfig:cannotfindfig','Error combfig: cannot find figure %g', fignrs(i));
   else      
      ax=findall(fignrs(i),'type','axes');
      hs = copyobj(ax, h);
      for j = 1:length(hs)
         pos2 = get(hs(j), 'Position');
         pos2(1:2) = pos(1:2) + scale(1:2) .* (pos2(1:2) - 0.13);
         pos2(3:4) = scale(1:2) .* pos2(3:4);
         set(hs(j), 'Position', pos2);
      end;
   end;
end;



