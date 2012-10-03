%NEXTPEN   Set next pen for phase plane
%   Assign the next color for g_grind.pen (the color of the next trajectory in the
%   phase plane).
%
%   See also setpen

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function pen = nextpen(apen)
global g_grind;
if nargin == 0
   if isempty(g_grind)
      apen=[];
   else
      apen = g_grind.pen;
   end;
end;
ncolors = 7;
colororder = [0, 0, 1; ...
   0, 0.5, 0;
1, 0, 0;
0, 0.75, 0.75;
0.75, 0, 0.75;
0.75, 0.75, 0;
0 0 0];
if isempty(apen)
   pen.pen = '-';
   pen.pen2 = '-';
   pen.markersize=5;
   pen.nextpen = 0;
   pen.i = 1;
   pen.color = [0, 0, 1];
   pen.color2 = [0, 0, 1]; % color2 is always cycling also if nextp = n=0
   pen.drawcolor = [0.75, 0, 0];
   %   pen.colormap=[0,0,1;1,0,0];
   pen.colormap = [0 1 1; 1 0 1];
   if ~isempty (g_grind)&&isfield(g_grind,'pen')
      if ~isempty(g_grind.pen)
      pen.linewidth = g_grind.pen.linewidth;
      pen.fontsize = g_grind.pen.fontsize;
      pen.tickdir = g_grind.pen.tickdir;
      pen.box = g_grind.pen.box;
      end;
   else
      pen.linewidth = 1;
      pen.fontsize = 16;
      pen.tickdir = 'out';
      pen.box = 'on';
   end;
else
   pen = apen;
   pen.i = pen.i + 1;
   while pen.i <= 0
      pen.i = pen.i + ncolors;
   end;
   if pen.i > ncolors
      npen=mod(floor(pen.i/ncolors),4);
      switch npen
      case 0
         pen.pen2='-';
      case 1
         pen.pen2=':';
      case 2
         pen.pen2='-.';
      case 3
         pen.pen2='--';
      end;
      if pen.i == 0
         pen.i = ncolors;
      end;
   end;
   ipen=mod(pen.i,ncolors);
   if ipen==0
      ipen=ncolors;
   end;
   pen.color2 = colororder(ipen, :);
   if pen.nextpen
      pen.color = pen.color2;
      pen.pen = pen.pen2;
   end;
end;
if nargout == 0
   g_grind.pen = pen;
end;
