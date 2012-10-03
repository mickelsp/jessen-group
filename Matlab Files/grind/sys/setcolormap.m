%SETCOLORMAP   set the colormap to a color gradient
%    Use this function to adapt the colorbar in 3D graphs and to set the default 
%    colormap in GRIND. (It can also be used outside GRIND) You set the bottom 
%    color (ie the color of the lowest value) and the top color and the range 
%    is made based on a gradient in the hue-saturation-value of both colors.
%   
%    Usage: 
%    SETCOLORMAP - sets the color map of the current figure. First you are asked to select
%    the bottom color, thereafter the top color.
%    SETCOLORMAP N - sets the size of the color map to N (default N=64).
%    SETCOLORMAP [R1 G1 B1] [R2 G2 B2] - sets the colormap to a gradient in 2 RGB colors.
%    Res=SETCOLORMAP(N,[R1 G1 B1],[R2 G2 B2]) - saves the colormap (size N) to Res 
%
%    See also viewcells, surf, colormap, colorbar, graph3d

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function result = setcolormap(m, bottomcolor, topcolor)
%hsv = 1;
g = get(0, 'currentfigure');
global g_grind;
if isempty(g)
   if ~isempty(g_grind)
      cm = i_grindcolormap;
   else
      cm = zeros(64, 3);
   end;
else
   cm = get(gcf, 'colormap');
end;
if (nargin < 1)
   m = size(cm, 1);
else
   m = i_checkstr(m);
end;
if (length(m)==1)&&(m<2)&&(nargin<3)
   m=[m,m,m]; %grays assumed;
end;
if length(m) == 3
   if nargin < 2
      topcolor = uisetcolor(cm(m, :), 'Select the top color');
   else
      topcolor = i_checkstr(bottomcolor);
      if (length(topcolor)==1)&&topcolor<2
         topcolor=[topcolor,topcolor,topcolor];
      end;
   end;
   bottomcolor = m;
   m = size(get(gcf, 'colormap'), 1);
else
   if nargin < 2
      bottomcolor = uisetcolor(cm(1, :), 'Select the bottom color');
   else
      bottomcolor = i_checkstr(bottomcolor);
      if (length(bottomcolor)==1)&&bottomcolor<2
         bottomcolor=[bottomcolor,bottomcolor,bottomcolor];
      end;
   end;
   if nargin < 3
      topcolor = uisetcolor(cm(m, :), 'Select the top color');
   else
      topcolor = i_checkstr(topcolor);
      if (length(topcolor)==1)&&topcolor<2
         topcolor=[topcolor,topcolor,topcolor];
      end;
   end;
end;
col = [bottomcolor; topcolor];
if ~isempty(g_grind) && isfield(g_grind, 'pen')
   g_grind.pen.colormap = col;
end;
map = i_grindcolormap(m, col);
if (nargout < 1)&& ~isempty(g)
   colormap(map);
else
   if nargout<1
      disp('No model selected, cannot select colormap for grind');
   end;
   result = map;
end;
