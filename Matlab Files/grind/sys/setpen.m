%SETPEN   Set the pen for the plots
%   The user is prompted for the following information:
%   - Color of the pen - select from list
%   - Cycle pen in plots, if 'Yes' is selected each trajectory gets another
%     color
%   - Line style - select from list
%   - Style of marker of data points - select from list
%
%   Usage:
%   SETPEN - User is prompted for information.
%   SETPEN COLOR LINESTYLE MARKERSTYLE CYCLE - if there are arguments, this
%   information is used. (Use the MATLAB abbreviations as specified below).
%   CYCLE can be 1 (=true) or 0 (false).
%
%      color               marker style             line style
%      y     yellow        .     point              -     solid
%      m     magenta       o     circle             :     dotted
%      c     cyan          x     x-mark             -.    dashdot 
%      r     red           +     plus               --    dashed   
%      g     green         *     star
%      b     blue          s     square
%      w     white         d     diamond
%      k     black         v     triangle (down)
%                          ^     triangle (up)
%                          <     triangle (left)
%                          >     triangle (right)
%                          p     pentagram
%
%                          h     hexagram
%
% 
%   See also RU, PHAS, NEXTPEN 

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function setpen(acolor, linestyle, markerstyle,cycle )
global g_grind;
i_parcheck;
penchanged = 0;
if (nargin == 2)
   if ischar(linestyle)
      i = i_checkstr(linestyle);
   end;
   if ~isempty(i)
      cycle = i;
      g_grind.pen.nextpen = i;
      linestyle = g_grind.pen.pen;
      penchanged = 1;
   end;
end;
if (nargin == 3)
   if ischar(markerstyle)
      i = i_checkstr(markerstyle);
   end;
   if ~isempty(i)
      cycle = i;
      g_grind.pen.nextpen = i;
      markerstyle = '';
   end;
end;
colorlist={'blue', 'green', 'red', 'light blue', 'purple', 'ochre', 'black', 'other'};
linestylelist = {'    no line','-     solid',':     dotted','-.    dashdot','--    dashed'};
markstylelist = {'      no symbol', '.     point', 'o     circle', 'x     x-mark',...
   '+     plus', '*     star', 's     square', 'd     diamond', 'v     triangle (down)',...
   '^     triangle (up)', '<     triangle (left)', '>     triangle (right)', 'p     pentagram',...
   'h     hexagram'};
if nargin  == 0
   icolor = listdlg('ListString', colorlist, ...
      'SelectionMode', 'single', ...
      'PromptString', 'Select a color', ...
      'InitialValue', mod(g_grind.pen.i,8) + 1);
else
   acolor = lower(acolor);
   icolor = find(strcmp(acolor, colorlist));
   if isempty(icolor) || (length(icolor) > 1)
      if acolor == 'b'
         icolor = 1;
      elseif acolor == 'g'
         icolor = 2;
      elseif acolor == 'r'
         icolor = 3;
      elseif acolor == 'c'
         icolor = 4;
      elseif acolor == 'm'
         icolor = 5;
      elseif acolor == 'y'
         icolor = 6;
      elseif acolor == 'k'
         icolor = 7;
      else
         icolor = 8;
      end;
   end;
end;
if ~isempty(icolor)
   if icolor ~= 8
      oldnext = g_grind.pen.nextpen;
      g_grind.pen.i = icolor - 1;
      g_grind.pen.nextpen = 1;
      nextpen;
      g_grind.pen.nextpen = oldnext;
   else
      g_grind.pen.i = 8;
      c = uisetcolor;
      if length(c) > 1
         g_grind.pen.color = c;
      end;
   end;
end;
if nargin  == 0
   iline = listdlg('ListString', linestylelist, ...
      'SelectionMode','single',...
      'PromptString','Select a line style',...
      'InitialValue', 2);
   if iline > 0
      penchanged = 1;
      g_grind.pen.pen = strtrim(linestylelist{iline}(1:2));
   end;
elseif nargin  >= 2
   g_grind.pen.pen = linestyle;
   penchanged = 1;
end;
if nargin == 0
   i = listdlg('ListString', markstylelist, ...
      'SelectionMode', 'single', ...
      'PromptString', 'Select a data point symbol', ...
      'InitialValue', 1);
   if i > 0
      markerstyle = markstylelist{i}(1);
   end;
elseif nargin < 3
   i = 0;
else
   i = 1;
end;
if i > 0
   if penchanged
      g_grind.pen.pen = strtrim([g_grind.pen.pen markerstyle]);
   else
      g_grind.pen.pen = strtrim(markerstyle);
   end;
end;
if nargin  == 0
   c=questdlg('Cycle colors in plots?','GRIND');
   if ~strcmp(c, 'Cancel')
      g_grind.pen.nextpen = strcmp(c, 'Yes');
   end;
elseif nargin  >= 4
   g_grind.pen.nextpen = i_checkstr(cycle);
end;
g_grind.pen

