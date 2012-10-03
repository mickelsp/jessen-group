%WHERE   Show initial conditions
%   display the current initial condition in the current plot.
%
%   Usage:
%   WHERE - shows the initial conditions as triangle (^)
%   WHERE START - shows the initial conditions as triangle (^)
%   WHERE MARK - shows the initial conditions with a marker MARK
%   . = point; o = circle; x = x-mark; + = plus; * = star; s = square;
%   d = diamond; v = triangle (down); ^ = triangle (up); < = triangle (left);
%   > = triangle (right); p = pentagram; h = hexagram
%   WHERE END - shows the final conditions as triangle (v)
%   WHERE END MARK - shows the final conditions with a marker MARK
%   
%   
%   See also null, ru
%
%

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function where(apen, apen1)
global g_Y g_grind t;
wherefinish = 0;
blink = 0;
if (nargin == 0) || (strcmpi(apen, 'start'))
   if nargin == 2
      apen = apen1;
   else
      apen = '^';
   end;
elseif strcmpi(apen, 'sblink')
   blink = 1;
   if nargin == 2
      apen = apen1;
   else
      apen = '^';
   end;
elseif strcmpi(apen, 'end')
   wherefinish = 1;
   if nargin  == 2
      apen = apen1;
   else
      apen = 'v';
   end;
end;
i_parcheck;
oldY = g_Y;
if wherefinish
   N0 = i_initvar;
   if i_settingschanged(N0)
      i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
   end;
   oldY=g_Y;
   g_Y=g_Y(end,:); %#ok
else
   g_Y = transpose(i_initvar); %#ok
end;
oldpen = g_grind.pen;
try
   g_grind.pen.markersize=5;
   if (length(apen)==2)   
       g_grind.pen.pen=apen(2);
       if (apen(1)=='f')
         fcolor=[0 0 0];
       elseif apen(1)=='g'
         fcolor=[0.8 0.8 0.8];
       else
         fcolor=[1 1 1];
       end;
   else
      fcolor='none';
      g_grind.pen.pen = apen;
   end;
   if blink
      g_grind.pen.color = [1 0 0];
      i_phas(gcf, 0);
      ch=get(gca,'children');
      set(ch(1),'markerfacecolor',[1 0 0]);
      pause(0.5);
      set(ch(1),'markerfacecolor',fcolor,'color',[0 0 0]);
  else
      g_grind.pen.color = [0 0 0];
      i_phas(gcf, 0);
      ch=get(gca,'children');
      set(ch(1),'markerfacecolor',fcolor); 
  end;
   g_grind.pen = oldpen;
   g_Y = oldY;
catch err
   g_grind.pen = oldpen;
   g_Y = oldY; %#ok
   rethrow(err);
end;

