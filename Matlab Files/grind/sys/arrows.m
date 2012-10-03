%ARROWS   Draw arrows in a figure to show directions of lines
%   This command works on all simple line plots. The arrows
%   indicate the direction in which a line is drawn. 
%   NOTE: in 2D phase planes, arrows can be added manually by pointing with the 
%   mouse and pressing Ctrl-Y. Similarly, arrows can be deleted by pressing Shift-Y in 
%   a phase plane.
%
%   Usage:
%   ARROWS - Adds arrows in the current figure (last series that was added, or 
%   a series that is selected)
%   ARROWS delete - delete the arrows.
%   ARROWS nearest [X, Y] find the closest point from the point [X,Y] in a line 
%   and add an arrow.
%   ARROWS delnearest [X, Y] delete the nearest arrow from the point [X,Y]
%
%
%   See also null, phasclick

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:25 $
function arrows(narrow, arg2)
nearest = [];
if nargin < 1
   narrow = 5;
else
   switch narrow
    case 'nearest'
      nearest = i_checkstr(arg2);
      narrow = 1;
    case 'delete'
       ars=findall(gcf,'tag','Arrow');
      delete(ars);
      return;
    case 'delnearest'
      ars=findall(gcf,'tag','Arrow');
      if ~isempty(ars)
         nearest = i_checkstr(arg2);
         mindist = 1E30;
         iarr = 1;
         for i = 1:length(ars)
            ud = get(ars(i), 'userdata');
            dist = min(sqrt((ud.x(1) - nearest(1))^2 + (ud.y(1) - nearest(2))^2),...
            sqrt((ud.x(2) - nearest(1))^2 + (ud.y(2) - nearest(2))^2));
            if dist < mindist
               mindist = dist;
               iarr = i;
            end;
         end;
         delete(ars(iarr));
      end;
      return;
    otherwise
      narrow = i_checkstr(narrow);
   end;
end;
ax = get(0, 'CurrentFigure');
if isempty(ax)
   error('GRIND:arrows:NoFigure','No figure available');
else
   ax = get(ax, 'CurrentAxes');
end;
iser = 1;
series = findobj(ax,'type','line') ;
if ~isempty(series)
   for i = 1:length(series)
      if strcmp(get(series(i),'selected'),'on')
         iser = i;
      end;
   end;
   X = get(series(iser), 'xdata');
   Y = get(series(iser), 'ydata');
   Z = get(series(iser), 'zdata');
   n = length(X);
   if n < narrow
      narrow = n;
   end;
   if isempty(Z)
      if ~isempty(nearest)
         mindist = 1E30;
         iarr = 1;
         for jser = 1:length(series)
            X = get(series(jser), 'xdata');
            Y = get(series(jser), 'ydata');
            for i = 1:length(X)
               dist =  sqrt((X(i) - nearest(1))^2 + (Y(i) - nearest(2))^2);
               if dist < mindist
                  mindist = dist;
                  iarr = i;
                  iser = jser;
               end;
            end;
         end;
         X = get(series(iser), 'xdata');
         Y = get(series(iser), 'ydata');
         col = get(series(iser), 'color');
         if iarr == size(X)
            iarr = iarr - 1;
         end;
         make_arrow([X(iarr), Y(iarr)], [X(iarr+1), Y(iarr + 1)],col);
      else
         col = get(series(iser), 'color');
         for i = 1:narrow
            iarr = floor((i-1) / narrow * (n - 1)) + 1;
            make_arrow([X(iarr), Y(iarr)], [X(iarr+1), Y(iarr + 1)],col);
         end;
      end;
   else
      axis(axis);
      col = get(series(iser), 'color');
      for i = 1:narrow
         iarr = floor((i-1) / narrow * (n - 1)) + 1;
         make_arrow([X(iarr), Y(iarr), Z(iarr)], [X(iarr+1), Y(iarr + 1), Z(iarr+1)],col);
      end;
   end
end
function make_arrow(X, Y, col)

%h = arrow(X,Y,15, 'BaseAngle', 30);
%set(h, 'facecolor', col);
%set(h, 'edgecolor', col);
%
axannotation('arrow',[X(1),Y(1)],[X(2),Y(2)],'Color', col, 'Tag','Arrow', 'HeadLength', 12, 'HeadWidth', 9);
%
%A=arrowline([X(1),Y(1)],[X(2),Y(2)],'arrowsize', 500);
%A=struct(A);
%set(A.arrowhead,'FaceColor', [0 0 0], 'EdgeColor',[0 0 0])
%delete(A.line);
%delete(A.fullline);
