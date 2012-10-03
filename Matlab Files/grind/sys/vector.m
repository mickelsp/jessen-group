%VECTOR   direction vectors
%   Show the direction of change as vectors in the 2D phase plane.
%
%   Usage:
%   VECTOR - shows 20x20 arrows.
%   VECTOR N - shows NxN arrows.
%   VECTOR N -contour - makes a contourplot of the speed of change.
%   VECTOR N -contour -relative - makes a contourplot of the relative speed of change.
%
%   See also phas, null

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function vector(igrid, opt, opt2)
i_parcheck;
dorelative=0;
docontour=0;
if nargin == 0
   docontour=0;
   igrid = 20;
else
   if nargin >= 2
      docontour = strncmpi(opt, '-contour', 2);
      dorelative = strncmpi(opt, '-relative', 2);
   end;
   if nargin==3
      docontour = docontour|strncmpi(opt2, '-contour', 2);
      dorelative = dorelative|strncmpi(opt2, '-relative', 2);
   end;
   if ischar(igrid)&&(igrid(1)=='-')
      docontour = docontour|strncmpi(igrid, '-contour', 2);
      dorelative = dorelative|strncmpi(igrid, '-relative', 2);
      if docontour
          igrid=50;
      else
          igrid=20;
      end;
   else
      igrid = i_checkstr(igrid);
   end;
end;
global g_grind;
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
if ~(iX.isvar||iX.ispar||iX.isext) || ~(iY.isvar||iY.ispar||iX.isext)
   ax('?');
   errordlg('Cannot create vector field if there are no state variables on the axes.');
   error('GRIND:vector:NoStatevars','Cannot create vector field if there are no state variables on the axes');
end
%if ~(isempty(iX.no) || isempty(iY.no))
[X, Y, Vect] = i_vector(igrid, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim, [],dorelative);

[H, new] = i_makefig('phase2');
%  for i=1:size(Vect,3);
%     Vect(:,:,i)=ln(Vect(:,:,i));
%  end;
if new
set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
end;
plotedit('off');
set(H, 'Name', 'Phase plane');
oldhold = ishold;
hold on;
if docontour
   speed = Vect(:, :, 1).^2;
   for i = 2:size(Vect, 3)
      speed = speed + Vect(:, :, i).^2;
   end;
   speed = ln(sqrt(speed));
   contourf(X, Y, speed,20);
   h=gca;
 %  shading flat;
   h1=colorbar;
   ylab = get(h1,'ylabel');
   set(ylab, 'FontAngle',  get(h, 'FontAngle'), ...
          'FontName',   get(h, 'FontName'), ...
          'FontSize',   get(h, 'FontSize'), ...
          'FontWeight', get(h, 'FontWeight'), ...
          'string',     'ln speed of change');
   i_plotdefaults(h1);
 %  axes(h);
else
   if iX.isvar && iY.isvar
      h = myquiver(X, Y, Vect(:, :, iX.no), Vect(:, :, iY.no));
   elseif iY.isvar
      h = myquiver(X, Y, zeros(size(Vect, 1), size(Vect, 2)), Vect(:, :, iY.no));
   elseif iX.isvar
      h = myquiver(X, Y, Vect(:, :, iX.no), zeros(size(Vect, 1), size(Vect, 2)));
      %else everything is zero
   end;
   set(h, 'Color', [0.5 0.5 0.5]);
   h=gca;
end;
set(h, 'XLim', g_grind.xaxis.lim);
set(h, 'YLim', g_grind.yaxis.lim);
if g_grind.statevars.dim > 2
   title(['Valid for ' i_othervars(i_initvar, iX.no, iY.no)]);
end;
if ~oldhold
   hold off;
end;
%plotedit on;
%end;
function  [h] = myquiver(X, Y, Vx, Vy)
maxlen = -1;
minlen = 999999;
for i = 1:size(Vx, 1)
   for j = 1:size(Vx, 2)
      len = norm([Vx(i, j), Vy(i, j)]);
      if len < minlen
         minlen = len;
      end;
      if len > maxlen
         maxlen = len;
      end;
   end;
end;
offset = (maxlen - minlen) / 10;
for i = 1:size(Vx, 1)
   for j = 1:size(Vx, 2)
      len = norm([Vx(i, j), Vy(i, j)]);
      Vx(i, j) = Vx(i, j) * (offset + len) / len;
      Vy(i, j) = Vy(i, j) * (offset + len) / len;
   end;
end;
[h] = quiver(X, Y, Vx, Vy);

