function status = i_odephas(t, y, flag)
%ODEPHAS2 2-D phase plane ODE output function.
%   When the string 'odephas2' is passed to an ODE solver as the 'OutputFcn'
%   property, i.e. options = odeset('OutputFcn','odephas2'), the solver
%   calls ODEPHAS2(T,Y) after every timestep.  The ODEPHAS2 routine plots
%   the first two components of the solution it is passed as it is computed,
%   adapting the axis limits of the plot dynamically.  To plot two
%   particular components, specify their indices in the 'OutputSel' property
%   passed to the ODE solver.
%
%   At the start of integration, a solver calls ODEPHAS2(TSPAN,Y0,'init') to
%   initialize the output function.  After each integration step to new time
%   point T with solution vector Y the solver calls STATUS = ODEPHAS2(T,Y).
%   If the solver's 'Refine' property is greater than one (see ODESET), then
%   T is a column vector containing all new output times and Y is an array
%   comprised of corresponding column vectors.  The STATUS return value is 1
%   if the STOP button has been pressed and 0 otherwise.  When the
%   integration is complete, the solver calls ODEPHAS2([],[],'done').
%
%   See also ODEPLOT, ODEPHAS3, ODEPRINT, ODE45, ODE15S, ODESET.

%   Mark W. Reichelt and Lawrence F. Shampine, 3-24-94
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.18 $  $Date: 1997/11/21 23:31:07 $
global g_grind;
status = 0;                             % Assume stop button wasn't pushed.
chunk = 128;                            % Memory is allocated in chunks.
nvar = g_grind.statevars.dim;
f = get(0, 'CurrentFigure');
if nargin < 3 || isempty(flag)           % odephas2(t, y)
   if ~isempty(f)
      ud = get(f, 'UserData');
      if isempty(ud) || (length(fieldnames(ud)) ~= 10)
         ud = [];
         i_odephas([], [], 'done');
      else
         % Append y to ud.y, allocating if necessary.
         nt = length(t);
         chunk = max(chunk, nt);
         rows = size(ud.y, 1);
         oldi = ud.i;
         newi = oldi + nt;
         if newi > rows
            ud.y = [ud.y; zeros(chunk, nvar)];
         end
         ud.y(oldi + 1:newi, :) = y.';
         ud.ncalled = ud.ncalled + 1;
         ud.i = newi;
         if g_grind.truncate
            minx = min(ud.y(1:newi, ud.iX));
            maxx = max(ud.y(1:newi, ud.iX));
            miny = min(ud.y(1:newi, ud.iY));
            maxy = max(ud.y(1:newi, ud.iY));
            if (minx < g_grind.xaxis.lim(1)) || (maxx > g_grind.xaxis.lim(2)) || ...
                  (miny < g_grind.yaxis.lim(1)) || (maxy > g_grind.yaxis.lim(2))
               ud.stop = 1;
            end;
         end;
         set(f, 'UserData', ud);
         if ud.stop == 1                       % Has stop button been pushed?
            status = 1;
            i_odephas([], [], 'done');
         else
            % Rather than redraw all of the data every timestep, we will simply move
            % the line segments for the new data, not erasing.  But if the data has
            % moved out of the axis range, we redraw everything.
            %    xlim = get(gca,'xlim');
            %    ylim = get(gca,'ylim');
            % Replot everything if out of axis range or if just initialized.
            if (oldi < 2)
               %|| ...
               %          (min(y(1,:)) < xlim(1)) || (xlim(2) < max(y(1,:))) || ...
               %          (min(y(2,:)) < ylim(1)) || (ylim(2) < max(y(2,:)))
               if ~ud.p3
                  set(ud.lines, ...
                     'Xdata',ud.y(1:newi,ud.iX), ...
                     'Ydata', ud.y(1:newi, ud.iY));
                  set(ud.line, ...
                     'Xdata',ud.y(oldi:newi,ud.iX), ...
                     'Ydata', ud.y(oldi:newi, ud.iY));
               else
                  set(ud.lines, ...
                     'Xdata',ud.y(1:newi,ud.iX), ...
                     'Ydata', ud.y(1:newi, ud.iY), ...
                     'Zdata', ud.y(1:newi, ud.iZ));
                  set(ud.line, ...
                     'Xdata',ud.y(oldi:newi,ud.iX), ...
                     'Ydata', ud.y(oldi:newi, ud.iY), ...
                     'Zdata', ud.y(oldi:newi, ud.iZ));
               end;
            else
               % Plot only the new data.
               %   set(ud.line,'Color',[1,0,0]...     % "erase" old segment
               %       'Xdata',ud.y(oldi:newi,1), ...
               %       'Ydata',ud.y(oldi:newi,2), ...
               if ~ud.p3
                  set(ud.lines, ...
                     'Xdata',ud.y(oldi:newi,ud.iX), ...
                     'Ydata',ud.y(oldi:newi,ud.iY), ...
                     'Color', g_grind.pen.drawcolor);
                  if newi > 25
                     set(ud.line, ...
                        'Xdata',ud.y(newi - 25:newi,ud.iX), ...
                        'Ydata',ud.y(newi - 25:newi,ud.iY), ...
                        'Color', [0, 0, 0]);
                  else
                     set(ud.line, ...
                        'Xdata',ud.y(1:newi,ud.iX), ...
                        'Ydata',ud.y(1:newi,ud.iY), ...
                        'Color', [0, 0, 0]);
                  end;
               else
                  set(ud.lines, ...
                     'Xdata',ud.y(oldi:newi,ud.iX), ...
                     'Ydata',ud.y(oldi:newi,ud.iY), ...
                     'Zdata', ud.y(oldi:newi, ud.iZ), ...
                     'Color', g_grind.pen.drawcolor);
                  if newi > 25
                     set(ud.line, ...
                        'Xdata',ud.y(newi - 25:newi,ud.iX), ...
                        'Ydata',ud.y(newi - 25:newi,ud.iY), ...
                        'Zdata', ud.y(newi - 25:newi, ud.iZ), ...
                        'Color', [0, 0, 0]);
                  end;
               end
            end;
         end
      end
   end;
else
   switch(flag)
    case 'init'                           % odephas2(tspan,y0,'init')
      ud = get(gcf, 'userdata');
      ud.y = zeros(chunk, nvar);
      ud.i = 1;
      ud.ncalled = 1;
      ud.y(1, :) = y.';
      ud.iX = i_varno(g_grind.xaxis.var);
      ud.iY = i_varno(g_grind.yaxis.var);
      ud.iZ = i_varno(g_grind.zaxis.var);
      % Rather than redraw all data at every timestep, we will simply move
      % the last line segment along, not erasing it.
      if ~isempty(f);
         ax = get(f, 'CurrentAxes');
         
         isnew = length(get(f, 'Children')) < 2;
         if isnew
            set(ax, 'Xlim', g_grind.xaxis.lim);
            set(ax, 'Ylim', g_grind.yaxis.lim);
         end;
         ud.p3=f == i_figno('phase3');
         hold on;
         if ud.p3
            ud.lines = plot3(y(ud.iX),y(ud.iY),y(ud.iZ),'-','Color',g_grind.pen.drawcolor,'EraseMode','none');
            ud.line = plot3(y(ud.iX),y(ud.iY),y(ud.iZ),'-','Color',[0,0,0],'EraseMode','xor');
            set(gca,'DrawMode', 'fast');
            if isnew
               set(ax, 'Zlim', g_grind.zaxis.lim);
               set(ax,'View',[322.5, 30]);
            end;
         else
            ud.lines = plot(y(ud.iX),y(ud.iY),'-','Color',g_grind.pen.drawcolor,'EraseMode','none');
            ud.line = plot(y(ud.iX),y(ud.iY),'-','Color',[0,0,0],'EraseMode','xor');
         end;
         hold off
         set(ax,'DrawMode','fast');
         % The STOP button.
         h = findobj(f,'Tag','stop');
         if isempty(h)
            ud.stop = 0;
            pos = get(0, 'DefaultUicontrolPosition');
            pos(1) = pos(1) - 15;
            pos(2) = pos(2) - 15;
            str = 'ud=get(gcf,''UserData''); ud.stop=1; set(gcf,''UserData'',ud);set(findobj(gcf,''Tag'',''stop''),''Visible'',''off'')';
            uicontrol( ...
               'Style','push', ...
               'String','Stop', ...
               'Position',pos, ...
               'Callback',str, ...
               'Tag','stop');
            set(f, 'DeleteFcn', str);
         else
            set(h,'Visible','on');            % make sure it's visible
            if ishold
               oud = get(f, 'UserData');
               ud.stop = oud.stop;             % don't change old ud.stop status
            else
               ud.stop = 0;
            end
         end
         set(f, 'UserData', ud);
      end;
    case 'done'                           % odephas2([],[],'done')
      if ~isempty(f)
         ud = get(f, 'UserData');
         %  ud.y = ud.y(1:ud.i,:);
         %  set(f,'UserData',ud);
         %  set(ud.lines, ...
         %      'Xdata',ud.y(:,1), ...
         %      'Ydata',ud.y(:,2));
         if ~isempty(ud) && isfield(ud, 'line')
            if ishandle(ud.lines)
               delete(ud.lines);
            end;
            if ishandle(ud.line)
               delete(ud.line);
            end;
            if ~ishold
               set(findobj(f,'Tag','stop'),'Visible','off');
               refresh;                          % redraw figure to remove marker frags
            end
         end;
      end;
      %
   end
end
if exist('ud','var')&&~isempty(ud)
   if (g_grind.slowdown < 0.01)
      if imod(ud.ncalled, 5) == 0
         pause(g_grind.slowdown)
      else
         drawnow;
      end;
   elseif (g_grind.slowdown < 0.02)
      if imod(ud.ncalled, 5) < 2
         pause(g_grind.slowdown / 2)
      else
         drawnow;
      end;
   else
      pause(g_grind.slowdown / 10);
   end;
   if isfield(ud, 'savingAVI') && ud.savingAVI
      if imod(ud.ncalled, 10) == 0
         i_makeavi('getframe');
      end;
   end;
end;
function m = imod(x, y)
%fast simple mod (for integers)
m = x - y .* floor(x ./ y);
