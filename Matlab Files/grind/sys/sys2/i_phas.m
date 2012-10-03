% i_phas initial function to create phase plane
function [] = i_phas(curplot, newrun)
global g_grind g_t g_Y;
iX = i_varno(g_grind.xaxis.var);
iY = i_varno(g_grind.yaxis.var);
i = i_figno('phase1');
done = 0;
if curplot == i
   done = 1;
   [H, new] = i_makefig('phase1');
   if new
      set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
      set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
   end;
   oldhold = ishold;
   hold on;
   if g_grind.solver.isdiffer
      if new
         close(H);
         itermap
         hold on;
      end;
      set(H, 'Name', 'Iteration map');
      %difference equation: cobwebbing
      if size(g_Y, 1) == 1
         plot(g_Y, g_Y, g_grind.pen.pen, 'Color', g_grind.pen.color, 'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize) %plot equilibrium only
      else
         %cobwebbing plotting method for difference equations
         if ~isempty(iX)
            cobwebbing(g_Y(:, iX), g_grind.drawnow);
         else
            cobwebbing(i_getoutfun(g_grind.xaxis.var), g_grind.drawnow);
         end;
      end;
   else
      %1D differential equation: draw on x axis
      if new
         close(H);
         plotdiff;
         hold on;
      end;
      set(H, 'Name','Plot of 1D differential equation');
      if ~isempty(iX)
         if g_grind.drawnow
            i_drawslowly(g_Y(:, iX), zeros(1, size(g_Y, 1)), 1);
         end;
         plot(g_Y(:, iX), zeros(1, size(g_Y, 1)), g_grind.pen.pen, 'Color', g_grind.pen.color2, 'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
      else
         YY = i_getoutfun(g_grind.xaxis.var);
         if g_grind.drawnow
            i_drawslowly(YY, zeros(1, size(YY, 1)), 1);
         end;
         set(H, 'Name','Plot of 1D differential equation')
         plot(YY, zeros(1, size(YY, 1)), g_grind.pen.pen, 'Color', g_grind.pen.color2,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
      end
   end;
   xlabel(i_disptext(g_grind.xaxis.var));
   if ~oldhold
      hold off;
   end;
elseif newrun
   if ishandle(i)
      set(i, 'visible', 'off');
   end;
end;
i = i_figno('phase2');
if curplot == i
   done = 1;
   [H, new] = i_makefig('phase2');
   if new
      set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
      set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
   end;
   oldhold = ishold;
   hold on;
   if isempty(g_grind.yaxis.var)
      close(H);
      i_phas(i_figno('phase1'), newrun);
      return;
   elseif ~(isempty(iX) || isempty(iY))
      plot(g_Y(:, iX), g_Y(:, iY), g_grind.pen.pen, 'Color', g_grind.pen.color, 'linewidth',g_grind.pen.linewidth,'Erase', 'none', 'markersize',g_grind.pen.markersize);
   elseif ~isempty(iX)
      plot(g_Y(:, iX), i_getoutfun(g_grind.yaxis.var), g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
   elseif ~isempty(iY)
      plot(i_getoutfun(g_grind.xaxis.var), g_Y(:, iY), 'Color', g_grind.pen.color2,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
   else
      plot(i_getoutfun(g_grind.xaxis.var), i_getoutfun(g_grind.yaxis.var), 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize)
   end;
   set(H, 'Name', 'Phase plane');
   xlabel(i_disptext(g_grind.xaxis.var));
   ylabel(i_disptext(g_grind.yaxis.var));
   if ~oldhold
      hold off;
   end;
elseif newrun
   if ishandle(i)
      set(i, 'visible', 'off');
   end;
end;
if ~isempty(g_grind.zaxis.var)
   i = i_figno('phase3');
   if curplot == i
      done = 1;
      [H, fignew] = i_makefig('phase3');
      if fignew
        set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
        set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
      end;
      set(H, 'Name', '3D phase space');
      oldhold = ishold;
      hold on;
      %plotedit off;
      iX = i_varno(g_grind.xaxis.var);
      iY = i_varno(g_grind.yaxis.var);
      iZ = i_varno(g_grind.zaxis.var);
      if ~(isempty(iX) || isempty(iY) || isempty(iZ))
         plot3(g_Y(:, iX), g_Y(:, iY), g_Y(:, iZ), g_grind.pen.pen, 'Color', g_grind.pen.color, 'Erase', 'none', 'markersize',g_grind.pen.markersize);
      else
         plot3(i_getoutfun(g_grind.xaxis.var), i_getoutfun(g_grind.yaxis.var), i_getoutfun(g_grind.zaxis.var) ...
            , g_grind.pen.pen, 'Color', g_grind.pen.color, 'Erase', 'none', 'markersize',g_grind.pen.markersize)
      end;
      set(gca,'DrawMode', 'fast');
      box on;
      xlabel(i_disptext(char(g_grind.xaxis.var)));
      ylabel(i_disptext(char(g_grind.yaxis.var)));
      zlabel(i_disptext(char(g_grind.zaxis.var)));
      if ~oldhold
         hold off;
      end;
      if fignew
         set(gca, 'View', [322.5, 30]);
      end;
   elseif newrun
      if ishandle(i)
         set(i, 'visible', 'off');
      end;
   end;
end;

i = i_figno('poinmap');
if curplot == i
   figure(i);
   us = get(i, 'userdata');
   if ~isempty(us)
      poincaremap(us.avar1, us.avar, us.avalue, us.increasing);
   end;
end;
i = i_figno('poinsec');
if curplot == i
   figure(i);
   us = get(i, 'userdata');
   if ~isempty(us)
      poincaresect(us.avar, us.avalue, us.increasing);
   end;
end;
i = i_figno('lorenzmap');
if curplot == i
   figure(i);
   us = get(i, 'userdata');
   if ~isempty(us)
      lorenzmap(us);
   end;
end;
i = i_figno('torus');
if curplot == i
   done = 1;
   figure(i);
   us = get(i, 'userdata');
   if ~isempty(us)
      torus(us.period, us.increment);
   else
      torus
   end;
end;

i = i_figno('dirfield');
if curplot == i
   done = 1;
   figure(i);
   hold on;
   plot(g_t, g_Y(:, iX), g_grind.pen.pen, 'Color', g_grind.pen.color, 'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
   hold off;
end;
i = i_figno('potent3');
if curplot == i
   done = 1;
   figure(i);
   hold on;
   plot3(g_Y(:, iX), g_Y(:, iY), zeros(size(g_Y, 1)), g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
   set(gca,'DrawMode', 'fast');
   hold off;
end;
i = i_figno('potent2');
if curplot == i
   done = 1;
   figure(i);
   hold on;
   plot(g_Y(:, iX), g_Y(:, iY), g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
   hold off;
end;
i = i_figno('potent1');
if curplot == i
   done = 1;
   figure(i);
   hold on;
   plot(g_Y(:, iX), 0, g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
   plot(g_Y(size(g_Y,1), iX), 0, 'o', 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
   hold off;
end;
i = i_figno('growths');
if ishandle(i + 1)
   done = 1;
   for j = 1:g_grind.statevars.dim
      if curplot == i + j
         figure(i + j);
         hold on;
         plot3(g_Y(:, iX), g_Y(:, iY), zeros(size(g_Y, 1)), g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
         set(gca,'DrawMode', 'fast');
         hold off;
      end;
   end;
end;
i = i_figno('varcontour');
if ishandle(i + 1)
   done = 1;
   varcontour;
end;

if ~done
   if isempty(g_grind.zaxis.var)
      i_phas(i_figno('phase2'), newrun)
   elseif isempty(g_grind.yaxis.var)
      i_phas(i_figno('phase1'), newrun)
   else
      i_phas(i_figno('phase3'), newrun)
   end;
end;
if g_grind.pen.nextpen
   nextpen;
end;


function cobwebbing(g_Y, drawnow)
global g_grind;
cobweb = zeros(size(g_Y, 1) * 2 - 1, 2);
cobweb(1, 1) = g_Y(1);
cobweb(1, 2) = 0;
for i = 2:size(g_Y, 1)
   cobweb((i - 1) * 2, 1) = g_Y(i - 1);
   cobweb((i - 1) * 2, 2) = g_Y(i);
   cobweb((i - 1) * 2 + 1, 1) = g_Y(i);
   cobweb((i - 1) * 2 + 1, 2) = g_Y(i);
end
if drawnow
   i_drawslowly(cobweb(:, 1), cobweb(:, 2));
end;
plot(cobweb(:, 1), cobweb(:, 2), g_grind.pen.pen, 'Color', g_grind.pen.color2,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);

