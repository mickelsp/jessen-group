function i_movie(flag, par1, par2, par3, par4)
global g_Y g_grind g_t t;
h0 = gcbf;
if isempty(h0)
   h0 = get(0, 'currentfigure');
end;
switch flag
 case 'movplot'
   ud.varno = par1;
   t_pause = 0.01;
   tframe = [];
   makeavi = [];
   if isempty(h0)
       H = i_makefig('varcontour', 1);
   else 
       H=h0;
   end;
   set(H,'doublebuffer','on')
   set(H, 'Name', 'Movie plot');
   %  i = 1;
   if isempty(ud.varno.Zs) && isempty(ud.varno.Ys);
      hcolor = feval(ud.varno.type, ud.varno.Xs(end, :));
   elseif isempty(ud.varno.Zs)
      hcolor = feval(ud.varno.type, ud.varno.Xs(end, :), ud.varno.Ys(end, :));
   else
      hcolor = feval(ud.varno.type, ud.varno.Xs(end, :), ud.varno.Ys(end, :), ud.varno.Zs(end, :));
   end;
   if isfield(ud.varno,'xlabel')
       xlabel(ud.varno.xlabel);
   end;
   if isfield(ud.varno,'ylabel')
       ylabel(ud.varno.ylabel);
   end;
   if isfield(ud.varno,'zlabel')
       zlabel(ud.varno.zlabel);
   end;
   set(gca,'drawmode','fast');
   set(hcolor,'erasemode','normal');
   i_plotdefaults;
   i_movie('init', ud.varno, t_pause, hcolor, 'movplot');
   hb=findobj(H,'Tag','TimeSlider');
   set(hb,'Value',1);
   if ~isempty(ud.varno.Ys)
       ylim([0, max(max(ud.varno.Ys))]);
   end;
   if isfield(ud.varno,'Zs')&&~isempty(ud.varno.Zs)
      zlim([0, 1.1*max(max(ud.varno.Zs))]);
   end;
   set(hb,'min',1, 'max', size(ud.varno.Xs, 1));
   if ~isempty(tframe)
      i_movie('sett', tframe);
   elseif ~isempty(makeavi)
      i_movie('writeavi');
   else
      i_movie('activate', H);
   end;
 case 'lines'
   statevars = g_grind.statevars.vectnames;
   t_pause = par1;
   tframe = par2;
   i = par3;
   makeavi = [];
   H = zeros(length(statevars), 1);
   %   for i = 1:length(statevars)
   H(i) = i_makefig('varcontour', i);
   set(H(i),'doublebuffer','on')
   set(H(i), 'Name', ['Vector plot of ', statevars{i}]);
   [dim1] = i_getvardims(statevars{i});
   [sfrom, sto] = i_statevarnos(statevars{i});
   if sto > sfrom
      A1 = reshape(g_Y(1,sfrom:sfrom + dim1 - 1)', 1,dim1);
      xbins = getxbins(i);
      if isempty(xbins)
         hcolor = plot(A1);
      else
         hcolor = plot(xbins, A1);
      end;
      maxY = max(max(g_Y(:, sfrom:sfrom + dim1 - 1)));
      xlim([1 sto - sfrom])
      ylim([0 maxY]);
      xlabel('Column');
      %     set(gca,'CLimMode','manual','CLim',[minY maxY]);
      set(gca,'drawmode','fast','ydir','reverse');
      set(hcolor,'erasemode','normal');
      %   shading flat;
      i_plotdefaults;
      i_movie('init', i, t_pause, hcolor, 'lines');
      if ~isempty(tframe)
         i_movie('sett', tframe);
      elseif ~isempty(makeavi)
         i_movie('writeavi');
      else
         i_movie('activate', H(i));
      end;
   end;
   %  end case 'bars'
 case 'bars'
   statevars = g_grind.statevars.vectnames;
   t_pause = par1;
   tframe = par2;
   i = par3;
   makeavi = [];
   H = zeros(length(statevars), 1);
   %   for i = 1:length(statevars)
   H(i) = i_makefig('varcontour', i);
   set(H(i),'doublebuffer','on')
   set(H(i), 'Name', ['Histogram of ', statevars{i}]);
   [dim1] = i_getvardims(statevars{i});
   [sfrom, sto] = i_statevarnos(statevars{i});
   if sto > sfrom
      A1 = reshape(g_Y(1,sfrom:sfrom + dim1 - 1)', 1,dim1);
      xbins = getxbins(i);
      if isempty(xbins)
         hcolor = bar(A1);
      else
         hcolor = bar(xbins, A1);
      end;
      maxY = max(max(g_Y(:, sfrom:sfrom + dim1 - 1)));
      xlim([1 sto - sfrom])
      ylim([0 maxY]);
      xlabel('Column');
      %     set(gca,'CLimMode','manual','CLim',[minY maxY]);
      set(gca,'drawmode','fast','ydir','reverse');
      set(hcolor,'erasemode','normal','CDataMapping','scaled','BackFaceLighting','unlit');
      %   shading flat;
      i_plotdefaults;
      i_movie('init', i, t_pause, hcolor, 'bars');
      if ~isempty(tframe)
         i_movie('sett', tframe);
      elseif ~isempty(makeavi)
         i_movie('writeavi');
      else
         i_movie('activate', H(i));
      end;
   end;
   %  end
   
 case 'varcontour'
   statevars = g_grind.viewcells.vars;
   t_pause = par1;
   tframe = par2;
   colmap = par3;
   makeavi = par4;
   H = zeros(length(statevars), 1);
   for i = 1:length(statevars)
      H(i) = i_makefig('varcontour', i);
      set(H(i),'doublebuffer','on')
      set(H(i), 'Name', ['Grid of ', statevars{i}.name]);
      g_grind.viewcells.y{i}=i_getoutfun(statevars{i}.name);
      maxY = max(max(g_grind.viewcells.y{i}));
      minY =  min(min(g_grind.viewcells.y{i}));
      A1 = reshape(g_grind.viewcells.y{i}(i,:), statevars{i}.dims(1),statevars{i}.dims(2));
      if statevars{i}.dims(2) == 1
         A1 = repmat(A1, 1, 3);
      end;
      A1(1) = maxY;
      A1(2) = minY;
      hcolor = pcolor(A1);
      daspect([1 1 1]);
      xlabel('Column');
      ylabel('Row');
      set(gca,'CLimMode','manual','CLim',[minY maxY]);
      set(gca,'drawmode','fast','ydir','reverse');
      set(hcolor,'erasemode','normal','CDataMapping','scaled','BackFaceLighting','unlit');
      shading flat;
      colormap(colmap);
      h = colorbar;
      if minY == maxY
         maxY = minY + 1;
      end;
      set(h, 'ylim', [minY, maxY]);
      i_plotdefaults;
      i_movie('init', i, t_pause, hcolor, 'varcontour');
      if ~isempty(tframe)
         i_movie('sett', tframe);
      elseif ~isempty(makeavi)
         i_movie('writeavi', makeavi);
      else
         i_movie('activate', H(i));
      end;
   end
 case 'init'
   if isempty(findobj(h0,'Tag','TimeSlider'))
      l = length(g_t);
      if l == 0
         l = 100;
      end;
      uicontrol('Parent',h0, ...
         	'Units','points', ...
         	'ListboxTop',0, ...
         	'Callback','i_movie(''sliderset'')', ...
         	'Position',[40 1 230 10], ...
         'max', l, ...
         'min', 1, ...
         'value', 1, ...
         	'Style','slider', ...
         'Tag','TimeSlider');
      uicontrol('Parent',h0, ...
         	'Units','points', ...
         	'Callback','i_movie(''play'',-1)', ...
         	'ListboxTop',0, ...
         	'Position',[271 1 10 10], ...
         'String','<<', ...
         'Tooltipstring','Play backwards',...
         'Tag','hb2');
      uicontrol('Parent',h0, ...
         	'Units','points', ...
         	'Callback','i_movie(''step'',-1)', ...
         	'ListboxTop',0, ...
         	'Position',[282 1 10 10], ...
         'String','<', ...
         'Tooltipstring','One step backward',...
         	'Tag','hb3');
      uicontrol('Parent',h0, ...
         	'Units','points', ...
         	'Callback','i_movie(''deactivate'')', ...
         	'ListboxTop',0, ...
         	'Position',[293 1 20 10], ...
         	'String','<>', ...
         'Tooltipstring','Stop',...
         	'Tag','hb4');
      uicontrol('Parent',h0, ...
         	'Units','points', ...
         	'Callback','i_movie(''step'',1)', ...
         	'ListboxTop',0, ...
         	'Position',[314 1 10 10], ...
         	'String','>', ...
         'Tooltipstring','One step forward',...
         'Tag','hb5');
      uicontrol('Parent',h0, ...
         	'Units','points', ...
         	'Callback','i_movie(''play'',1)', ...
         	'ListboxTop',0, ...
         	'Position',[325 1 10 10], ...
         	'String','>>', ...
         'Tooltipstring','Play',...
         'Tag','hb6');
      uicontrol('Parent',h0, ...
         	'Units','points', ...
         	'Callback','i_movie(''hide'',''off'')', ...
         	'ListboxTop',0, ...
         'Position',[336 1 10 10], ...
         'String','x', ...
         'Tooltipstring','Cancel movie',...
         'Tag','hb7');
      i_makeavi('addmenu', 13)
   else
      hb=findobj(h0,'Tag','TimeSlider');
      set(hb,'max',length(g_t), 'min', 1);
   end;
   ud.active = 1;
   ud.step = 1;
   ud.savingAVI = 0;
   ud.varno = par1;
   ud.pause = par2;
   ud.handle = par3;
   ud.type = par4;
   set(gcf, 'userdata',ud);
 case 'sliderset'
   ud = get(h0, 'userdata');
   hb=findobj(h0,'tag','TimeSlider');
   t1 = round(get(hb, 'value'));
   maxt1 = get(hb, 'max');
   mint1 = get(hb, 'min');
   if t1 >= maxt1
      t1 = maxt1;
      set(hb, 'value', t1);
   end;
   if t1 < mint1
      t1 = mint1;
      set(hb, 'value', t1);
   end;
   if ~ud.active
      i_movie('activate');
   end;
 case 'rewind'
   hb=findobj(h0,'tag','TimeSlider');
   set(hb, 'value', 1);
   ud = get(h0, 'userdata');
   if ~ud.active
      i_movie('activate');
   end;
 case 'stopavi'
   i_movie('deactivate');
   i_makeavi('stop/save');
 case 'writeavi'
   if nargin < 2
      par1 = [];
   end;
   i_makeavi('init', par1);
   ud = get(h0, 'userdata');
   hb=findobj(h0,'tag','TimeSlider');
   set(hb, 'value', ud.from);
   i_movie('play', 1);
   i_movie('stopavi');
 case 'deactivate'
   ud = get(h0, 'userdata');
   ud.active = 0;
   ud.step = 0;
   set(h0, 'userdata',ud);
 case 'hide'
   i_movie('deactivate');
   hidebar(par1, h0);
   hmenu=findobj('tag','hidebuttons');
   if strcmpi(par1, 'off')
      set(hmenu,'label','Show buttons','callback','i_movie(''hide'',''on'')');
   else
      set(hmenu,'label','Hide buttons','callback','i_movie(''hide'',''off'')');
   end;
 case 'play'
   ud = get(h0, 'userdata');
   ud.step = par1;
   set(h0, 'userdata',ud);
   if ~ud.active
      i_movie('activate');
   end;
 case 'step'
   ud = get(h0, 'userdata');
   ud.step = 0;
   set(h0, 'userdata',ud);
   if ~ud.active
      hb=findobj(h0,'tag','TimeSlider');
      t1 = get(hb, 'value');
      t1 = t1 + par1;
      set(hb, 'value', t1);
      i_movie('activate');
   end;
 case 'sett'
   ud = get(h0, 'userdata');
   ud.step = 0;
   set(h0, 'userdata',ud);
   hb=findobj(h0,'tag','TimeSlider');
   t1 = round((par1 - t) / (g_grind.ndays / g_grind.tstep)) + 1;
   set(hb, 'value', t1)
   i_movie('activate')
 case 'activate'
   if nargin == 1
      par1 = h0;
   end;
   ud = get(par1, 'userdata');
   ud.active = 1;
   set(par1, 'userdata',ud);
   ud.step = 1;
   hb=findobj(par1,'tag','TimeSlider');

   t1 = round(get(hb, 'value'));
   if isempty(t1) || (t1 <= 0)
      t1 = 1;
   end;

   while ud.step ~= 0
      switch ud.type
       case 'movplot'
          set(hb,'max',size(ud.varno.Xs,1), 'min', 1);
          if t1>size(ud.varno.Xs,1)
             t1=size(ud.varno.Xs,1);
          end;
          xlim1 = get(gca, 'xlim');
          ylim1 = get(gca, 'ylim');
          zlim1 = get(gca, 'zlim');
%          h1 = get(gca, 'children');
%          ls = myget(h1, 'LineStyle');
%          lw = myget(h1, 'LineWidth');
%          mr = myget(h1, 'Marker');
%          cl = myget(h1, 'Color');
%          xl=get(get(gca, 'xlabel'),'string');
%          yl=get(get(gca, 'ylabel'),'string');
          ch=get(gca,'children');

         if (length(ch)>1)&&strcmp(ud.varno.type,'stem3') %R11 problem, you cannot set xdata etc
           x= ud.varno.Xs(t1,:);
           y= ud.varno.Ys(t1,:);
           z= ud.varno.Zs(t1,:);
           x1= [x;x;NaN+zeros(size(x))];
           y1 = [y;y;NaN+zeros(size(y))];
           z1 = [zeros(size(z));z;NaN+zeros(size(z))];
           set(ch(1),'xdata',x);
           set(ch(2),'xdata',x1(:)');
           set(ch(1),'ydata',y);
           set(ch(2),'ydata',y1(:)');
           set(ch(1),'zdata',z);
           set(ch(2),'zdata',z1(:)');
        elseif (length(ch)>1)&&strcmp(ud.varno.type,'stem') %R11 problem, you cannot set xdata etc
           x= ud.varno.Xs(t1,:);
           y= ud.varno.Ys(t1,:);
           x1= [x;x;NaN+zeros(size(x))];
           y1 = [zeros(size(y));y;NaN+zeros(size(y))];
           set(ch(2),'xdata',x);
           set(ch(1),'xdata',x1(:)');
           set(ch(2),'ydata',y);
           set(ch(1),'ydata',y1(:)');
         else
         if ~isempty(ud.varno.Xs)
             set(ch,'xdata', ud.varno.Xs(t1,:));
         end;
         if ~isempty(ud.varno.Ys)
             set(ch,'ydata', ud.varno.Ys(t1,:));
         end;
         if ~isempty(ud.varno.Zs)
             set(ch,'zdata', ud.varno.Zs(t1,:));
         end;
         end;
%          if isempty(ud.varno.Zs) & isempty(ud.varno.Ys);
%             feval(ud.varno.type, ud.varno.Xs(t1, :));
%          elseif isempty(ud.varno.Zs)
%             feval(ud.varno.type, ud.varno.Xs(t1, :), ud.varno.Ys(t1, :), 'o');
%          else
%             feval(ud.varno.type, ud.varno.Xs(t1, :), ud.varno.Ys(t1, :), ud.varno.Zs(t1, :));
%          end;
%          i_plotdefaults;
%          h1 = get(gca, 'children');
%          xlabel(xl);
%          ylabel(yl);
%          myset(h1, 'LineStyle', ls);
%          myset(h1, 'LineWidth', lw);
%          myset(h1, 'Marker', mr);
%          myset(h1, 'Color', cl);
          set(gca, 'xlim', xlim1);
          set(gca, 'ylim', ylim1);
          if ~isempty(ud.varno.Zs)
             set(gca, 'zlim', zlim1);
          end;
       case 'varcontour'
         % g_'grind.viewcells.vars{ud.varno}(ie 'varcontour',:)
         %[sfrom, sto] = i_statevarnos(g_grind.statevars.vectnames{ud.varno});
         %[dim1, dim2] = i_getvardims(g_grind.statevars.vectnames{ud.varno});
         A1 = reshape(g_grind.viewcells.y{ud.varno}(t1, :)', g_grind.viewcells.vars{ud.varno}.dims(1), ...
            g_grind.viewcells.vars{ud.varno}.dims(2));
         if g_grind.viewcells.vars{ud.varno}.dims(2) == 1
            A1 = repmat(A1, 1, 3);
         end;
         set(ud.handle, 'cdata', A1);
       case 'lines'
         [sfrom] = i_statevarnos(g_grind.statevars.vectnames{ud.varno});
         [dim1] = i_getvardims(g_grind.statevars.vectnames{ud.varno});
         xl = get(gca, 'xlim');
         yl = get(gca, 'ylim');
         A1 = reshape(g_Y(t1,sfrom:sfrom + dim1 - 1)',1, dim1);
         xbins = getxbins(ud.varno);
         if isempty(xbins)
            plot(A1);
         else
            plot(xbins, A1);
         end;
         xlim(xl);
         ylim(yl);
         i_plotdefaults;
       case 'bars'
         [sfrom] = i_statevarnos(g_grind.statevars.vectnames{ud.varno});
         [dim1] = i_getvardims(g_grind.statevars.vectnames{ud.varno});
         xl = get(gca, 'xlim');
         yl = get(gca, 'ylim');
         A1 = reshape(g_Y(t1,sfrom:sfrom + dim1 - 1)',1, dim1);
         xbins = getxbins(ud.varno);
         if isempty(xbins)
            bar(A1);
         else
            bar(xbins, A1);
         end;
         xlim(xl);
         ylim(yl);
         i_plotdefaults;
      end;
      if t1 < 1 || t1 > length(g_t)
         set(hb,'tooltipstring',sprintf('step = %d',t1));
      else
         set(hb,'tooltipstring',sprintf('t = %0.5g',g_t(t1)));
      end;
      if ud.savingAVI
         i_makeavi('getframe');
      end;
      pause(ud.pause);
      if ishandle(hb)
         t1 = round(get(hb, 'value'));
         maxt1 = get(hb, 'max');
         ud = get(par1, 'userdata');
         t1 = t1 + ud.step;
         if t1 < 1
            t1 = 1;
            ud.step = 0;
            i_movie('deactivate');
         elseif (t1 > maxt1)  || (isfield(ud, 'to') && (t1>ud.to))
            t1 = maxt1;
            ud.step = 0;
            i_movie('deactivate');
         else
            set(hb, 'value', t1);
         end;
      else
         return;
      end;
   end;
   ud.active = 0;
   set(par1, 'userdata',ud);
end;
function hidebar(par1, h0)
hb(1)=findobj(h0,'tag','TimeSlider');
hb(2)=findobj(h0,'tag','hb2');
hb(3)=findobj(h0,'tag','hb3');
hb(4)=findobj(h0,'tag','hb4');
hb(5)=findobj(h0,'tag','hb5');
hb(6)=findobj(h0,'tag','hb6');
hb(7)=findobj(h0,'tag','hb7');
set(hb, 'visible', par1);
function xbins = getxbins(i)
global g_grind;
if isfield(g_grind.statevars, 'xbins')
   xbins=evalin('base',evalin('base',sprintf('g_grind.statevars.xbins{%d}',i)));
else
   xbins = [];
end;
% function res = myget(varargin)
% try
%    res = get(varargin{:});
% catch
%    res = [];
% end;
% function myset(varargin)
% try
%    set(varargin{:});
% catch
% end;




