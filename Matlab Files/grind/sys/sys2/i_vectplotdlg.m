function fig = i_vectplotdlg(flag)
global g_grind;
if nargin == 0
   flag = 'init';
end;
switch flag
case 'init'
   if ~isempty(g_grind)&&~isfield(g_grind,'vectplot')
      vectplot('-defaults');
   end;      
   if ~isempty(g_grind)
       N=length(g_grind.vectplot.vars);
       while (N>0)&&isempty(g_grind.vectplot.vars{N})
          N=N-1;
       end;
    end;
    if N<1
       N=1;
    end;
    plotlst = cell(1,N);
    for i = 1:N
      plotlst{i} = sprintf('Plot No. %d', i);
   end;
   chartypes={'bar','bar3','pcolor','plot','plot3','scatter','scatter3','stem','stem3','surface'};
h0 = figure('Units','points', ...
	'Color',[0.914 0.914 0.914], ...
	'MenuBar','none', ...
	'Name','Edit Plots for Vectplot', ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperType','A4', ...
	'PaperUnits','points', ...
	'Position',[305 302 366 245], ...
	'Tag','vectplotdlg', ...
	'ToolBar','none',...
    'CreateFcn',@(h,evnt)movegui(gcbf, 'center'));
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.914 0.914 0.914], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[17.25 138.75 303.75 15], ...
	'String','Y-axis: t or one auxiliary or state variable', ...
	'Style','text', ...
	'Tag','StaticText2');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Max',30, ...
	'Position',[17 126 317 17], ...
	'Style','edit', ...
	'Tag','EditY', ...
	'TooltipString','Enter here a list of variables or functions for the Y-axis');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.914 0.914 0.914], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[17 206 72.75 15], ...
	'String','Vectplot No.', ...
	'Style','text', ...
	'Tag','StaticText1');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','i_vectplotdlg(''changeplot'')', ...
	'ListboxTop',0, ...
	'Position',[105 205.5 113.25 22.5], ...
	'String',plotlst, ...
	'Style','popupmenu', ...
	'Tag','PlotNo', ...
	'TooltipString','Select Plot number here', ...
	'Value',1);
 uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_vectplotdlg(''newplot'')', ...
	'ListboxTop',0, ...
	'Position',[243 205.5 93.75 22.5], ...
	'String','Add Vector Plot', ...
	'Tag','AddBtn', ...
	'TooltipString','Add new plot to the list of vectplots');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_vectplotdlg(''ok'')', ...
	'ListboxTop',0, ...
	'Position',[31.5 8.25 93.75 22.5], ...
	'String','OK', ...
	'Tag','OKBtn', ...
	'TooltipString','Close and save');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_vectplotdlg(''cancel'')', ...
	'ListboxTop',0, ...
	'Position',[141 8.25 93.75 22.5], ...
	'String','Cancel', ...
	'Tag','CancelBtn', ...
	'TooltipString','Close and undo changes');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','commands vectplot', ...
	'ListboxTop',0, ...
	'Position',[250.5 8.25 93.75 22.5], ...
	'String','Help', ...
	'Tag','HelpBtn', ...
	'TooltipString','Open help window');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.914 0.914 0.914], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[17 180 303.75 15], ...
	'String','X-axis: t or an auxliiary or state variable', ...
	'Style','text', ...
	'Tag','StaticText2');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[17 168 318 16.5], ...
	'Style','edit', ...
	'Tag','EditX', ...
	'TooltipString','Edit the variable for the X axis here');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.914 0.914 0.914], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[17 97.5 303.75 15], ...
	'String','Z-axis, one auxiliary or state variable', ...
	'Style','text', ...
	'Tag','StaticText2');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Max',30, ...
	'Position',[17 84.75 317.25 16.5], ...
	'Style','edit', ...
	'Tag','EditZ', ...
	'TooltipString','Enter here a variable or function for the Z-axis');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[106.5 42 113.25 22.5], ...
	'String',chartypes, ...
	'Style','popupmenu', ...
	'Tag','TypeNo', ...
	'TooltipString','Select Time Plot here', ...
	'Value',1);
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.914 0.914 0.914], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[17 48 72.75 15], ...
	'String','Kind of plot', ...
	'Style','text', ...
	'Tag','StaticText1');
   %if N==1
  %   set(hplotno,'Enable','off');
  % end;
%    global g_grind; 
    if ~isempty(g_grind)
       ud.oudlist=g_grind.vectplot;
    else
       ud.oudlist={};
    end;
    ud.curr=g_grind.vectplot.currno;
    set(h0,'userdata',ud);
   setplot(ud.curr,h0);
   if nargout > 0, fig = h0; end
 case 'ok'
%    global g_grind;
    ud=get(gcbf,'userdata');
    savecurrent(ud.curr);
    delete(gcbf);
    disp(' ');
    vectplot('-list');
    vectplot;
 case 'cancel'
    ud=get(gcbf,'userdata');
%    global g_grind;
    g_grind.vectplot=ud.oudlist;
    delete(gcbf);
    disp(' ');
    vectplot('-list');
 case 'newplot'
   h=findobj(gcbf,'Tag','PlotNo');
   plotlst = get(h,'string');
   No=length(plotlst)+1;
   plotlst{No} = sprintf('Plot No. %d', No);
   if No>1
     set(h,'Enable','on');
   end;
   set(h,'string',plotlst);
   set(h,'value',No);
   i_vectplotdlg('changeplot');
case 'changeplot'
    ud=get(gcbf,'userdata');
    savecurrent(ud.curr);
    h=findobj(gcbf,'tag','PlotNo');
  v=get(h,'Value');
  setplot(v,gcbf);
  ud.curr=v;
  set(gcbf,'userdata',ud);
end;

function savecurrent(No)
global g_grind;
f=gcbf;
h=findobj(f,'tag','EditX');
g_grind.vectplot.vars{No}.xax=outf('changeshortcut', get(h,'string'));
h=findobj(f,'tag','EditY');
g_grind.vectplot.vars{No}.yax=outf('changeshortcut', get(h,'string'));
h=findobj(f,'tag','EditZ');
g_grind.vectplot.vars{No}.zax=outf('changeshortcut', get(h,'string'));
h=findobj(f,'tag', 'TypeNo');
typenames=get(h,'string');
i=get(h,'value');
g_grind.vectplot.vars{No}.type=typenames{i};
g_grind.vectplot.currno=No;

function setplot(No,f)
global g_grind;
if isempty(g_grind)||(No>length(g_grind.vectplot.vars))
   type='';
   x='t';
   y='';
   z='';
else
   x = g_grind.vectplot.vars{No}.xax;
   y = g_grind.vectplot.vars{No}.yax;
   z = g_grind.vectplot.vars{No}.zax;
   type = g_grind.vectplot.vars{No}.type;
end;
h=findobj(f,'tag','EditX');
set(h,'string',x);
h=findobj(f,'tag','EditY');
set(h,'string',y);
h=findobj(f,'tag','EditZ');
set(h,'string',z);
h=findobj(f,'tag', 'PlotNo');
set(h,'value',No);
h=findobj(f,'tag', 'TypeNo');
typenames=get(h,'string');
i=1;
while i<=length(typenames)&&~strcmpi(typenames{i},type)
   i=i+1;
end;
if i>length(typenames)
   i=1;
end;
set(h,'value',i);


