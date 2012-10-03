function fig = i_moddlg(doclear, docheck, isvisible, doera)
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
%
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.
% 'ButtonDownFcn','i_cmoddlg(''keydown'')',...
global g_grind;
g=grindpath(2);
if isempty(g_grind)&&~strcmp(pwd, g)
   cd(g);
   fprintf('changed directory to %s\n',g);
end;
fontsize=11;
font='Courier New';
if nargin == 0
   doclear = 0;
end;
if nargin <2
   docheck = 1;
end;
h0=i_figno('dialog');
for i=1:10
   if ishandle(h0+i)&&strcmp(get(h0+i,'tag'),'GRIND Model')
      close(h0+i);
   end;
end;
if ishandle(h0)
   close(h0);
   while ishandle(h0)
      h0=h0+1;
   end;   
end; 
h0= figure(h0);
if ~isvisible
   set(h0,'visible','off');
end;

set(h0, 'WindowButtonDown', '');
set(h0, 'WindowButtonMotionFcn', '');
set(h0,'Color',[0.941 0.941 0.941], ...
	'MenuBar','none', ...
	'Name','GRIND Model', ...
	'NumberTitle','off', ...
	'PaperPosition',[1296 12960 41472 31104], ...
   'PaperUnits','points', ...
   'ResizeFcn','i_cmoddlg(''resize'')',...
	'Position',get(0,'defaultfigurepos'), ...
	'Tag','GRIND Model', ...
   'ToolBar','none');
set(h0,'KeypressFcn','i_cmoddlg(''keydown'')')
set(h0,	'Units','points');
pos=round(get(h0,'position'));
pos(2)=pos(2)+25;
pos(3)=430;
pos(4)=320;
set(h0,'position',pos);
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.941 0.941 0.941], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[11 281 106 16], ...
	'String','Name of model:', ...
	'Style','text', ...
	'Tag','namelabel');
 uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','i_cmoddlg(''updated'')', ...
	'HorizontalAlignment','left', ...
   'ListboxTop',0, ...
 	'FontName',font, ...
	'FontSize',fontsize, ...
	'Position',[130 284 181 16], ...
	'Style','edit', ...
	'Tag','FileName', ...
   'TooltipString','Enter here the name of the ini file');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_cmoddlg(''cd'')', ...
	'ListboxTop',0, ...
	'Position',[314 284 14 16], ...
	'String','...', ...
	'Tag','CDButton', ...
   'TooltipString','Select directory');
%328-14=314
if ~isempty(g_grind)&&(length(g_grind.scheme)>10)
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.941 0.941 0.941], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
   'Position',[11 255 295 16], ...
	'String','Model equations (use vismod to change)', ...
	'Style','text', ...
   'Tag','modellabel');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_cmoddlg(''clearScheme'')', ...
	'ListboxTop',0, ...
	'Position',[200 260 70 16], ...
	'String','Clear scheme', ...
	'Tag','ClearSchemeButton', ...
	'TooltipString','Clear the vismod scheme');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','i_cmoddlg(''updated'')', ...
	'FontName',font, ...
	'FontSize',fontsize, ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Max',100, ...
   'Position',[11 143 317 116], ...
   'Enable','off', ...
	'Style','edit', ...
	'Tag','Model', ...
	'TooltipString','Enter here the differential equations (N''=..)');
else
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.941 0.941 0.941], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
   'Position',[11 255 295 16], ...
	'String','Model equations (diffential/difference equations) ', ...
	'Style','text', ...
   'Tag','modellabel');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','i_cmoddlg(''updated'')', ...
	'FontName',font, ...
	'FontSize',fontsize, ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Max',100, ...
   'Position',[11 143 317 116], ...
	'Style','edit', ...
	'Tag','Model', ...
	'TooltipString','Enter here the differential equations (N''=..)');
end
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.941 0.941 0.941], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[11 120 295 16], ...
	'String','Parameters / Initial conditions / Commands', ...
	'Style','text', ...
	'Tag','parameterslabel');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','i_cmoddlg(''updated'')', ...
	'FontName',font, ...
	'FontSize',fontsize, ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Max',100, ...
	'Position',[11 6 318 118], ...
	'Style','edit', ...
	'Tag','Parameters', ...
	'TooltipString','Assign default values to the parameters here');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_cmoddlg(''extract'')', ...
	'ListboxTop',0, ...
	'Position', [339 145 83 25], ...
	'String','Extract parameters', ...
	'Tag','ExtractButton', ...
	'TooltipString','Extract a list of parameters from the model eq.');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_cmoddlg(''cancel'')', ...
	'ListboxTop',0, ...
	'Position',[339 46 83 25], ...
	'String','Cancel', ...
	'Tag','CancelButton', ...
	'TooltipString','Cancel, don''t change/use model');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_cmoddlg(''ok'')', ...
	'ListboxTop',0, ...
	'Position', [339 13 83 25], ...
	'String','OK', ...
	'Tag','OKButton', ...
	'TooltipString','Save file and create model');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_cmoddlg(''openfile'')', ...
	'ListboxTop',0, ...
 	'Position',[339 276 83 25], ...
	'String','Load model', ...
	'Tag','OpenFileButton', ...
   'TooltipString','Load an existing model');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','i_cmoddlg(''newmodel'')', ...
	'ListboxTop',0, ...
	'Position',[339 111  83 25], ...
	'String','New model', ...
	'Tag','NewButton', ...
	'TooltipString','Clear all edit fields');
uicontrol('Parent',h0, ...
	'Units','points', ...
	'ButtonDownFcn',')', ...
	'Callback','commands modelpanel', ...
	'ListboxTop',0, ...
	'Position',[339 78 83 25], ...
	'String','Help', ...
	'Tag','Helpbutton');
if nargout > 0, fig = h0; end
i_cmoddlg('init');
ud = get(h0, 'Userdata');
ud.clear = doclear;
ud.check = docheck;
ud.era=doera;
ud.debug = 0;
set(gcf, 'Userdata', ud);

