function [m] = i_table(flag, A)
switch flag
 case 'init'
   if ~exist('uitable','builtin')
      varcopy(A);
      m = A;
      disp('a table is copied to the clipboard');
   else
      f = figure('MenuBar','none', ...
          'Color',[0.914 0.914 0.914], ...
         'NumberTitle','off', ...
         'Name','Table',...
    'CreateFcn',@(h,evnt)movegui(gcbf, 'center'));
      uitable(f,'Units','normalized','Position',...
         [0 0.1 1 0.9], 'Data', A, ...
         'ColumnEditable', true, ...
         'Tag','Table');
      uicontrol('Parent',f, ...
         'Units','normalized', ...
         'Callback','i_table(''cancel'')', ...
         'ListboxTop',0, ...
         'Position',[0.4,0.02,0.09,0.06], ...
         'String','Cancel', ...
         'Tag','Cancelbutton');
      uicontrol('Parent',f, ...
         'Units','normalized', ...
         'Callback','i_table(''ok'')', ...
         'ListboxTop',0, ...
         'Position',[0.5,0.02,0.09,0.06], ...
         'String','OK', ...
         'Tag','OKbutton');
      uicontrol('Parent',f, ...
         'Units','normalized', ...
         'Callback','i_table(''clipboard'')', ...
         'ListboxTop',0, ...
         'Position',[0.6,0.02,0.09,0.06], ...
         'String','Clipboard', ...
         'Tag','Clipboardbutton');
      uicontrol('Parent',f, ...
         'Units','normalized', ...
         'Callback','i_table(''paste'')', ...
         'ListboxTop',0, ...
         'Position',[0.7,0.02,0.09,0.06], ...
         'String','Paste', ...
         'Tag','Pastebutton');
      uicontrol('Parent',f, ...
         'Units','normalized', ...
         'Callback','i_table(''addcol'')', ...
         'ListboxTop',0, ...
         'Position',[0.8,0.02,0.09,0.06], ...
         'String','Add column', ...
         'Tag','AddColbutton');
      uicontrol('Parent',f, ...
         'Units','normalized', ...
         'Callback','i_table(''addrow'')', ...
         'ListboxTop',0, ...
         'Position',[0.9,0.02,0.09,0.06], ...
         'String','Add row', ...
         'Tag','AddRowbutton');
     if nargout>0
      set(gcf,'userdata','uiwait');
      uiwait;
      ud = get(f, 'userdata');
      if ~isempty(ud) && strcmp(ud, 'OK')
         m=get(findobj(f,'Tag','Table'),'data');
      else
         m = [];
      end;
      close(f);
     end;
   end;
 case 'addcol'
   dat=get(findobj(gcbf,'Tag','Table'),'data');
   if iscell(dat)
       dat{1,end+1}=[];
   else
       dat(1,end+1)=0;
   end;
   set(findobj(gcbf,'Tag','Table'),'data',dat)
 case 'addrow'
   dat=get(findobj(gcbf,'Tag','Table'),'data');
   if iscell(dat)
       dat{end+1,1}=[];
   else
       dat(end+1,1)=0;
   end;
   set(findobj(gcbf,'Tag','Table'),'data',dat)
 case 'clipboard'
   dat=get(findobj(gcbf,'Tag','Table'),'data');
   varcopy(dat);
 case 'paste'
   set(findobj(gcbf,'Tag','Table'),'data',varpaste);
 case 'ok'
   ud=get(gcbf,'userdata');
   set(gcbf,'userdata','OK');
   if strcmp(ud,'uiwait')
     uiresume;
   else
      close;
   end;
 case 'cancel'
   ud=get(gcbf,'userdata');
   set(gcbf,'userdata','cancel');
   if strcmp(ud,'uiwait')
     uiresume;
   else
      close;
   end;
end;