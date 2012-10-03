function cls
ax = get(0,'CurrentFigure');
if isempty(ax)
   error('GRIND:cls:noFigure','No figure available');
else
   ax=get(ax,'CurrentAxes');
end;
if ~isempty(ax)
series = get(ax, 'children');
if ~isempty(series)>0
   found=0;
   for i=1:length(series)
      if strcmp(get(series(i),'selected'),'on')
         delete(series(i));
         found=1;
         break;
      end;
   end;
   if ~found
      delete(series(1));
   end;
   refresh;
else
   error('GRIND:cls:NoSeries','No series');
end;
else
   error('GRIND:cls:NoAxis','Figure has no axis');
end;
