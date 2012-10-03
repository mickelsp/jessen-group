%FIGDATA   Extract xyz data from a figure
%  Use this function to export data from the currently selected 
%  figure to a matrix. 
%
%
%  Usage:
%  M=FIGDATA - copies the first series to matrix M. If there is a series 
%  selected, that series is copied.
%  FIGDATA -MERGE - merges two selected series.
%  FIGDATA -ALL - gets one matrix of all series.
%  FIGDATA ? - gives an overview of the series of the current figure
%  M=FIGDATA(N) - copies the N'th series to M.
%  FIGDATA N - shows the N'th series.
%  
%  See also copyfig, varcopy

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function A = figdata(serno)
if nargin == 0
   serno=1;
   flag='';
else
   if ischar(serno)&&~strcmp(serno,'?')&&~strcmp(serno,'-merge')&&~strcmp(serno,'-all');
      serno = i_checkstr(serno);
      flag='';
   else 
      flag=serno;
   end;
end;
ax = get(0,'CurrentFigure');
if isempty(ax)
   error('GRIND:figdata:NoFig','No figure available');
else
   ax=get(ax,'CurrentAxes');
end;
series = get(ax, 'children');
if nargin==0
   sel=strcmp(get(series,'selected'),'on');
   if any(sel)
      serno=find(sel,1);
      fprintf('Using selected series, series %d\n',serno);
   else
      disp('Using series 1');
   end;
end;
if strcmp(flag,'-merge');
   s=find(strcmp(get(series,'selected'),'on'));
   if length(s)<2
      error('GRIND:figdata:SelSeries','select two series to merge');
   end;
   X1=getv(series(s(1)),'xdata');
   Y1=getv(series(s(1)),'ydata');
   Z1=getv(series(s(1)),'zdata');
   X2=getv(series(s(2)),'xdata');
   Y2=getv(series(s(2)),'ydata');
   Z2=getv(series(s(2)),'zdata');
   if X1(1)==X2(1)
      X1=flipud(X1);
      Y1=flipud(Y1);
      Z1=flipud(Z1);
   end;
   set(series(s(1)),'xdata',[X1;X2]); 
   set(series(s(1)),'ydata',[Y1;Y2]); 
   set(series(s(1)),'zdata',[Z1;Z2]);
   fprintf('series %d and %d merged\n',s(1),s(2));
   return
elseif strcmp(flag,'-all');
   if isempty(series)
      A=[];
      return;
   end;
   X1=getv(series(1),'xdata');
   Y1=getv(series(1),'ydata');
   Z1=getv(series(1),'zdata');
   YZ1=[Y1,Z1];
   for i=2:length(series)
      YZ2=[getv(series(i),'ydata'),getv(series(i),'zdata')];
      X2=getv(series(i),'xdata');
      [X1,YZ1]=i_concatdata(X1,YZ1,X2,YZ2);
   end;
   A=[X1,YZ1];
   return;   
elseif ~isempty(flag)&&(flag == '?')
   fprintf('The current figure has %d series\n', length(series));
   for i=1:length(series)
      X=getv(series(i),'xdata');
      Y=getv(series(i),'ydata');
      Z=getv(series(i),'zdata');
      s=sprintf('Length series%3d: X: %6d, Y: %6d, Z: %6d',i,length(X),length(Y),length(Z));
      if strcmp(get(series(i),'selected'),'on')
         s=sprintf('%s <== selected',s);
      end;
      disp(s);
   end;   
   return;
elseif (serno <= length(series)) && (serno > 0)
   curser = series(serno);
else
   error('GRIND:figdata:UnknownSeries','Series %d doesn''t exist',serno);
end;
X = getv(curser, 'xdata');
Y = getv(curser, 'ydata');
Z = getv(curser, 'zdata');
if isempty(X) && isempty(Y) && isempty(Z)
   A = [];
elseif ~isempty(X) && ~isempty(Y) && isempty(Z)
   A = [X,Y];
elseif ~isempty(X) && ~isempty(Y) && ~isempty(Z)
   A=[X,Y,Z];
end;
function Ser=getv(h,dataax)
Ser=get(h,dataax);
if size(Ser,1)==1
   Ser=Ser';
end;
