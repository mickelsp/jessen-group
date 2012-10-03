function i_draw2Darrow(x, y,sumlen,lenarrow)
global g_grind;
if nargin<3
    sumlen=0.15;
end;
if nargin<4
    lenarrow=10;
end;
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
if (isempty(iX.no) || isempty(iY.no))
   ax('?');
   errordlg('Cannot create arrows if there are no state variables on the axes, use "phas 2" instead.');
   error('GRIND:null:NoStatevars','null: Cannot create arrows if there are no state variables on the axes, use "phas 2" instead.');
end
if nargin==0
  NO=i_initvar;
  x=NO(iX.no);
  y=NO(iY.no);
end;
[dum, dum, Vect] = i_vector(1, iX, [x x], iY, [y y], [],0);
dx=Vect(1,1,1);
dy=Vect(2,2,2);
lenx=abs(dx/(g_grind.xaxis.lim(2)-g_grind.xaxis.lim(1)));
leny=abs(dy/(g_grind.yaxis.lim(2)-g_grind.yaxis.lim(1)));
minratio=2;
if minratio*lenx<leny
    dx=dx*leny/(minratio*lenx);
    lenx=abs(dx/(g_grind.xaxis.lim(2)-g_grind.xaxis.lim(1)));
elseif minratio*leny<lenx
    dy=dy*lenx/(minratio*leny);
    leny=abs(dy/(g_grind.yaxis.lim(2)-g_grind.yaxis.lim(1)));
end;
dx=dx*sumlen/(lenx+leny);
dy=dy*sumlen/(lenx+leny);

i_makefig('phase2');
hold on;
axannotation('arrow',[x x+dx],[y y], 'Tag','Arrow','HeadLength',lenarrow);
axannotation('arrow',[x x],[y y+dy], 'Tag','Arrow','HeadLength',lenarrow);


