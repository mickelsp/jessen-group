function rotatefig(h)
if nargin==0
    h=gca;
end;
v=get(h,'view');
fig=get(h,'parent');
figure(fig);
aviobj = avifile('rotatefig.avi');
aviobj.quality = 100;
for i=-90:90
    set(h,'view',[i,v(2)]);
    drawnow;
    pause(0.01);
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
end;
close(aviobj);
disp('written to rotatefig.avi');

