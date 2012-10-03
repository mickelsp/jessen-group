%INITGRIND   Initiate grind
%   initiate global variables. If you use use or model, you never
%   need to call this function directly. Only used in combination
%   with setodefile.
%
%   See also use, model, setodefile

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
global t g_Y g_t g_data g_grind g_paranal g_func;
addpath([grindpath filesep 'sys2']);
g_grind.commands= {};
g_grind.model= {};
g_grind.scheme = {};
g_grind.funcs = [];
g_grind.funcnames.names ={};
g_func=[];
g_grind.inifile='';
g_grind.odefile='';
g_grind.onstart.funs={};
g_grind.pars = {};
g_grind.ndays = 1000;
g_grind.tstep = NaN;
g_grind.Jacobian = {};
g_grind.timevars = {};
g_grind.xaxis.var = [];
g_grind.yaxis.var = [];
g_grind.zaxis.var = [];
g_grind.statevars.names = {};
g_grind.statevars.vectnames={};
g_grind.statevars.dims={};
g_grind.statevars.vector=0;
g_grind.boxcar.names={};
g_grind.boxcar.trains={};
g_grind.pen = nextpen([]);
g_grind.diffto='Time (t)';
g_grind.outt={};
g_grind.lastsettings = [];
g_grind.truncate = 0;
g_grind.externvars={};
g_data=[];
g_paranal.Y=[];
g_paranal.t=[];
g_paranal.p=[];
g_grind.solver.name = 'ode45';
g_grind.solver.isdiffer = 0;
g_grind.solver.opt = odeset;
g_grind.solver.opt.RelTol=5e-5;
g_grind.solver.opt.AbsTol=1e-7;
%g_grind.solver.opt.InitialStep=0.001;
g_grind.solver.iters=1;
g_grind.solver.backwards = 0;
g_grind.solver.addmode = 0;
g_grind.solver.hasevents = 0;
g_grind.solver.haslag = 0;
g_grind.solver.nonautonomous = 0;
g_grind.solver.reset_at_data = 0;
[l_r,l_asys]=getrelease;
%g_grind.version.isoctave=strcmpi(l_asys,'octave');
% if g_grind.version.isoctave
%   warning off; %#ok
% end;
warning('off','backtrace');
g_grind.version.matlabrelease=l_r;
clear l_r l_asys;
t = 0;
%g_noise=[];
g_Y = [];
g_t = [];
setdefaults(1);
