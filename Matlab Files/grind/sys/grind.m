%GRIND for MATLAB
%   GRIND is a convenient DOS program for analyzing sets of 
%   differential equations made by Rob de Boer.
%   GRIND for MATLAB is a MATLAB version of GRIND, which also supports
%   difference equations, delay differential equations (DDE),  vectors, 
%   matrices (cellular automata).
%   GRIND for MATLAB is a command based system, i.e. the user types commands
%   in the MATLAB command Window to do most analyses and make figures. 
%   These figures can be edited using standard MATLAB commands and menus. 
%   Additionally, a user can click in some figures for instance to show 
%   trajectories in a phase plane. 
%   The system includes commands for:
%   - simulating sets of n ordinary differential equations, including time 
%     delays (see lag) and vector- and matrix notations. 
%   - simulating sets of difference equations.
%   - creating null-isoclines (nullclines) in phase spaces of 2 or 3 dimensions. (see null, null3)
%   - creating nullclines in spaces spanned by state variables and parameters. 
%   - creating one-dimensional bifurcation plots by simulation. (see paranal)
%   - stability analysis of equilibria, using eigenvalues of the Jacobian matrix. (see eigen)
%   - creating plots of user-defined functions (including parameters from 
%     the current model)(see funplot)
%   - using (red) noise and interpolation of external variables (for example 
%     real temperature data) (see rednoise)
%   - automatic calibration of parameters by optimizing the sum of squares 
%     between observed and predicted values(see optimpars).
%   - various special analyses such as determining the Lyapunov coefficient
%     to detect chaos, Poincaré sections etc. (see lyapunov)
%   - discrete events within a continuous model (see setevent)
%   - stochastic differential equations (see dwiener)
%   - three kinds of "boxcartrains" (Goudriaan, 1989) to model stage structured populations (see boxcartrain)
%
%
%   Use the command model to create a model (which is saved to an ini
%   file). Use the upper panel in the window to enter differential 
%   equation or difference equations. You can also define the model as a Forrester diagram 
%   (see vismod).
%
%   Example:
%   Logistic differential equation:
%     N'=N*r*(1-N/K)
%   Logistic difference equation:
%      N(t+1)=N(t)*r*(1-N(t)/K)
%
%   Note that parameters and state variables are case-sensitive.
%   Parameter names can be all alphanumeric names except the
%   following reserved words:
%      t   time
%     pi   pi=3.1416
%     inf  Infinity
%     Inf  Infinity
%     nan  NaN (not-a-number)
%     NaN  NaN (not-a-number)
%     eps  Floating point relative accuracy eps=2.2204e-016 
%
%   In the lower panel default values of the parameters and initial
%   values of the state variables are entered.  
%   Furthermore default commands can be entered here.
%
%   Example:
%     N=0.01;
%     r=0.5;
%     ax x N [0 10];
%   (the semicolon is not required, but suppresses unnecessary
%   output) 
%
%   Author: 
%
%   Egbert van Nes (egbert.vanNes@wur.nl)
%   Wageningen University
%   Aquatic Ecology and Water Quality Management group
%   PO Box 47
%   6700 AA Wageningen 
%   The Netherlands
%      
%
%   GRIND for Matlab website:
%   http://www.aew.wur.nl/UK/grind/ 
%   See also installing, model

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:25 $
ver(fullfile('grind','sys'));
if isempty(strfind(pwd,'grind'))
   disp(' ');
   disp(['changed directory to: ' grindpath(2)]);
   cd(grindpath(2))
end;
addpath([grindpath filesep 'sys2']);
disp(' ');
disp('Type "commands" for an overview of the commands.');
disp('Type "model" to create a model.');