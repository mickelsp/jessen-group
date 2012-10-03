% i_back odefile for backward simulation
function Result = i_back(t, Input)
global g_grind;
Result = -feval(g_grind.odefile, t, Input);

