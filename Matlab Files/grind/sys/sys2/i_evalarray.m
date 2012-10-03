%script i_evalarray
%g_l_s1 is the local variable with the commands
g_l_s1=strrep(g_l_s1,'.*','*');
g_l_s1=strrep(g_l_s1,'*','.*');
g_l_s1=strrep(g_l_s1,'.^','^');
g_l_s1=strrep(g_l_s1,'^','.^');
g_l_s1=strrep(g_l_s1,'./','/');
g_l_s1=strrep(g_l_s1,'/','./');
eval(g_l_s1);
