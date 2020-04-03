% load Multi-Parametric Toolbox 3
% https://www.mpt3.org/
function load_mpt3
	addpath(strcat(pwd,'\tbxmanager'));
    tbxmanager restorepath;
    mpt_init;
    % % cite using MPT3:
    % % M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, pages 502?10, Zurich, Switzerland, July 17?9 2013.
end