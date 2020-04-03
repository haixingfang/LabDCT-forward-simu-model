
% load DIPimage toolbox
% http://www.diplib.org/
function load_diplib
	addpath(strcat(pwd,'\dipimage_2.9_win64\dip\common\dipimage'));
    dip_initialise;
    dipsetpref('imagefilepath',strcat('\dipimage_2.9_win64\dip\images'));
end
