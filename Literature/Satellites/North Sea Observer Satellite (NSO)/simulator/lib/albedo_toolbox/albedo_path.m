% ALBEDO_PATH Returns the path to albedo_path.m, under the assumption that this% is the albedo toolbox installation dir.
%
%  pathstr = albedo_path%% $Id: albedo_path.m,v 1.1 2006/03/17 12:13:50 mnkr02 Exp $function pathstr = albedo_path;[pathstr,name,ext,versn] = fileparts(which('albedo_path.m'));return
