function print_model(inputModelDir,inputModel,outputPath,orin,force)
% PRINT_MODEL - prints simulink models or library items to an EPS-file
%
%USAGE:
%PRINT_MODEL(inputModel,outputEps) opens the model to be printed, places an 
%EPS-file in the output path and file name specified in outputPath, and then
%closes the model again. 
%The model to be printed is the one specified in inputModel which is 
%located in inputModelDir. Remember to specify model name from top library, e.g. 
%'LCADlib/LiniarDynamics/A6-K1 -alpha Plant' with LCADlib as top library.
%The orinetation is determined by the string orin, e.g. 'portrait'
%To force the object open (if masked) set force = 1 else force = 0

%input check
if nargin < 3
    error('Need three inputs at least.');
end
if nargin == 3
    orin = 'portrait';
    force = 0;
end
if nargin == 4
   force = 0; 
end
%if orin ~= 'landscape' | orin ~= 'portrait' | orin ~= 'rotated' | orin ~= 'tall'
%    error('orin must be eigher: landscape, portrait, rotated, or tall.');
%end
if force ~= 0 && force ~= 1
    error('force must be eigher: 0 or 1.');
end
% function start
if force == 1
open_system(strcat(inputModelDir,inputModel),'force')
else
    open_system(strcat(inputModelDir,inputModel))
end
orient(inputModel,orin)
handle = get_param(gcs,'handle');
print(handle,'-depsc2',outputPath)
close_system(strcat(inputModel),0)
