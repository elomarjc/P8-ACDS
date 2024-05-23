% script used for compilation of mex-files (should be used when you got at
% version of Matlab newer than 6.5

% Make sure you use the one for your operating system and outcomment the
% other one.

%%% WINDOWS USERS %%%
% mex -outdir ..space_simulator_wiki_version/lib/utils/math_utils ..space_simulator_wiki_version/lib/utils/math_utils/qmul.c ..space_simulator_wiki_version/lib/utils/math_utils/qmul_wrapper.c 
% mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/ephemeris/igrfS.c 
% mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/ephemeris/sgp4S.c
% mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/math_utils/qestimator.c
% mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/scdynamics.c
% mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/zonal.c
% mex -outdir ..space_simulator_wiki_version/lib/utils/math_utils ..space_simulator_wiki_version/lib/utils/math_utils/svdmex.c

%%% MacOS USERS %%%
mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/qmul.c ../lib/utils/math_utils/qmul_wrapper.c 
mex -outdir ../lib/utils/ephemeris ../lib/utils/ephemeris/igrfS.c 
mex -outdir ../lib/utils/ephemeris ../lib/utils/ephemeris/sgp4S.c
mex -outdir ../lib/utils/ephemeris ../lib/utils/math_utils/qestimator.c
mex -outdir ../lib/utils/ephemeris ../lib/utils/scdynamics.c
mex -outdir ../lib/utils/ephemeris ../lib/utils/zonal.c
mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/svdmex.c