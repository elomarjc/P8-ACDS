% script used for compilation of mex-files (should be used when you got at
% version of Matlab newer than 6.5

mex -outdir ..space_simulator_wiki_version/lib/utils/math_utils ..space_simulator_wiki_version/lib/utils/math_utils/qmul.c ..space_simulator_wiki_version/lib/utils/math_utils/qmul_wrapper.c 
mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/ephemeris/igrfS.c 
mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/ephemeris/sgp4S.c
mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/math_utils/qestimator.c
mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/scdynamics.c
mex -outdir ..space_simulator_wiki_version/lib/utils/ephemeris ..space_simulator_wiki_version/lib/utils/zonal.c
mex -outdir ..space_simulator_wiki_version/lib/utils/math_utils ..space_simulator_wiki_version/lib/utils/math_utils/svdmex.c
