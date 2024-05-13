% script used for compilation of mex-files (should be used when you got at
% version of Matlab newer than 6.5

mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/qmul.c ../lib/utils/math_utils/qmul_wrapper.c 
mex -outdir ../lib/utils/ephemeris ../lib/utils/ephemeris/igrfS.c 
mex -outdir ../lib/utils/ephemeris ../lib/utils/ephemeris/sgp4S.c
mex -outdir ../lib/utils/ephemeris ../lib/utils/math_utils/qestimator.c
mex -outdir ../lib/utils/ephemeris ../lib/utils/scdynamics.c
mex -outdir ../lib/utils/ephemeris ../lib/utils/zonal.c
mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/svdmex.c


/space_simulator_Wiki_Version/