% script used for compilation of mex-files (should be used when you got at
% version of Matlab newer than 6.5

mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/qmul.c ../lib/utils/math_utils/qmul_wrapper.c 
mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/qestimator.c

mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/s_cross_product.c
mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/s_norm.c
mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/s_pow.c
mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/s_modulus.c
mex -outdir ../lib/utils/math_utils ../lib/utils/math_utils/s_sign.c

mex -outdir ../lib/utils/ephemeris ../lib/utils/ephemeris/igrfS.c 
mex -outdir ../lib/utils/ephemeris ../lib/utils/ephemeris/sgp4S.c 

mex -outdir ../lib/utils ../lib/utils/zonal.c
mex -outdir ../lib/utils ../lib/utils/zonalcatesian.c
mex -outdir ../lib/utils ../lib/utils/scdynamics.c

mex -outdir ../lib/utils/acs ../lib/utils/acs/desaturator.c
