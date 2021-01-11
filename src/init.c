#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#include <Rinternals.h> // for SEXP
#include <R_ext/RS.h>

void F77_NAME(aniawsim)(int* y, int* n1, int* n2, int* dv, double* ani,
  double* hakt, double* lambda, int* theta, double* bi, int* thnew, int* kern,
  int* skern, double* spmin, double* spmax, double* wght, double* swjy);
void F77_NAME(aniawsv)(int* y, int* n1, int* n2, int* dv, double* ani,
  double* vcoef, int* nvpar, double* meanvar, double* chcorr, double* hakt,
  double* lambda, int* theta, double* bi, int* thnew, int* kern, int* skern,
  double* spmin, double* spmax, double* wghts, double* swjy);
void F77_NAME(awsimg)(int* y, int* n1, int* n2, int* dv, double* hakt,
  int* thnew, double* bi, int* kern, double* lw, double* swjy);
void F77_NAME(awsimg0)(int* y, int* n1, int* n2, int* dv, double* hakt,
  int* thnew, double* bi, int* kern, double* lw, double* swjy);
void F77_NAME(awspimg)(int* y, int* n1, int* n2, int* dv, int* degr,
  double* hw, double* vcoef, int* nv, double* mvar, double* hakt,
  double* lambda, double* theta, double* bi, double* bi0, double* ai, int* kern,
  double* spmin, double* lw, double* w, double* slw, double* sw, double* wght,
  int* ind);
void F77_NAME(awsvimg0)(int* y, int* fix, int* n1, int* n2, int* n, int* dv,
   double* vcoef, int* nvpar, double* meanvar, double* chcorr, double* hakt,
   double* hhom, double* lambda, int* theta, double* bi, double* bi0,
   int* thnew, int* kern, double* spmin, double* wghts, double* lw, int* early,
   int* homogen);
void F77_NAME(cam2rgb)( int* theta, int* n, double* outcam, int* thetanew);
void F77_NAME(connect1)( int* segm, int* n1, int* n2, int* i1, int* i2,
  int* ind1, int* ind2, int* checked);
void F77_NAME(convolve)( double* img, double* kernel, double* eimg,
  int* height, int* width, int* ksize);
void F77_NAME(cumhist)( int* x, int* n, int* yi);
void F77_NAME(demmed16)( int*  sensor, int* theta, int* n1, int* n2, int* h1,
  int* h2, int* bayer);
void F77_NAME(demmed4)( int*  sensor, int* theta, int* n1, int* n2, int* h1,
  int* h2, int* bayer);
void F77_NAME(demmed4b)( int*  sensor, int* theta, int* n1, int* n2, int* h1,
  int* h2, int* bayer);
void F77_NAME(dhomogen)( int* img, int* n1, int* n2, int* imghom, int* bayer);
void F77_NAME(epsigmac)( int* y, int* n, int* dv, int* theta, double* bi,
  int* quant, double* varcoef, double* mvar, int* dp1);
void F77_NAME(epsigmal)( int* y, int* n, int* dv, int* theta, double* bi,
  int* quant, double* varcoef, double* mvar, int* dp1  );
void F77_NAME(esigmac)( int* y, int* n, int* dv, int* theta, double* bi,
  int* quant, double* varcoef, double* mvar);
void F77_NAME(esigmal)( int* y, int* n, int* dv, int* theta, double* bi,
  int* quant, double* varcoef, double* mvar);
void F77_NAME(esigmaq)( int* y, int* n, int* dv, int* theta, double* bi,
  int* quant, double* varcoef, double* mvar);
void F77_NAME(estcorr)( double* res, int* n1, int* n2, int* dv, double* scorr,
  double* chcorr);
void F77_NAME(fullsize)( int*  sensor, int* theta, int* n1, int* n2, int* h1,
  int* h2, int* bayer);
void F77_NAME(geth2)( double* x, double* y, int* kern, double* value,
  double* eps, double* bw);
void F77_NAME(getvofh2)( double* bw, int* kern, double* vol);
void F77_NAME(halfsize)( int*  sensor, int* theta, int* n1, int* n2, int* h1,
  int* h2, int* bayer);
void F77_NAME(hequalc)( int* x, int* n, int* y, int* yi);
void F77_NAME(hequalg)( int* x, int* n, int* y, int* yi);
void F77_NAME(ihequal)( int* x, int* n, int* y, int* yi);
void F77_NAME(ihequalc)( int* x, int* n, int* y, int* yi);
void F77_NAME(indemos4)( int* sensor, int* theta, int* n1, int* n2, int* bayer,
  int* bi, int* bi3);
void F77_NAME(mawsimg0)(int* y, int* fix, int* mask, int* n1, int* n2, int* dv,
  double* hakt, double* lambda, int* theta, double* bi, double* bi0, int* thnew,
  int* kern, double* spmin, double* lw, double* wght);
void F77_NAME(median1)(double* x, int* n, double* y, double* tol);
void F77_NAME(median3)(double* x, int* n, double* y, double* tol);
void F77_NAME(mpaws2)(int* n, int* dp1, int* dp2, double* ai, double* bi,
  double* theta, double* dmat, int* ind);
void F77_NAME(segment)(int* y, double* level, double* delta, int* n1, int* n2,
  double* hakt, double* lambda, int* theta, double* vcoef, int* nvc,
  double* meanvar, double* bi, double* s2i, int* thnew, int* kern, double* spmin,
  double* lw, double* pvalue, int* segm, double* thresh, double* fov,
  double* varest);
void F77_NAME(senvar)( int* s, int* n1, int* n2, int* shat, double* bi,
  int* bayer, double* vcoef, double* mvar, int* nothom);
void F77_NAME(shrinkc)(int* x, int* nx1, int* nx2, int* y, int* ny1, int* ny2,
  double* tol, double* z, int* nz, int* method, int* nc);
void F77_NAME(shrnkcsp)(double* img, int* nx, int* ny, int* dv, double* imgnew,
  int* nxnew, int* nynew, int* indx, int* indy, int* method);
void F77_NAME(shrinkg)(int* x, int* nx1, int* nx2, int* y, int* ny1, int* ny2,
  double* tol, double* z, int* nz, int* method, int* nc);
void F77_NAME(shrnkgr)(int* img, int* nx, int* ny, int* imgnew,
  int* nxnew, int* nynew, int* indx, int* indy, int* method);
void F77_NAME(shrnkrgb)(int* img, int* nx, int* ny, int* dv, int* imgnew,
  int* nxnew, int* nynew, int* indx, int* indy, int* method);
void F77_NAME(sm2dtens)(double* x, int* n1, int* n2, double* h, double* rho,
  double* xhat);
void F77_NAME(smsens0)(int* s, int* shat, double* bi, int* n1, int* n2,
  int* bayer);
void F77_NAME(smsensor)(int* s, int* shat, int* th, int* n1, int* n2, int* nt1,
  int* nt2, int* bayer, double* vcoef, double* meanvar, double* hakt,
  double* lambda, double* bi, int* kern, double* spmin, double* lw);
void F77_NAME(wbalance)( int* sensor, int* n1, int* n2, double* wb, int* bayer);

static R_NativePrimitiveArgType aniawsim_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType aniawsv_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType awsimg_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType awsimg0_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType awspimg_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP};
static R_NativePrimitiveArgType awsvimg0_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP,
  INTSXP};
static R_NativePrimitiveArgType cam2rgb_t[]={ INTSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType connect1_t[]={ INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType convolve_t[]={ REALSXP, REALSXP, REALSXP,
  INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType cumhist_t[]={ INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType demmed16_t[]={ INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType demmed4_t[]={ INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType demmed4b_t[]={ INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType dhomogen_t[]={ INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP};
static R_NativePrimitiveArgType epsigmac_t[]={ INTSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType epsigmal_t[]={ INTSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, INTSXP  };
static R_NativePrimitiveArgType esigmac_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
   REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType esigmal_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
   REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType esigmaq_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
   REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType estcorr_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, REALSXP};
static R_NativePrimitiveArgType fullsize_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType geth2_t[]={ REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP, REALSXP};
static R_NativePrimitiveArgType getvofh2_t[]={REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType halfsize_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType hequalc_t[]={ INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType hequalg_t[]={ INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType ihequal_t[]={ INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType ihequalc_t[]={ INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType indemos4_t[]={ INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType mawsimg0_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType median1_t[]={REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType median3_t[]={REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType mpaws2_t[]={INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType segment_t[]={INTSXP, REALSXP, REALSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType senvar_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType shrinkc_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType shrnkcsp_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType shrinkg_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType shrnkgr_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType shrnkrgb_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType sm2dtens_t[]={ REALSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP};
static R_NativePrimitiveArgType smsens0_t[]={INTSXP, INTSXP, REALSXP, INTSXP, INTSXP,
  INTSXP};
static R_NativePrimitiveArgType smsensor_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType wbalance_t[]={ INTSXP, INTSXP, INTSXP, REALSXP,
  INTSXP};

static const R_FortranMethodDef fmethods[] = {
            {"aniawsim", (DL_FUNC) &aniawsim_ ,16 , aniawsim_t},
            {"aniawsv", (DL_FUNC) &aniawsv_ , 20, aniawsv_t},
            {"awsimg", (DL_FUNC) &awsimg_ , 10, awsimg_t},
            {"awsimg0", (DL_FUNC) &awsimg0_ , 10, awsimg0_t},
            {"awspimg", (DL_FUNC) &awspimg_ , 23, awspimg_t},
            {"awsvimg0", (DL_FUNC) &awsvimg0_ , 23, awsvimg0_t},
            {"connect1", (DL_FUNC) &connect1_ , 8, connect1_t},
            {"cumhist", (DL_FUNC) &cumhist_ , 3, cumhist_t},
            {"demmed16", (DL_FUNC) &demmed16_ , 7, demmed16_t},
            {"demmed4", (DL_FUNC) &demmed4_ , 7, demmed4_t},
            {"demmed4b", (DL_FUNC) &demmed4b_ , 7, demmed4b_t},
            {"dhomogen", (DL_FUNC) &dhomogen_ , 5, dhomogen_t},
            {"esigmac", (DL_FUNC) &esigmac_ , 8, esigmac_t},
            {"esigmal", (DL_FUNC) &esigmal_ , 8, esigmal_t},
            {"esigmaq", (DL_FUNC) &esigmaq_ , 8, esigmaq_t},
            {"epsigmac", (DL_FUNC) &epsigmac_ , 9, epsigmac_t},
            {"epsigmal", (DL_FUNC) &epsigmal_ , 9, epsigmal_t},
            {"estcorr", (DL_FUNC) &estcorr_ , 6, estcorr_t},
            {"fullsize", (DL_FUNC) &fullsize_ , 7, fullsize_t},
            {"geth2", (DL_FUNC) &geth2_ , 6, geth2_t},
            {"getvofh2", (DL_FUNC) &getvofh2_ , 3, getvofh2_t},
            {"halfsize", (DL_FUNC) &halfsize_ , 7, halfsize_t},
            {"hequalc", (DL_FUNC) &hequalc_ , 4, hequalc_t},
            {"hequalg", (DL_FUNC) &hequalg_ , 4, hequalg_t},
            {"ihequal", (DL_FUNC) &ihequal_ , 4, ihequal_t},
            {"ihequalc", (DL_FUNC) &ihequalc_ , 4, ihequalc_t},
            {"indemos4", (DL_FUNC) &indemos4_ , 7, indemos4_t},
            {"mawsimg0", (DL_FUNC) &mawsimg0_ , 16, mawsimg0_t},
            {"median1", (DL_FUNC) &median1_ ,4, median1_t},
            {"median3", (DL_FUNC) &median3_ ,4, median3_t},
            {"mpaws2", (DL_FUNC) &mpaws2_ , 8, mpaws2_t},
            {"segment", (DL_FUNC) &segment_ , 22, segment_t},
            {"senvar", (DL_FUNC) &senvar_ , 9, senvar_t},
            {"shrinkc", (DL_FUNC) &shrinkc_ ,11 , shrinkc_t},
            {"shrinkg", (DL_FUNC) &shrinkg_ ,11 , shrinkg_t},
            {"sm2dtens", (DL_FUNC) &sm2dtens_ , 6, sm2dtens_t},
            {"wbalance", (DL_FUNC) &wbalance_ , 5, wbalance_t},
            {"cam2rgb", (DL_FUNC) &cam2rgb_ , 4, cam2rgb_t},
            {"convolve", (DL_FUNC) &convolve_ , 6, convolve_t},
            {"shrnkcsp", (DL_FUNC) &shrnkcsp_ , 10, shrnkcsp_t},
            {"shrnkgr", (DL_FUNC) &shrnkgr_ ,9, shrnkgr_t},
            {"shrnkrgb", (DL_FUNC) &shrnkrgb_ , 10, shrnkrgb_t},
            {"smsensor", (DL_FUNC) &smsensor_ , 16, smsensor_t},
            {"smsens0", (DL_FUNC) &smsens0_ , 6, smsens0_t},
            {NULL, NULL, 0,NULL}
};

void R_init_adimpro(DllInfo *dll)
         {
             R_registerRoutines(dll, NULL, NULL, fmethods , NULL);
             R_useDynamicSymbols(dll,FALSE);
             R_forceSymbols(dll,TRUE);
         }
