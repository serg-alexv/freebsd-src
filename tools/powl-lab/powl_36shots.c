#pragma STDC FENV_ACCESS ON
#include <dlfcn.h>
#include <fenv.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <mpfr.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
typedef long double (*powfn)(long double,long double);
static uint64_t sig(long double x){int e;long double f=frexpl(fabsl(x),&e);return (uint64_t)ldexpl(f,64);}
int main(int argc,char **argv){
 if(argc!=4 || LDBL_MANT_DIG!=64) return 2;
 void *h[3];powfn fn[3];const char *label[]={"system","PR","overlay"};
 for(int k=0;k<3;k++){
  h[k]=dlopen(argv[k+1],RTLD_NOW|RTLD_LOCAL);
  if(!h[k]){fprintf(stderr,"%s\n",dlerror());return 3;}
  fn[k]=(powfn)dlsym(h[k],"powl");Dl_info d;
  if(!fn[k] || !dladdr((void *)fn[k],&d)) return 4;
  printf("BIND %s %s\n",label[k],d.dli_fname);
  if(strcmp(d.dli_fname,argv[k+1])) return 5;
 }
 const long double x[]={0x3.439e117f80a71554p-13512L,0x7.5573852b2bddfbbp+16312L};
 const long double y[]={-0xe.7cc0e93023eb04p-16L,-0x6.c2fe326a74000008p-32L};
 mpfr_t a,b,z256,z512,g,delta,unit;
 mpfr_inits2(512,a,b,z512,g,delta,unit,(mpfr_ptr)0);mpfr_init2(z256,256);
 long double reference[2],first[2][3];int firstflags[2][3];
 if(fesetround(FE_TONEAREST)) return 6;
 for(int i=0;i<2;i++){
  mpfr_set_ld(a,x[i],MPFR_RNDN);mpfr_set_ld(b,y[i],MPFR_RNDN);
  mpfr_pow(z256,a,b,MPFR_RNDN);mpfr_pow(z512,a,b,MPFR_RNDN);
  reference[i]=mpfr_get_ld(z512,MPFR_RNDN);
  if(reference[i]!=mpfr_get_ld(z256,MPFR_RNDN)) return 7;
  printf("REF pair=%d x=%La y=%La ld80=%La sig=%016" PRIx64 " mpfr256_512_rounding_agree=1\n",i,x[i],y[i],reference[i],sig(reference[i]));
  mpfr_printf("REF512 pair=%d exact_approx=%Ra\n",i,z512);
 }
 unsigned n=0;
 for(int rep=0;rep<6;rep++)for(int i=0;i<2;i++)for(int k=0;k<3;k++){
  volatile long double vx=x[i],vy=y[i];
  if(fesetround(FE_TONEAREST)) return 8;
  feclearexcept(FE_ALL_EXCEPT);long double r=fn[k](vx,vy);
  int raw=fetestexcept(FE_ALL_EXCEPT),flags=raw&(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW|FE_UNDERFLOW|FE_INEXACT);
  if(!isfinite(r)) return 9;
  if(rep==0){first[i][k]=r;firstflags[i][k]=flags;}
  if(r!=first[i][k] || flags!=firstflags[i][k]) return 10;
  mpfr_set_ld(a,x[i],MPFR_RNDN);mpfr_set_ld(b,y[i],MPFR_RNDN);
  mpfr_pow(z512,a,b,MPFR_RNDN);mpfr_set_ld(g,r,MPFR_RNDN);
  mpfr_sub(delta,g,z512,MPFR_RNDN);
  mpfr_set_ui_2exp(unit,1,mpfr_get_exp(z512)-64,MPFR_RNDN);
  mpfr_div(delta,delta,unit,MPFR_RNDN);
  uint64_t rs=sig(r),es=sig(reference[i]),xr=rs^es;
  printf("SHOT=%u rep=%d pair=%d lib=%s r=%La flags=0x%x raw=0x%x sig=%016" PRIx64 " refsig=%016" PRIx64 " xor=%016" PRIx64 " common_leading_bits=%d differing_bits=%d ",++n,rep,i,label[k],r,flags,raw,rs,es,xr,xr?__builtin_clzll(xr):64,__builtin_popcountll(xr));
  mpfr_printf("signed_error_ulp=%.20Rg\n",delta);fflush(stdout);
 }
 printf("DONE shots=%u repetitions_stable=1\n",n);
 return 0;
}
