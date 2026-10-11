#pragma STDC FENV_ACCESS ON
#include <dlfcn.h>
#include <fenv.h>
#include <math.h>
#include <mpfr.h>
#include <stdio.h>
typedef long double (*powfn)(long double,long double);
struct test {const char *name; long double x,y; int flags;};
int main(int argc,char **argv){
 if(argc!=2) return 2;
 void *h=dlopen(argv[1],RTLD_NOW|RTLD_LOCAL);
 if(!h){fprintf(stderr,"%s\n",dlerror());return 3;}
 powfn fn=(powfn)dlsym(h,"powl");
 Dl_info info;
 if(!fn || !dladdr((void *)fn,&info)) return 4;
 printf("binding=%s\n",info.dli_fname);fflush(stdout);
 struct test tests[]={
 {"general-negative-odd-overflow",-0x1p-16445L,-1.L,FE_OVERFLOW|FE_INEXACT},
 {"general-negative-even-overflow",-0x1p-16445L,-2.L,FE_OVERFLOW|FE_INEXACT},
 {"general-negative-odd-underflow",-0.5L,32767.L,FE_UNDERFLOW|FE_INEXACT},
 {"general-negative-even-underflow",-0.5L,32766.L,FE_UNDERFLOW|FE_INEXACT},
 {"early-negative-odd-overflow",-2.L,32769.L,FE_OVERFLOW|FE_INEXACT},
 {"early-negative-odd-underflow",-2.L,-32769.L,FE_UNDERFLOW|FE_INEXACT},
 {"general-positive-overflow",0x1p-16445L,-1.L,FE_OVERFLOW|FE_INEXACT},
 {"general-positive-underflow",0.5L,32767.L,FE_UNDERFLOW|FE_INEXACT},
 {"integer-negative-odd-overflow",-2.L,32767.L,FE_OVERFLOW|FE_INEXACT},
 {"integer-negative-odd-underflow",-2.L,-32767.L,FE_UNDERFLOW|FE_INEXACT},
 {"integer-negative-even-overflow",-2.L,32766.L,FE_OVERFLOW|FE_INEXACT},
 {"integer-negative-even-underflow",-2.L,-32766.L,FE_UNDERFLOW|FE_INEXACT},
 {"integer-positive-odd-overflow",2.L,32767.L,FE_OVERFLOW|FE_INEXACT},
 {"integer-positive-odd-underflow",2.L,-32767.L,FE_UNDERFLOW|FE_INEXACT},
 {"integer-positive-even-overflow",2.L,32766.L,FE_OVERFLOW|FE_INEXACT},
 {"integer-positive-even-underflow",2.L,-32766.L,FE_UNDERFLOW|FE_INEXACT}
 };
 int modes[]={FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO};
 const char *names[]={"nearest","downward","upward","towardzero"};
 mpfr_rnd_t rnd[]={MPFR_RNDN,MPFR_RNDD,MPFR_RNDU,MPFR_RNDZ};
 mpfr_t a,b,z;mpfr_inits2(256,a,b,z,(mpfr_ptr)0);
 unsigned count=0;
 for(unsigned i=0;i<sizeof(tests)/sizeof(tests[0]);i++){
  mpfr_set_ld(a,tests[i].x,MPFR_RNDN);mpfr_set_ld(b,tests[i].y,MPFR_RNDN);
  mpfr_pow(z,a,b,MPFR_RNDN);
  for(int m=0;m<4;m++){
   if(fesetround(modes[m])) return 5;
   long double expected=mpfr_get_ld(z,rnd[m]);
   volatile long double vx=tests[i].x,vy=tests[i].y;
   feclearexcept(FE_ALL_EXCEPT);
   long double got=fn(vx,vy);
   int flags=fetestexcept(FE_ALL_EXCEPT)&(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW|FE_UNDERFLOW|FE_INEXACT);
   int ok=got==expected && !!signbit(got)==!!signbit(expected) && flags==tests[i].flags;
   printf("%s mode=%s got=%La expected=%La flags=0x%x expected_flags=0x%x %s\n",tests[i].name,names[m],got,expected,flags,tests[i].flags,ok?"PASS":"FAIL");
   fflush(stdout);count++;
   if(!ok){printf("STOP first_failure completed=%u\n",count);return 1;}
  }
 }
 printf("PASS completed=%u\n",count);
 mpfr_clears(a,b,z,(mpfr_ptr)0);
 return 0;
}
