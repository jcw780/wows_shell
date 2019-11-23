//#ifdef __cplusplus
//extern "C" {
//#endif  /* __cplusplus */

#if defined (__GNUC__) || defined (__INTEL_COMPILER) || defined (__clang__)
#ifdef INFINITY
#undef INFINITY
#endif

#ifdef NAN
#undef NAN
#endif

#define NAN __builtin_nan("")
#define NANf __builtin_nanf("")
#define INFINITY __builtin_inf()
#define INFINITYf __builtin_inff()
#else

//#include <bits/nan.h>
//#include <bits/inf.h>

#endif

//#if defined __cplusplus
//}; /* End "C" */
//#endif  /* defined __cplusplus */