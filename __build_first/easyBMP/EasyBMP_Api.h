#ifndef __EasyBMP_Api_H__
#define __EasyBMP_Api_H__

#ifndef __GNUC__
  #ifdef EASYBMP_EXPORTS
    #define EASYBMP_API __declspec(dllexport)
    #define EASYBMP_API_C extern "C" EASYBMP_API
  #else
    #define EASYBMP_API __declspec(dllimport)
    #define EASYBMP_API_C extern "C" EASYBMP_API
  #endif
#else
  #define EASYBMP_API
  #define EASYBMP_API_C
#endif

#endif
