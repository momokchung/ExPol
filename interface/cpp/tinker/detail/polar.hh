#pragma once

#include "macro.hh"

namespace tinker { namespace polar {
extern int& npolar;
extern int*& ipolar;
extern double*& polarity;
extern double*& thole;
extern double*& dirdamp;
extern double*& pdamp;
extern double*& udir;
extern double*& udirp;
extern double*& udirs;
extern double*& udirps;
extern double*& uind;
extern double*& uinp;
extern double*& uinds;
extern double*& uinps;
extern double*& uexact;
extern int*& douind;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(polar, npolar);
extern "C" int* TINKER_MOD(polar, ipolar);
extern "C" double* TINKER_MOD(polar, polarity);
extern "C" double* TINKER_MOD(polar, thole);
extern "C" double* TINKER_MOD(polar, dirdamp);
extern "C" double* TINKER_MOD(polar, pdamp);
extern "C" double* TINKER_MOD(polar, udir);
extern "C" double* TINKER_MOD(polar, udirp);
extern "C" double* TINKER_MOD(polar, udirs);
extern "C" double* TINKER_MOD(polar, udirps);
extern "C" double* TINKER_MOD(polar, uind);
extern "C" double* TINKER_MOD(polar, uinp);
extern "C" double* TINKER_MOD(polar, uinds);
extern "C" double* TINKER_MOD(polar, uinps);
extern "C" double* TINKER_MOD(polar, uexact);
extern "C" int* TINKER_MOD(polar, douind);

int& npolar = TINKER_MOD(polar, npolar);
int*& ipolar = TINKER_MOD(polar, ipolar);
double*& polarity = TINKER_MOD(polar, polarity);
double*& thole = TINKER_MOD(polar, thole);
double*& dirdamp = TINKER_MOD(polar, dirdamp);
double*& pdamp = TINKER_MOD(polar, pdamp);
double*& udir = TINKER_MOD(polar, udir);
double*& udirp = TINKER_MOD(polar, udirp);
double*& udirs = TINKER_MOD(polar, udirs);
double*& udirps = TINKER_MOD(polar, udirps);
double*& uind = TINKER_MOD(polar, uind);
double*& uinp = TINKER_MOD(polar, uinp);
double*& uinds = TINKER_MOD(polar, uinds);
double*& uinps = TINKER_MOD(polar, uinps);
double*& uexact = TINKER_MOD(polar, uexact);
int*& douind = TINKER_MOD(polar, douind);
#endif
} }
