#pragma once

#include "macro.hh"

namespace tinker { namespace resdue {
const int maxamino = 38;
const int maxnuc = 12;
extern int (&ntyp)[maxamino];
extern int (&catyp)[maxamino];
extern int (&ctyp)[maxamino];
extern int (&hntyp)[maxamino];
extern int (&otyp)[maxamino];
extern int (&hatyp)[maxamino];
extern int (&cbtyp)[maxamino];
extern int (&nntyp)[maxamino];
extern int (&cantyp)[maxamino];
extern int (&cntyp)[maxamino];
extern int (&hnntyp)[maxamino];
extern int (&ontyp)[maxamino];
extern int (&hantyp)[maxamino];
extern int (&nctyp)[maxamino];
extern int (&cactyp)[maxamino];
extern int (&cctyp)[maxamino];
extern int (&hnctyp)[maxamino];
extern int (&octyp)[maxamino];
extern int (&hactyp)[maxamino];
extern int (&o5typ)[maxnuc];
extern int (&c5typ)[maxnuc];
extern int (&h51typ)[maxnuc];
extern int (&h52typ)[maxnuc];
extern int (&c4typ)[maxnuc];
extern int (&h4typ)[maxnuc];
extern int (&o4typ)[maxnuc];
extern int (&c1typ)[maxnuc];
extern int (&h1typ)[maxnuc];
extern int (&c3typ)[maxnuc];
extern int (&h3typ)[maxnuc];
extern int (&c2typ)[maxnuc];
extern int (&h21typ)[maxnuc];
extern int (&o2typ)[maxnuc];
extern int (&h22typ)[maxnuc];
extern int (&o3typ)[maxnuc];
extern int (&ptyp)[maxnuc];
extern int (&optyp)[maxnuc];
extern int (&h5ttyp)[maxnuc];
extern int (&h3ttyp)[maxnuc];
extern char (&amino1)[maxamino][1];
extern char (&nuclz1)[maxnuc][1];
extern char (&amino)[maxamino][3];
extern char (&nuclz)[maxnuc][3];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(resdue, ntyp)[maxamino];
extern "C" int TINKER_MOD(resdue, catyp)[maxamino];
extern "C" int TINKER_MOD(resdue, ctyp)[maxamino];
extern "C" int TINKER_MOD(resdue, hntyp)[maxamino];
extern "C" int TINKER_MOD(resdue, otyp)[maxamino];
extern "C" int TINKER_MOD(resdue, hatyp)[maxamino];
extern "C" int TINKER_MOD(resdue, cbtyp)[maxamino];
extern "C" int TINKER_MOD(resdue, nntyp)[maxamino];
extern "C" int TINKER_MOD(resdue, cantyp)[maxamino];
extern "C" int TINKER_MOD(resdue, cntyp)[maxamino];
extern "C" int TINKER_MOD(resdue, hnntyp)[maxamino];
extern "C" int TINKER_MOD(resdue, ontyp)[maxamino];
extern "C" int TINKER_MOD(resdue, hantyp)[maxamino];
extern "C" int TINKER_MOD(resdue, nctyp)[maxamino];
extern "C" int TINKER_MOD(resdue, cactyp)[maxamino];
extern "C" int TINKER_MOD(resdue, cctyp)[maxamino];
extern "C" int TINKER_MOD(resdue, hnctyp)[maxamino];
extern "C" int TINKER_MOD(resdue, octyp)[maxamino];
extern "C" int TINKER_MOD(resdue, hactyp)[maxamino];
extern "C" int TINKER_MOD(resdue, o5typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, c5typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, h51typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, h52typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, c4typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, h4typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, o4typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, c1typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, h1typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, c3typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, h3typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, c2typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, h21typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, o2typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, h22typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, o3typ)[maxnuc];
extern "C" int TINKER_MOD(resdue, ptyp)[maxnuc];
extern "C" int TINKER_MOD(resdue, optyp)[maxnuc];
extern "C" int TINKER_MOD(resdue, h5ttyp)[maxnuc];
extern "C" int TINKER_MOD(resdue, h3ttyp)[maxnuc];
extern "C" char TINKER_MOD(resdue, amino1)[maxamino][1];
extern "C" char TINKER_MOD(resdue, nuclz1)[maxnuc][1];
extern "C" char TINKER_MOD(resdue, amino)[maxamino][3];
extern "C" char TINKER_MOD(resdue, nuclz)[maxnuc][3];

int (&ntyp)[maxamino] = TINKER_MOD(resdue, ntyp);
int (&catyp)[maxamino] = TINKER_MOD(resdue, catyp);
int (&ctyp)[maxamino] = TINKER_MOD(resdue, ctyp);
int (&hntyp)[maxamino] = TINKER_MOD(resdue, hntyp);
int (&otyp)[maxamino] = TINKER_MOD(resdue, otyp);
int (&hatyp)[maxamino] = TINKER_MOD(resdue, hatyp);
int (&cbtyp)[maxamino] = TINKER_MOD(resdue, cbtyp);
int (&nntyp)[maxamino] = TINKER_MOD(resdue, nntyp);
int (&cantyp)[maxamino] = TINKER_MOD(resdue, cantyp);
int (&cntyp)[maxamino] = TINKER_MOD(resdue, cntyp);
int (&hnntyp)[maxamino] = TINKER_MOD(resdue, hnntyp);
int (&ontyp)[maxamino] = TINKER_MOD(resdue, ontyp);
int (&hantyp)[maxamino] = TINKER_MOD(resdue, hantyp);
int (&nctyp)[maxamino] = TINKER_MOD(resdue, nctyp);
int (&cactyp)[maxamino] = TINKER_MOD(resdue, cactyp);
int (&cctyp)[maxamino] = TINKER_MOD(resdue, cctyp);
int (&hnctyp)[maxamino] = TINKER_MOD(resdue, hnctyp);
int (&octyp)[maxamino] = TINKER_MOD(resdue, octyp);
int (&hactyp)[maxamino] = TINKER_MOD(resdue, hactyp);
int (&o5typ)[maxnuc] = TINKER_MOD(resdue, o5typ);
int (&c5typ)[maxnuc] = TINKER_MOD(resdue, c5typ);
int (&h51typ)[maxnuc] = TINKER_MOD(resdue, h51typ);
int (&h52typ)[maxnuc] = TINKER_MOD(resdue, h52typ);
int (&c4typ)[maxnuc] = TINKER_MOD(resdue, c4typ);
int (&h4typ)[maxnuc] = TINKER_MOD(resdue, h4typ);
int (&o4typ)[maxnuc] = TINKER_MOD(resdue, o4typ);
int (&c1typ)[maxnuc] = TINKER_MOD(resdue, c1typ);
int (&h1typ)[maxnuc] = TINKER_MOD(resdue, h1typ);
int (&c3typ)[maxnuc] = TINKER_MOD(resdue, c3typ);
int (&h3typ)[maxnuc] = TINKER_MOD(resdue, h3typ);
int (&c2typ)[maxnuc] = TINKER_MOD(resdue, c2typ);
int (&h21typ)[maxnuc] = TINKER_MOD(resdue, h21typ);
int (&o2typ)[maxnuc] = TINKER_MOD(resdue, o2typ);
int (&h22typ)[maxnuc] = TINKER_MOD(resdue, h22typ);
int (&o3typ)[maxnuc] = TINKER_MOD(resdue, o3typ);
int (&ptyp)[maxnuc] = TINKER_MOD(resdue, ptyp);
int (&optyp)[maxnuc] = TINKER_MOD(resdue, optyp);
int (&h5ttyp)[maxnuc] = TINKER_MOD(resdue, h5ttyp);
int (&h3ttyp)[maxnuc] = TINKER_MOD(resdue, h3ttyp);
char (&amino1)[maxamino][1] = TINKER_MOD(resdue, amino1);
char (&nuclz1)[maxnuc][1] = TINKER_MOD(resdue, nuclz1);
char (&amino)[maxamino][3] = TINKER_MOD(resdue, amino);
char (&nuclz)[maxnuc][3] = TINKER_MOD(resdue, nuclz);
#endif
} }
