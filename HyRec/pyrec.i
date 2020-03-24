%module pyrec
%{
#include "pyrec.h"
%}

%typemap(in) double[3](double temp[3]) {
  if (PyTuple_Check($input)) {
    if (!PyArg_ParseTuple($input,"ddd",temp,temp+1,temp+2)) {
      PyErr_SetString(PyExc_TypeError,"tuple must have 3 elements");
      SWIG_fail;
    }
    $1 = &temp[0];
  } else {
    PyErr_SetString(PyExc_TypeError,"expected a tuple.");
    SWIG_fail;
  }
 }

%typemap(in) double[63](double temp[63]) {
  if (PyTuple_Check($input)) {
    if (!PyArg_ParseTuple($input,"ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd",
			  temp,temp+1,temp+2,temp+3,temp+4,temp+5,temp+6,temp+7,temp+8,temp+9,
			  temp+10,temp+11,temp+12,temp+13,temp+14,temp+15,temp+16,temp+17,temp+18,temp+19,
			  temp+20,temp+21,temp+22,temp+23,temp+24,temp+25,temp+26,temp+27,temp+28,temp+29,
			  temp+30,temp+31,temp+32,temp+33,temp+34,temp+35,temp+36,temp+37,temp+38,temp+39,
			  temp+40,temp+41,temp+42,temp+43,temp+44,temp+45,temp+46,temp+47,temp+48,temp+49,
			  temp+50,temp+51,temp+52,temp+53,temp+54,temp+55,temp+56,temp+57,temp+58,temp+59,
			  temp+60,temp+61,temp+62)) {
      PyErr_SetString(PyExc_TypeError,"tuple must have 63 elements");
      SWIG_fail;
    }
    $1 = &temp[0];
  } else {
    PyErr_SetString(PyExc_TypeError,"expected a tuple.");
    SWIG_fail;
  }
 }


extern void rec_build_history_wrap(double tcmb, double obh2, double odmh2, double okh2, double odeh2, 
				   double w0, double wa, double yp, double nnu, double mnu[3],
				   double INJion[63],double INJexc[63],double INJheat[63]);
extern double hyrec_xe(double a);
extern double hyrec_tm(double a);
