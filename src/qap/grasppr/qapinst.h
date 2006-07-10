#ifndef _QAPINST_H
#define _QAPINST_H

typedef struct {
   int n, **F, **D;
} qap_inst;

// ***********************************************
// Commented out: function to create new QAP
//                instance reading from file
//
// qap_inst *qi_read_new(char *fname);
//
// Added: function to create QAP instance
//        gathering data from one-dimensional
//        arrays
//
qap_inst *qi_new(int dim, int *f, int *d);
//
// -Sergio A. de Carvalho Jr.
// ***********************************************

void qi_delete(qap_inst *q);

// ***********************************************
// Commented out: functions not needed
//
// int qi_sparsity(qap_inst *q);
// int qi_is_symmetric(qap_inst *q);
//
// -Sergio A. de Carvalho Jr.
// ***********************************************

#endif
