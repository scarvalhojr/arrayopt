#include "qapinst.h"
#include "common.h"

// ******************************************************************
// Commented out: no use of following libraries here
//
// #include "randgen.h"
// #include <stdio.h>
// #include <stdlib.h>
//
// Added: use of memset function in string library was implicit
//
#include <string.h>
//
// Code writing to the standard output was deliberately commented out
// since, when integrated into arrayopt, no output will be seen.
//
// -Sergio A. de Carvalho Jr.
// ******************************************************************

/* --------------------------------------------------------------
 * The QAP instance structure and functions
 * --------------------------------------------------------------*/

// ******************************************************************
// Commented out: no need for input/output functions
//
// /* read a instance from the file fname and returns
//  * an instance structure filled with information. */
// qap_inst *qi_read_new(char *fname)
// {
//    int i, j;
//    FILE *f;
//    qap_inst *qi = New(qap_inst);
//    if (!qi) return 0;
//    f =  fopen(fname, "r");
//    if (!f) {
//       fprintf(stderr, "error: cannot open file %s\n", fname);
//       return 0;
//    }
//    /* read the instance */
//    fscanf(f, "%d", &qi->n);
//    if (qi->n < 1) {
//       fprintf(stderr, "error: ilegal dimension %d\n", qi->n);
//       return 0;
//    }
//    NewVec2zero(qi->F, int, qi->n, qi->n);
//    if (!qi->F) return 0;
//    NewVec2zero(qi->D, int, qi->n, qi->n);
//    if (!qi->D) return 0;
//    for (i=0; i<qi->n; i++)
//       for (j=0; j<qi->n; j++)
// 	 fscanf(f, "%d", &qi->D[i][j]);
//    for (i=0; i<qi->n; i++)
//       for (j=0; j<qi->n; j++)
// 	 fscanf(f, "%d", &qi->F[i][j]);
//    fclose(f);
//    return qi;
// }
//
// -Sergio A. de Carvalho Jr.
// ***********************************************

// ******************************************************************
// Added: new function that only creates an instance of a QAP and,
//        instead of reading data from a file, it gathers data from
//        one dimensional arrays; it substitutes the qi_read_new
//        function that has been commented out
//
qap_inst *qi_new(int dim, int f[], int d[])
{
	int i, j;

	qap_inst *qi = New(qap_inst);
	if (!qi) return 0;

	qi->n = dim;

	NewVec2zero(qi->F, int, qi->n, qi->n);
	if (!qi->F) return 0;
	NewVec2zero(qi->D, int, qi->n, qi->n);
	if (!qi->D) return 0;

	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
	 		qi->D[i][j] = d[i * dim + j];

	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			qi->F[i][j] = f[i * dim + j];

	return qi;
}
//
// -Sergio A. de Carvalho Jr.
// ***********************************************

void qi_delete(qap_inst *q)
{
   if (q->D)
      DelVec2(q->D, q->n);
   if (q->F)
      DelVec2(q->F, q->n);
   Delete(q);
}

// ******************************************************************
// Commented out: functions not needed for arrayopt
//
// /* returns the % of sparcity of this instance */
// int qi_sparsity(qap_inst *q)
// {
//    int cd=0, cf=0, i, j, n = q->n, max;
//    for (i=0; i<n; i++)
//       for (j=0; j<n; j++) {
// 	 if (q->D[i][j] == 0) cd++;
// 	 if (q->F[i][j] == 0) cf++;
//       }
//    max = cd > cf ? cd : cf;
//    return (int)(100*(max/(float)(n*n)));
// }
//
// /* 1 if the instance is symmetric (matrices F and D).
//  * 0 otherwise */
// int qi_is_symmetric(qap_inst *q)
// {
//    int i, j, n = q->n;
//    for (i=0; i<=n/2; i++)
//       for (j=0; j<=n/2; j++) {
// 	 if (q->D[i][j] != q->D[j][i]) return 0;
// 	 if (q->F[i][j] != q->F[j][i]) return 0;
//       }
//    return 1;
// }
//
// -Sergio A. de Carvalho Jr.
// ***********************************************
