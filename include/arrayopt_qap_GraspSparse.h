/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class arrayopt_qap_GraspSparse */

#ifndef _Included_arrayopt_qap_GraspSparse
#define _Included_arrayopt_qap_GraspSparse
#ifdef __cplusplus
extern "C" {
#endif
#undef arrayopt_qap_GraspSparse_DEFAULT_MAX_ITERACTIONS
#define arrayopt_qap_GraspSparse_DEFAULT_MAX_ITERACTIONS 100L
#undef arrayopt_qap_GraspSparse_DEFAULT_ALPHA
#define arrayopt_qap_GraspSparse_DEFAULT_ALPHA 0.25f
#undef arrayopt_qap_GraspSparse_DEFAULT_BETA
#define arrayopt_qap_GraspSparse_DEFAULT_BETA 0.5f
#undef arrayopt_qap_GraspSparse_DEFAULT_SEED
#define arrayopt_qap_GraspSparse_DEFAULT_SEED 270001L
/*
 * Class:     arrayopt_qap_GraspSparse
 * Method:    qap_grasps
 * Signature: (IIFFI[I[I[I[I)J
 */
JNIEXPORT jlong JNICALL Java_arrayopt_qap_GraspSparse_qap_1grasps
  (JNIEnv *, jobject, jint, jint, jfloat, jfloat, jint, jintArray, jintArray, jintArray, jintArray);

#ifdef __cplusplus
}
#endif
#endif
