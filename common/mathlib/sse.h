//========= Copyright Valve Corporation, All rights reserved. ============//
//
// Purpose: 
//
//=====================================================================================//

#ifndef _SSE_H
#define _SSE_H

#include <aa_luse.h>

extern float _SSE_Sqrt( float x );
extern float _SSE_RSqrtAccurate( float a );
extern float _SSE_RSqrtFast( float x );
extern float FASTCALL _SSE_VectorNormalize( LVector3& vec );
extern void FASTCALL _SSE_VectorNormalizeFast( LVector3& vec );
extern float _SSE_InvRSquared( const float* v );
extern void _SSE_SinCos( float x, float* s, float* c );
extern float _SSE_cos( float x );
extern void _SSE2_SinCos( float x, float* s, float* c );
extern float _SSE2_cos( float x );
extern void VectorTransformSSE( const float *in1, const LMatrix4& in2, float *out1 );
extern void VectorRotateSSE( const float *in1, const LMatrix4& in2, float *out1 );

#endif // _SSE_H