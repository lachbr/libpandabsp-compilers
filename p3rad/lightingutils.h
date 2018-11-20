/**
 * PANDA3D BSP TOOLS
 * Copyright (c) CIO Team. All rights reserved.
 *
 * @file lightingutils.h
 * @author Brian Lach
 * @date July 31, 2018
 *
 */

#ifndef LIGHTINGUTILS_H
#define LIGHTINGUTILS_H

#include <aa_luse.h>

#include "bspfile.h"
#include "qrad.h"
#include "bsptools.h"

class LightSurface : public BaseBSPEnumerator
{
public:
        LightSurface( int thread, bspdata_t *data );

        virtual bool enumerate_node( int node_id, const Ray &ray,
                                     float f, int context );

        virtual bool enumerate_leaf( int leaf_id, const Ray &ray, float start,
                                     float end, int context );

        virtual bool find_intersection( const Ray &ray );

private:
        bool test_point_against_surface( const LVector3 &point, dface_t *face, texinfo_t *tex );

        bool test_point_against_sky_surface( const LVector3 &point, dface_t *face );

public:
        int _thread;
        dface_t *_surface;
        float _hit_frac;
        LTexCoordf _luxel_coord;
        bool _has_luxel;
};

extern float trace_leaf_brushes( int leaf_id, const LVector3 &start, const LVector3 &end, Trace &trace_out );
extern void clip_box_to_brush( Trace *trace, const LPoint3 &mins, const LPoint3 &maxs,
                               const LPoint3 &p1, const LPoint3 &p2, dbrush_t *brush );

extern void calc_ray_ambient_lighting( int thread, const LVector3 &start,
                                       const LVector3 &end, float tan_theta,
                                       LVector3 *color );

extern void compute_ambient_from_surface( dface_t *face, directlight_t *skylight,
                                          LRGBColor &color );
extern void compute_lightmap_color_from_average( dface_t *face, directlight_t *skylight,
                                                 float scale, LVector3 *colors );
extern void compute_lightmap_color_point_sample( dface_t *face, directlight_t *skylight, LTexCoordf &coord, float scale, LVector3 *colors );

extern void ComputeIndirectLightingAtPoint( const LVector3 &vpos, const LNormalf &vnormal, LVector3 &color, bool ignore_normals );
extern void ComputeDirectLightingAtPoint( const LVector3 &vpos, const LNormalf &vnormal, LVector3 &color );

#endif // LIGHTINGUTILS_H