#ifndef LIGHTINGUTILS_H
#define LIGHTINGUTILS_H

#include <aa_luse.h>

#include "bspfile.h"
#include "qrad.h"

#define TEST_EPSILON 0.03125
#define DIST_EPSILON 0.03125

struct Ray
{
        Ray( const LPoint3 &s, const LPoint3 &e,
             const LPoint3 &mi, const LPoint3 &ma )
        {
                start = s;
                end = e;
                mins = mi;
                maxs = ma;
                delta = end - start;

                is_swept = delta.length_squared() != 0;
                extents = maxs - mins;
                is_ray = extents.length_squared() < 1e-6;
                start_offset = ( mins + maxs ) / 2;
                start = ( start + start_offset );
                start_offset /= -1;
        }

        LPoint3 start;
        LPoint3 end;
        LPoint3 mins;
        LPoint3 maxs;
        LVector3 extents;
        bool is_swept;
        bool is_ray;
        LVector3 start_offset;
        LVector3 delta;
};

struct Trace
{
        LPoint3 start_pos;
        LPoint3 end_pos;
        dplane_t plane;
        PN_stdfloat fraction;
        int contents;
        bool all_solid;
        bool start_solid;
        bool is_point;
        LPoint3 mins;
        LPoint3 maxs;
        LVector3 extents;
        texinfo_t *surface;
};

class LightSurface
{
public:
        LightSurface( int thread );

        bool enumerate_node( int node_id, const Ray &ray,
                             float f, int context );

        bool enumerate_leaf( int leaf_id, const Ray &ray, float start,
                             float end, int context );

        bool find_intersection( const Ray &ray );

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

extern bool r_enumerate_nodes_along_ray( int node_id, const Ray &ray, float start,
                                         float end, LightSurface *surf, int context );

extern bool enumerate_nodes_along_ray( const Ray &ray, LightSurface *surf, int context );

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