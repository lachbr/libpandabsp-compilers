#if 1

#include "leaf_ambient_lighting.h"
#include "qrad.h"
#include "anorms.h"

#include <aa_luse.h>
#include <randomizer.h>
#include <plane.h>

#define DIST_EPSILON 0.03125

FORCEINLINE PN_stdfloat inv_r_squared( const LVector3 &v )
{
        return 1.f / std::max( 1.f, v[0] * v[0] + v[1] * v[1] + v[2] * v[2] );
}

static LVector3i box_directions[6] = {
        LVector3i( 1,  0,  0 ),
        LVector3i( -1,  0,  0 ),
        LVector3i( 0,  1,  0 ),
        LVector3i( 0, -1,  0 ),
        LVector3i( 0,  0,  1 ),
        LVector3i( 0,  0, -1 )
};

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

struct AmbientSample
{
        LVector3 pos;
        LVector3 cube[6];
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

#define NEVER_UPDATED -9999

void clip_box_to_brush( Trace *trace, const LPoint3 &mins, const LPoint3 &maxs,
                        const LPoint3 &p1, const LPoint3 &p2, dbrush_t *brush )
{
        dplane_t *plane, *clip_plane;
        float dist;
        LVector3 ofs;
        float d1, d2;
        float f;
        dbrushside_t *side, *leadside;

        if ( !brush->numsides )
        {
                return;
        }

        float enter_frac = NEVER_UPDATED;
        float leave_frac = 1.0;
        clip_plane = nullptr;

        bool getout = false;
        bool startout = false;
        leadside = nullptr;

        for ( int i = 0; i < brush->numsides; i++ )
        {
                side = &g_dbrushsides[brush->firstside + i];
                plane = &g_dplanes[side->planenum];

                if ( !trace->is_point )
                {
                        ofs[0] = ( plane->normal[0] < 0 ) ? maxs[0] : mins[0];
                        ofs[1] = ( plane->normal[1] < 0 ) ? maxs[1] : mins[1];
                        ofs[2] = ( plane->normal[2] < 0 ) ? maxs[2] : mins[2];
                        dist = DotProduct( ofs, plane->normal );
                        dist = plane->dist - dist;

                }
                else
                {
                        if ( side->bevel == 1 )
                        {
                                continue;
                        }

                        dist = plane->dist;
                }

                d1 = DotProduct( p1, plane->normal ) - dist;
                d2 = DotProduct( p2, plane->normal ) - dist;

                if ( d1 > 0 && d2 > 0 )
                {
                        return;
                }

                if ( d2 > 0 )
                {
                        getout = true;
                }
                if ( d1 > 0 )
                {
                        startout = true;
                }

                if ( d1 <= 0 && d2 <= 0 )
                {
                        continue;
                }

                if ( d1 > d2 )
                {
                        f = ( d1 - DIST_EPSILON ) / ( d1 - d2 );
                        if ( f > enter_frac )
                        {
                                enter_frac = f;
                                clip_plane = plane;
                                leadside = side;
                        }

                }
                else
                {
                        f = ( d1 + DIST_EPSILON ) / ( d1 - d2 );
                        if ( f < leave_frac )
                        {
                                leave_frac = f;
                        }
                }
        }

        if ( !startout )
        {
                trace->start_solid = true;
                if ( !getout )
                {
                        trace->all_solid = true;
                }

                return;
        }

        if ( enter_frac < leave_frac )
        {
                if ( enter_frac > NEVER_UPDATED && enter_frac < trace->fraction )
                {
                        if ( enter_frac < 0 )
                        {
                                enter_frac = 0;
                        }

                        trace->fraction = enter_frac;
                        trace->plane.dist = clip_plane->dist;
                        VectorCopy( clip_plane->normal, trace->plane.normal );
                        trace->plane.type = clip_plane->type;
                        if ( leadside->texinfo != -1 )
                        {
                                trace->surface = &g_texinfo[leadside->texinfo];
                        }
                        else
                        {
                                trace->surface = 0;
                        }
                        trace->contents = brush->contents;
                }
        }
}

float trace_leaf_brushes( int leaf_id, const LVector3 &start, const LVector3 &end, Trace &trace_out )
{
        dleaf_t *leaf = &g_dleafs[leaf_id];
        Trace trace;
        memset( &trace, 0, sizeof( Trace ) );
        trace.is_point = true;
        trace.start_solid = false;
        trace.fraction = 1.0;

        for ( int i = 0; i < leaf->numleafbrushes; i++ )
        {
                int brushnum = g_dleafbrushes[leaf->firstleafbrush + i];
                dbrush_t *b = &g_dbrushes[brushnum];
                if ( b->contents != CONTENTS_SOLID )
                {
                        continue;
                }
                LVector3 extents( 0 );
                clip_box_to_brush( &trace, extents, extents, start, end, b );
                if ( trace.fraction != 1.0 || trace.start_solid )
                {
                        if ( trace.start_solid )
                        {
                                trace.fraction = 0.0;
                        }
                        trace_out = trace;
                        return trace.fraction;
                }
        }

        trace_out = trace;
        return 1.0;
}

typedef pvector<AmbientSample> vector_ambientsample;
typedef pvector<dplane_t> vector_dplane;
pvector<vector_ambientsample> leaf_ambient_samples;

class LeafSampler
{
public:
        LeafSampler( int thread ) :
                _thread( thread )
        {
        }

        bool cast_ray_in_leaf( const LVector3 &start, const LVector3 &end,
                               int leaf_id, float *fraction, LNormalf &normal )
        {
                fraction[0] = 1.0;

                Trace trace;
                if ( trace_leaf_brushes( leaf_id, start, end, trace ) != 1.0 )
                {
                        fraction[0] = trace.fraction;
                        normal[0] = trace.plane.normal[0];
                        normal[1] = trace.plane.normal[1];
                        normal[2] = trace.plane.normal[2];
                }
                else
                {
                        nassertr( !trace.start_solid && !trace.all_solid, false );
                }

                return fraction[0] != 1.0;
        }

        void generate_leaf_sample_position( int leaf_id, const vector_dplane &leaf_planes,
                                            LVector3 &sample_pos )
        {
                dleaf_t *leaf = &g_dleafs[leaf_id];

                float dx = leaf->maxs[0] - leaf->mins[0];
                float dy = leaf->maxs[1] - leaf->mins[1];
                float dz = leaf->maxs[2] - leaf->mins[2];

                bool valid = false;
                for ( int i = 0; i < 1000 && !valid; i++ )
                {
                        sample_pos[0] = leaf->mins[0] + _random.random_real( dx );
                        sample_pos[1] = leaf->mins[1] + _random.random_real( dy );
                        sample_pos[2] = leaf->mins[2] + _random.random_real( dz );
                        valid = true;
                        for ( int j = (int)leaf_planes.size() - 1; j >= 0 && valid; j-- )
                        {
                                float dist = DotProduct( leaf_planes[j].normal, sample_pos ) - leaf_planes[j].dist;
                                if ( dist < DIST_EPSILON )
                                {
                                        // Not inside the leaf, try again.
                                        valid = false;
                                        break;
                                }
                        }

                        if ( !valid )
                        {
                                continue;
                        }

                        for ( int j = 0; j < 6; j++ )
                        {
                                LVector3 start = sample_pos;
                                int axis = j % 3;
                                start[axis] = ( j < 3 ) ? leaf->mins[axis] : leaf->maxs[axis];
                                float t;
                                LNormalf normal;
                                cast_ray_in_leaf( sample_pos, start, leaf_id, &t, normal );
                                if ( t == 0.0 )
                                {
                                        valid = false;
                                        break;
                                }

                                if ( t != 1.0 )
                                {
                                        LVector3 delta = start - sample_pos;
                                        if ( delta.dot( normal ) > 0 )
                                        {
                                                valid = false;
                                                break;
                                        }
                                }
                        }
                }

                if ( !valid )
                {
                        // Didn't generate a valid sample point, just use the center of the leaf bbox
                        sample_pos = ( LPoint3( leaf->mins[0], leaf->mins[1], leaf->mins[2] ) + LPoint3( leaf->maxs[0], leaf->maxs[1], leaf->maxs[2] ) ) / 2.0;
                }
        }

private:
        int _thread;
        Randomizer _random;
};

bool is_leaf_ambient_surface_light( directlight_t *dl )
{
        static const float directlight_min_emit_surface = 0.005f;
        static const float directlight_min_emit_surface_distance_ratio = inv_r_squared( LVector3( 0, 0, 512 ) );

        if ( dl->type != emit_surface )
        {
                return false;
        }

        if ( dl->style != 0 )
        {
                return false;
        }

        PN_stdfloat intensity = std::max( dl->intensity[0], dl->intensity[1] );
        intensity = std::max( intensity, dl->intensity[2] );

        return ( intensity * directlight_min_emit_surface_distance_ratio ) < directlight_min_emit_surface;
}

void get_leaf_boundary_planes( vector_dplane &list, int leaf_id )
{
        list.clear();
        int node_id = leafparents[leaf_id];
        int child = ~leaf_id;

        while ( node_id >= 0 )
        {
                dnode_t *node = &g_dnodes[node_id];
                dplane_t *node_plane = &g_dplanes[node->planenum];
                if ( node->children[0] == child )
                {
                        // Front side
                        list.push_back( *node_plane );
                }
                else
                {
                        dplane_t plane;
                        plane.dist = -node_plane->dist;
                        plane.normal[0] = -node_plane->normal[0];
                        plane.normal[1] = -node_plane->normal[1];
                        plane.normal[2] = -node_plane->normal[2];
                        plane.type = node_plane->type;
                }
                child = node_id;
                node_id = nodeparents[child];
        }
}

static directlight_t *find_ambient_sky_light()
{
        static directlight_t *found_light = nullptr;

        // So we don't have to keep finding the same sky light
        if ( found_light == nullptr )
        {
                for ( int i = 0; i < numdlights; i++ )
                {
                        directlight_t *dl = directlights[i];
                        if ( dl == nullptr )
                        {
                                continue;
                        }
                        if ( dl->type == emit_skylight )
                        {
                                found_light = dl;
                                break;
                        }
                }
        }

        return found_light;
}

class LightSurface;
bool enumerate_nodes_along_ray( const Ray &ray, LightSurface *surf, int context );

class LightSurface
{
public:
        LightSurface( int thread ) :
                _thread( thread ), _surface( nullptr ),
                _hit_frac( 1.0 ), _has_luxel( false )
        {
        }

        bool enumerate_node( int node_id, const Ray &ray,
                             float f, int context )
        {
                dface_t *sky_surface = nullptr;

                LVector3 pt;
                VectorMA( ray.start, f, ray.delta, pt );

                dnode_t *node = &g_dnodes[node_id];
                
                for ( int i = 0; i < node->numfaces; i++ )
                {
                        dface_t *face = &g_dfaces[node->firstface + i];
                        // Don't take into account faces that are in a leaf
                        if ( face->on_node == 0 )
                        {
                                continue;
                        }

                        // TODO: Don't test displacement faces

                        texinfo_t *tex = &g_texinfo[face->texinfo];
                        if ( tex->flags & TEX_SPECIAL )
                        {
                                if ( test_point_against_sky_surface( pt, face ) )
                                {
                                        sky_surface = face;
                                }

                                continue;
                        }

                        if ( test_point_against_surface( pt, face, tex ) )
                        {
                                _hit_frac = f;
                                _surface = face;
                                _has_luxel = true;
                                return false;
                        }
                }

                // if we hit a sky surface, return it
                _surface = sky_surface;
                return ( _surface == nullptr );
        }

        bool enumerate_leaf( int leaf_id, const Ray &ray, float start,
                             float end, int context )
        {
                bool hit = false;
                dleaf_t *leaf = &g_dleafs[leaf_id];
                for ( int i = 0; i < leaf->nummarksurfaces; i++ )
                {
                        dface_t *face = &g_dfaces[leaf->firstmarksurface + i];

                        if ( face->on_node == 1 )
                        {
                                continue;
                        }

                        texinfo_t *tex = &g_texinfo[face->texinfo];
                        dplane_t *plane = &g_dplanes[face->planenum];

                        // Backface cull
                        if ( DotProduct( plane->normal, ray.delta ) > 0 )
                        {
                                continue;
                        }

                        float start_dot_n = DotProduct( ray.start, plane->normal );
                        float delta_dot_n = DotProduct( ray.delta, plane->normal );

                        float front = start_dot_n + start * delta_dot_n - plane->dist;
                        float back = start_dot_n + end * delta_dot_n - plane->dist;

                        int side = front < 0;
                        // Blow it off if it doesn't split the plane
                        if ( (int)( back < 0 ) == side )
                        {
                                continue;
                        }

                        // Dont't test a surface that is further away from the closest found intersection.
                        float f = front / ( front - back );
                        float mid = start * ( 1.0 - f ) + end * f;
                        if ( mid >= _hit_frac )
                        {
                                continue;
                        }

                        LVector3 pt;
                        VectorMA( ray.start, mid, ray.delta, pt );

                        if ( test_point_against_surface( pt, face, tex ) )
                        {
                                _hit_frac = mid;
                                _surface = face;
                                hit = true;
                                _has_luxel = true;
                        }  
                }

                return !hit;
        }

        bool find_intersection( const Ray &ray )
        {
                return !enumerate_nodes_along_ray( ray, this, 0 );
        }

private:
        bool test_point_against_surface( const LVector3 &point, dface_t *face, texinfo_t *tex )
        {
                // Specials don't have lightmaps
                if ( tex->flags & TEX_SPECIAL )
                {
                        return false;
                }

                float s, t;
                s = DotProduct( point, tex->lightmap_vecs[0] ) + tex->lightmap_vecs[0][3];
                t = DotProduct( point, tex->lightmap_vecs[1] ) + tex->lightmap_vecs[1][3];

                if ( s < face->lightmap_mins[0] || t < face->lightmap_mins[1] )
                {
                        // Not in bounds of lightmap
                        return false;
                }

                float ds = s - face->lightmap_mins[0];
                float dt = t - face->lightmap_mins[1];
                if ( ds < face->lightmap_size[0] || dt < face->lightmap_size[1] )
                {
                        // Doesn't lie in rectangle.
                        return false;
                }

                _luxel_coord.set( ds, dt );
                return true;
        }

        bool test_point_against_sky_surface( const LVector3 &point, dface_t *face )
        {
                Winding winding( *face );

                dplane_t plane;
                winding.getPlane( plane );

                vec3_t vpoint;
                vpoint[0] = point[0];
                vpoint[1] = point[1];
                vpoint[2] = point[2];
                return point_in_winding( winding, plane, vpoint );
        }

public:
        int _thread;
        dface_t *_surface;
        float _hit_frac;
        LTexCoordf _luxel_coord;
        bool _has_luxel;
};

#define TEST_EPSILON	(0.03125)

bool r_enumerate_nodes_along_ray( int node_id, const Ray &ray, float start,
                                  float end, LightSurface *surf, int context )
{
        float front, back;
        float start_dot_n, delta_dot_n;

        while ( node_id >= 0 )
        {
                dnode_t *node = &g_dnodes[node_id];
                dplane_t *plane = &g_dplanes[node->planenum];

                if ( plane->type == plane_z )
                {
                        start_dot_n = ray.start[plane->type];
                        delta_dot_n = ray.delta[plane->type];

                }
                else
                {
                        start_dot_n = DotProduct( ray.start, plane->normal );
                        delta_dot_n = DotProduct( ray.delta, plane->normal );
                }

                front = start_dot_n + start * delta_dot_n - plane->dist;
                back = start_dot_n + end * delta_dot_n - plane->dist;

                if ( front <= -TEST_EPSILON && back <= -TEST_EPSILON )
                {
                        node_id = node->children[1];

                }
                else if ( front >= TEST_EPSILON && back >= TEST_EPSILON )
                {
                        node_id = node->children[0];

                }
                else
                {
                        // test the front side first
                        bool side = front < 0;

                        float split_frac;
                        if ( delta_dot_n == 0.0 )
                        {
                                split_frac = 1.0;

                        }
                        else
                        {
                                split_frac = ( plane->dist - start_dot_n ) / delta_dot_n;
                                if ( split_frac < 0.0 )
                                {
                                        split_frac = 0.0;
                                }
                                else if ( split_frac > 1.0 )
                                {
                                        split_frac = 1.0;
                                }
                        }

                        bool r = r_enumerate_nodes_along_ray( node->children[side], ray, start,
                                                              split_frac, surf, context );
                        if ( !r )
                        {
                                return r;
                        }

                        // Visit the node...
                        if ( !surf->enumerate_node( node_id, ray, split_frac, context ) )
                        {
                                return false;
                        }

                        return r_enumerate_nodes_along_ray( node->children[!side], ray, split_frac,
                                                            end, surf, context );
                }
        }

        // Visit the leaf...
        return surf->enumerate_leaf( ~node_id, ray, start, end, context );
}

bool enumerate_nodes_along_ray( const Ray &ray, LightSurface *surf, int context )
{
        return r_enumerate_nodes_along_ray( 0, ray, 0.0, 1.0, surf, context );
}

void compute_ambient_from_surface( dface_t *face, directlight_t *skylight,
                                   LRGBColor &color )
{
        texinfo_t *texinfo = &g_texinfo[face->texinfo];
        if ( texinfo )
        {
                if ( texinfo->flags & TEX_SPECIAL )
                {
                        if ( skylight )
                        {
                                // Add in sky ambient
                                vec3_t linear;
                                linear[0] = linear[1] = linear[2] = 255.0f;
                                VectorDivide( skylight->diffuse_intensity, linear, color );
                        }

                } else
                {
                        vec3_t one;
                        one[0] = one[1] = one[2] = 1.0f;
                        VectorMultiply( color, one, color );
                }
        }
}

static void compute_lightmap_color_from_average( dface_t *face, directlight_t *skylight,
                                                 float scale, LVector3 *colors )
{
        texinfo_t *tex = &g_texinfo[face->texinfo];
        if ( tex->flags & TEX_SPECIAL )
        {
                if ( skylight )
                {
                        LVector3 amb( skylight->diffuse_intensity[0] / 255.0,
                                      skylight->diffuse_intensity[1] / 255.0,
                                      skylight->diffuse_intensity[2] / 255.0 );
                        colors[0] += amb * scale;
                }
                return;
        }

        for ( int maps = 0; maps < MAXLIGHTMAPS && face->styles[maps] != 0xFF; maps++ )
        {
                LRGBColor avg_color = dface_AvgLightColor( face, maps );
                LRGBColor color = avg_color / 255.0f;

                compute_ambient_from_surface( face, skylight, color );

                int style = face->styles[maps];
                colors[style] += color * scale;
        }
}

void compute_lightmap_color_point_sample( dface_t *face, directlight_t *skylight, LTexCoordf &coord, float scale, LVector3 *colors )
{
        // Face unaffected by light
        if ( face->lightofs == -1 )
        {
                return;
        }

        int smax = face->lightmap_size[0] + 1;
        int tmax = face->lightmap_size[1] + 1;

        int ds = clamp( (int)coord[0], 0, smax - 1 );
        int dt = clamp( (int)coord[1], 0, tmax - 1 );

        int offset = smax * tmax;
        // bumped lightmaps todo!

        unsigned char *lightmap = (unsigned char *)&g_dlightdata[face->lightofs];
        lightmap += dt * smax + ds;
        for ( int maps = 0; maps < MAXLIGHTMAPS && face->styles[maps] != 0xFF; maps++ )
        {
                int style = face->styles[maps];

                LRGBColor color;
                color[0] = ( lightmap[0] / 255.0 );
                color[1] = ( lightmap[1] / 255.0 );
                color[2] = ( lightmap[2] / 255.0 );

                compute_ambient_from_surface( face, skylight, color );
                colors[style] += color * scale;

                lightmap += offset;
        }
}

void calc_ray_ambient_lighting( int thread, const LVector3 &start,
                                const LVector3 &end, float tan_theta,
                                LVector3 *color )
{
        directlight_t *skylight = find_ambient_sky_light();

        LightSurface surf( thread );
        Ray ray( start, end, LVector3::zero(), LVector3::zero() );
        if ( !surf.find_intersection( ray ) )
        {
                return;
        }

        // compute the approximate radius of a circle centered around the intersection point
        float dist = ray.delta.length() * tan_theta * surf._hit_frac;

        // until 20" we use the point sample, then blend in the average until we're covering 40"
        // This is attempting to model the ray as a cone - in the ideal case we'd simply sample all
        // luxels in the intersection of the cone with the surface.  Since we don't have surface 
        // neighbor information computed we'll just approximate that sampling with a blend between
        // a point sample and the face average.
        // This yields results that are similar in that aliasing is reduced at distance while 
        // point samples provide accuracy for intersections with near geometry
        float scale_avg = RemapValClamped( dist, 20, 40, 0.0f, 1.0f );

        if ( !surf._has_luxel )
        {
                scale_avg = 1.0;
        }

        float scale_sample = 1.0f - scale_avg;
        if ( scale_avg != 0 )
        {
                compute_lightmap_color_from_average( surf._surface, skylight, scale_avg, color );
        }
        if ( scale_sample != 0 )
        {
                compute_lightmap_color_point_sample( surf._surface, skylight, surf._luxel_coord, scale_sample, color );
        }
}

float worldlight_angle( const directlight_t *dl, const LVector3 &lnormal,
                        const LVector3 &snormal, const LVector3 &delta )
{
        float dot, dot2;
        if ( dl->type != emit_surface )
        {
                return 1.0;
        }

        dot = DotProduct( snormal, delta );
        if ( dot < 0 )
        {
                return 0;
        }

        dot2 = -DotProduct( delta, lnormal );
        if ( dot2 <= ON_EPSILON / 10 )
        {
                return 0; // behind light surface
        }

        return dot * dot2;
}

float worldlight_distance_falloff( const directlight_t *dl, const vec3_t &delta )
{
        if ( dl->type != emit_surface )
        {
                return 1.0;
        }

        LVector3 vec;
        VectorCopy( delta, vec );
        return inv_r_squared( vec );
}

void add_emit_surface_lights( const LVector3 &start, LVector3 *cube )
{
        float fraction_visible[4];
        memset( fraction_visible, 0, sizeof( fraction_visible ) );

        vec3_t vstart;
        VectorCopy( start, vstart );

        for ( int i = 0; i < numdlights; i++ )
        {
                directlight_t *dl = directlights[i];

                if ( dl == nullptr )
                {
                        continue;
                }

                if ( !( dl->flags & DLF_in_ambient_cube ) )
                {
                        continue;
                }

                if ( dl->type != emit_surface )
                {
                        continue;
                }

                if ( TestLine( vstart, dl->origin ) == CONTENTS_SOLID )
                {
                        continue;
                }

                // Add this light's contribution.
                vec3_t delta;
                VectorSubtract( dl->origin, vstart, delta );
                float distance_scale = worldlight_distance_falloff( dl, delta );

                LVector3 deltanorm;
                VectorCopy( delta, deltanorm );
                deltanorm.normalize();
                LVector3 dlnorm;
                VectorCopy( dl->normal, dlnorm );
                float angle_scale = worldlight_angle( dl, dlnorm, deltanorm, deltanorm );

                float ratio = distance_scale * angle_scale * 1.0;
                if ( ratio == 0 )
                {
                        continue;
                }

                for ( int i = 0; i < 6; i++ )
                {
                        float t = DotProduct( box_directions[i], deltanorm );
                        if ( t > 0 )
                        {
                                LVector3 intensity;
                                VectorCopy( dl->intensity, intensity );
                                intensity *= ( t * ratio );
                                VectorAdd( cube[i], intensity, cube[i] );
                        }
                }
        }
}

#define MAX_COORD_INTEGER (16384)
#define COORD_EXTENT (2 * MAX_COORD_INTEGER)

void compute_ambient_from_spherical_samples( int thread, const LVector3 &sample_pos,
                                             LVector3 *cube )
{
        LVector3 radcolor[NUMVERTEXNORMALS];
        float tan_theta = std::tan( VERTEXNORMAL_CONE_INNER_ANGLE );

        for ( int i = 0; i < NUMVERTEXNORMALS; i++ )
        {
                LVector3 end = sample_pos + g_anorms[i] * ( COORD_EXTENT * 1.74 );

                LVector3 light_style_colors[MAX_LIGHTSTYLES];
                light_style_colors[0].set( 0, 0, 0 );
                calc_ray_ambient_lighting( thread, sample_pos, end, tan_theta, light_style_colors );
                radcolor[i] = light_style_colors[0];
        }

        // accumulate samples into radiant box
        for ( int j = 5; j >= 0; j-- )
        {
                float t = 0;
                cube[j].set( 0, 0, 0 );
                for ( int i = 0; i < NUMVERTEXNORMALS; i++ )
                {
                        float c = DotProduct( g_anorms[i], box_directions[j] );
                        if ( c > 0 )
                        {
                                t += c;
                                cube[j] += radcolor[i] * c;
                        }
                }

                cube[j] *= ( 1 / t );
        }

        // Now add direct light from the emit_surface lights. These go in the ambient cube because
        // there are a ton of them and they are often so dim that they get filtered out by r_worldlightmin.
        add_emit_surface_lights( sample_pos, cube );
}

void add_sample_to_list( vector_ambientsample &list, const LVector3 &sample_pos, LVector3 *cube )
{
        const size_t max_samples = 16;

        AmbientSample sample;
        sample.pos = sample_pos;
        for ( int i = 0; i < 6; i++ )
        {
                sample.cube[i] = cube[i];
        }

        list.push_back( sample );

        if ( list.size() <= max_samples )
        {
                return;
        }

        int nearest_neighbor_idx = 0;
        float nearest_neighbor_dist = FLT_MAX;
        float nearest_neighbor_total = 0;
        for ( int i = 0; i < list.size(); i++ )
        {
                int closest_idx = 0;
                float closest_dist = FLT_MAX;
                float total_dc = 0;
                for ( int j = 0; j < list.size(); j++ )
                {
                        if ( j == i )
                        {
                                continue;
                        }

                        float dist = ( list[i].pos - list[j].pos ).length();
                        float max_dc = 0;
                        for ( int k = 0; k < 6; k++ )
                        {
                                // color delta is computed per-component, per cube side
                                for ( int s = 0; s < 3; s++ )
                                {
                                        float dc = std::fabs( list[i].cube[k][s] - list[j].cube[k][s] );
                                        max_dc = std::max( max_dc, dc );
                                }
                                total_dc += max_dc;
                        }

                        // need a measurable difference in color or we'll just rely on position
                        if ( max_dc < 1e-4f )
                        {
                                max_dc = 0;

                        }
                        else if ( max_dc > 1.0f )
                        {
                                max_dc = 1.0f;
                        }

                        float distance_factor = 0.1f + ( max_dc * 0.9f );
                        dist *= distance_factor;

                        // find the "closest" sample to this one
                        if ( dist < closest_dist )
                        {
                                closest_dist = dist;
                                closest_idx = j;
                        }
                }

                // the sample with the "closest" neighbor is rejected
                if ( closest_dist < nearest_neighbor_dist || ( closest_dist == nearest_neighbor_dist && total_dc < nearest_neighbor_total ) )
                {
                        nearest_neighbor_dist = closest_dist;
                        nearest_neighbor_idx = i;
                }
        }

        list.erase( list.begin() + nearest_neighbor_idx );
}

void compute_ambient_for_leaf( int thread, int leaf_id,
                               vector_ambientsample &list )
{
        if ( g_dleafs[leaf_id].contents == CONTENTS_SOLID )
        {
                // Don't generate any samples in solid leaves.
                // NOTE: We copy the nearest non-solid leaf sample pointers into this leaf at the end.
                std::cout << "Leaf " << leaf_id << " is solid" << std::endl;
                return;
        }

        vector_dplane leaf_planes;
        LeafSampler sampler( thread );
        get_leaf_boundary_planes( leaf_planes, leaf_id );
        list.clear();

        int xsize = ( g_dleafs[leaf_id].maxs[0] - g_dleafs[leaf_id].mins[0] ) / 32;
        int ysize = ( g_dleafs[leaf_id].maxs[1] - g_dleafs[leaf_id].mins[1] ) / 32;
        int zsize = ( g_dleafs[leaf_id].maxs[2] - g_dleafs[leaf_id].mins[2] ) / 64;
        xsize = std::max( xsize, 1 );
        ysize = std::max( ysize, 1 );
        zsize = std::max( zsize, 1 );

        int volume_count = xsize * ysize * zsize;
        // Don't do any more than 128 samples
        int sample_count = clamp( volume_count, 1, 128 );
        LVector3 cube[6];
        for ( int i = 0; i < 6; i++ )
        {
                LVector3 sample_pos;
                sampler.generate_leaf_sample_position( leaf_id, leaf_planes, sample_pos );
                compute_ambient_from_spherical_samples( thread, sample_pos, cube );
                add_sample_to_list( list, sample_pos, cube );
        }
}

static void thread_compute_leaf_ambient( int thread )
{
        vector_ambientsample list;
        while ( true )
        {
                int leaf_id = GetThreadWork();
                if ( leaf_id == -1 )
                {
                        break;
                }

                list.clear();
                compute_ambient_for_leaf( thread, leaf_id, list );

                leaf_ambient_samples[leaf_id].resize( list.size() );
                for ( size_t i = 0; i < list.size(); i++ )
                {
                        leaf_ambient_samples[leaf_id][i] = list[i];
                }
        }
}

// maps a float to a byte fraction between min & max
static byte fixed_8_fraction( float t, float tMin, float tMax )
{
        if ( tMax <= tMin )
                return 0;

        float frac = RemapValClamped( t, tMin, tMax, 0.0f, 255.0f );
        return byte( frac + 0.5f );
}

void LeafAmbientLighting::
compute_per_leaf_ambient_lighting()
{
        Log( "Computing per leaf ambient lighting...\n" );

        // Figure out which ambient lights should go in the per-leaf ambient cubes.
        int in_ambient_cube = 0;
        int surface_lights = 0;

        for ( int i = 0; i < numdlights; i++ )
        {
                directlight_t *dl = directlights[i];

                if ( dl == nullptr )
                {
                        continue;
                }

                if ( is_leaf_ambient_surface_light( dl ) )
                {
                        dl->flags |= DLF_in_ambient_cube;

                }
                else
                {
                        dl->flags &= ~DLF_in_ambient_cube;
                }

                if ( dl->type == emit_surface )
                {
                        surface_lights++;
                }

                if ( dl->flags & DLF_in_ambient_cube )
                {
                        in_ambient_cube++;
                }
        }

        leaf_ambient_samples.resize( g_numleafs );

        RunThreadsOn( g_numleafs, true, thread_compute_leaf_ambient );

        // now write out the data :)
        g_leafambientindex.clear();
        g_leafambientlighting.clear();
        g_leafambientindex.resize( g_numleafs );
        g_leafambientlighting.reserve( g_numleafs * 4 );
        for ( int leaf_id = 0; leaf_id < g_numleafs; leaf_id++ )
        {
                const vector_ambientsample &list = leaf_ambient_samples[leaf_id];
                g_leafambientindex[leaf_id].num_ambient_samples = list.size();

                if ( list.size() == 0 )
                {
                        g_leafambientindex[leaf_id].first_ambient_sample = 0;
                }
                else
                {
                        g_leafambientindex[leaf_id].first_ambient_sample = g_leafambientlighting.size();

                        for ( int i = 0; i < list.size(); i++ )
                        {
                                dleafambientlighting_t light;
                                light.x = fixed_8_fraction( list[i].pos[0], g_dleafs[leaf_id].mins[0], g_dleafs[leaf_id].maxs[0] );
                                light.y = fixed_8_fraction( list[i].pos[1], g_dleafs[leaf_id].mins[1], g_dleafs[leaf_id].maxs[1] );
                                light.z = fixed_8_fraction( list[i].pos[2], g_dleafs[leaf_id].mins[2], g_dleafs[leaf_id].maxs[2] );
                                light.pad = 0;
                                for ( int side = 0; side < 6; side++ )
                                {
                                        LRGBColor temp;
                                        VectorCopy( list[i].cube[side], temp );
                                        VectorToColorRGBExp32( list[i].cube[side], light.cube.color[side] );
                                }

                                g_leafambientlighting.push_back( light );
                        }
                }
        }

        
        for ( int i = 0; i < g_numleafs; i++ )
        {
                if ( g_leafambientindex[i].num_ambient_samples == 0 )
                {
                        if ( g_dleafs[i].contents == CONTENTS_SOLID )
                        {
                                Warning( "Bad leaf ambient for leaf %d\n", i );
                        }

                        //int ret_leaf = nearest_neighbor_with_light( i );
                        //g_leafambientindex[i].num_ambient_samples = 0;
                        //g_leafambientindex[i].first_ambient_sample = ret_leaf;
                }
        }
        
}

#endif