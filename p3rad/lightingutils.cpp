#include "lightingutils.h"
#include "winding.h"
#include "anorms.h"
#include "halton.h"

#define NEVER_UPDATED -9999

LightSurface::LightSurface( int thread ) :
        _thread( thread ), _surface( nullptr ),
        _hit_frac( 1.0 ), _has_luxel( false )
{
}

bool LightSurface::enumerate_node( int node_id, const Ray &ray,
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

bool LightSurface::enumerate_leaf( int leaf_id, const Ray &ray, float start,
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

bool LightSurface::find_intersection( const Ray &ray )
{
        return !enumerate_nodes_along_ray( ray, this, 0 );
}

bool LightSurface::test_point_against_surface( const LVector3 &point, dface_t *face, texinfo_t *tex )
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

bool LightSurface::test_point_against_sky_surface( const LVector3 &point, dface_t *face )
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

                }
                else
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

        colorrgbexp32_t *lightmap = &g_dlightdata[face->lightofs];
        lightmap += dt * smax + ds;
        for ( int maps = 0; maps < MAXLIGHTMAPS && face->styles[maps] != 0xFF; maps++ )
        {
                int style = face->styles[maps];

                LVector3 color( 0 );
                ColorRGBExp32ToVector( *lightmap, color );
                color /= 255;

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

void ComputeIndirectLightingAtPoint( const LVector3 &vpos, const LNormalf &vnormal, LVector3 &color, bool ignore_normals )
{
        LightSurface surf( 0 );

        int samples = NUMVERTEXNORMALS;
        if ( g_fastmode )
        {
                samples /= 4;
        }

        PN_stdfloat total_dot = 0;
        DirectionalSampler_t sampler;
        for ( int i = 0; i < samples; i++ )
        {
                LVector3 sampling_normal = sampler.NextValue();
                PN_stdfloat dot;
                if ( ignore_normals )
                {
                        dot = 0.7071 / 2;
                }
                else
                {
                        dot = DotProduct( vnormal, sampling_normal );
                }
                if ( dot <= EQUAL_EPSILON )
                {
                        // reject angles behind our plane
                        continue;
                }
                total_dot += dot;

                // trace to determine surface
                LVector3 end;
                VectorScale( sampling_normal, MAX_TRACE_LENGTH, end );
                VectorAdd( vpos, end, end );
                Ray ray( vpos, end, LVector3::zero(), LVector3::zero() );
                if ( !surf.find_intersection( ray ) )
                {
                        continue;
                }

                // get color from surface lightmap
                texinfo_t *tex = &g_texinfo[surf._surface->texinfo];
                if ( !tex || tex->flags & TEX_SPECIAL )
                {
                        // ignore contribution from sky or non lit textures
                        // sky ambient already accounted for during direct pass
                        continue;
                }

                if ( surf._surface->styles[0] == 255 || surf._surface->lightofs == -1 )
                {
                        // no light affects this face
                        continue;
                }

                LVector3 lightmap_col;
                if ( !surf._has_luxel )
                {
                        lightmap_col = dface_AvgLightColor( surf._surface, 0 );
                }
                else
                {
                        // get color from the luxel itself
                        int smax = surf._surface->lightmap_size[0] + 1;
                        int tmax = surf._surface->lightmap_size[1] + 1;

                        // luxelcoord is in the space of the accumulated lightmap page; we need to convert
                        // it to be in the space of the surface
                        int ds = clamp( (int)surf._luxel_coord[0], 0, smax - 1 );
                        int dt = clamp( (int)surf._luxel_coord[1], 0, tmax - 1 );

                        colorrgbexp32_t *lightmap = &g_dlightdata[surf._surface->lightofs];
                        lightmap += dt * smax + ds;
                        ColorRGBExp32ToVector( *lightmap, lightmap_col );
                }

                //VectorMultiply( lightmap_col, 1.0, lightmap_col );
                VectorAdd( color, lightmap_col, color );
        }

        if ( total_dot )
        {
                VectorScale( color, 1.0 / total_dot, color );
        }
}

void ComputeDirectLightingAtPoint( const LVector3 &vpos, const LNormalf &vnormal, LVector3 &color )
{
        vec3_t pos;
        VectorCopy( vpos, pos );
        vec3_t normal;
        VectorCopy( vnormal, normal );

        dleaf_t *leaf = PointInLeaf( pos );
        byte pvs[( MAX_MAP_LEAFS + 7 ) / 8];
        DecompressVis( &g_dvisdata[leaf->visofs], pvs, sizeof( pvs ) );
        byte styles[ALLSTYLES];
        memset( styles, 255, sizeof( styles ) );
        styles[0] = 0;
        vec3_t sample[ALLSTYLES];
        memset( sample, 0, sizeof( sample ) );

        GatherSampleLight( pos, pvs, normal, sample, styles, 0, -1, -1 );

        VectorCopy( sample[0], color );
}