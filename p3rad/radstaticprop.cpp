#include "radstaticprop.h"
#include <loader.h>
#include "bspfile.h"
#include <collisionNode.h>
#include <geomNode.h>
#include <geom.h>
#include <geomVertexReader.h>
#include <geomVertexData.h>
#include <geomPrimitive.h>
#include <nodePathCollection.h>
#include <lightMutexHolder.h>
#include <look_at.h>
#include <virtualFileSystem.h>

#include "threads.h"
#include "qrad.h"
#include "lightingutils.h"
#include "lightmap.h"

TypeHandle RADCollisionPolygon::_type_handle;

pvector<TestGroup> g_test_groups;
pvector<RADStaticProp *> g_static_props;

bool g_collisions_loaded = false;
LightMutex g_prop_lock( "RadStaticPropRayCastLock" );

void LoadStaticProps()
{
        if ( !g_collisions_loaded )
        {
                for ( int thread = 0; thread < g_numthreads; thread++ )
                {
                        PT( CollisionSegment ) seg = new CollisionSegment;

                        TestGroup tg;
                        tg.seg = seg;

                        g_test_groups.push_back( tg );
                }


                g_collisions_loaded = true;
        }

        Loader *loader = Loader::get_global_ptr();

        for ( size_t i = 0; i < g_bspdata->dstaticprops.size(); i++ )
        {
                dstaticprop_t *prop = &g_bspdata->dstaticprops[i];
                string mdl_path = prop->name;

                NodePath propnp( loader->load_sync( Filename( mdl_path ) ) );
                if ( !propnp.is_empty() )
                {
                        propnp.set_scale( prop->scale[0] * PANDA_TO_HAMMER, prop->scale[1] * PANDA_TO_HAMMER, prop->scale[2] * PANDA_TO_HAMMER );
                        propnp.set_pos( prop->pos[0], prop->pos[1], prop->pos[2] );
                        propnp.set_hpr( prop->hpr[1] - 90, prop->hpr[0], prop->hpr[2] );

                        // Which leaf does the prop reside in?
                        dleaf_t *leaf = PointInLeaf( prop->pos );

                        bool shadow_caster = prop->flags & STATICPROPFLAGS_LIGHTMAPSHADOWS;

                        RADStaticProp *sprop = new RADStaticProp;
                        sprop->shadows = shadow_caster;
                        sprop->leafnum = leaf - g_bspdata->dleafs;
                        sprop->mdl = propnp;
                        sprop->propnum = (int)i;

                        // Apply all transforms and attribs to the vertices so they are now in world space,
                        // but do not remove any nodes or rearrange the vertices.
                        propnp.clear_model_nodes();
                        propnp.flatten_light();

                        NodePathCollection npc = propnp.find_all_matches( "**/+GeomNode" );
                        for ( int i = 0; i < npc.get_num_paths(); i++ )
                        {
                                NodePath geomnp = npc.get_path( i );
                                LMatrix4 mat_to_world = geomnp.get_net_transform()->get_mat();
                                GeomNode *gn = DCAST( GeomNode, geomnp.node() );
                                for ( int j = 0; j < gn->get_num_geoms(); j++ )
                                {
                                        PT( Geom ) geom = gn->get_geom( j )->decompose();
                                        const GeomVertexData *vdata = geom->get_vertex_data();
                                        GeomVertexReader vreader( vdata, InternalName::get_vertex() );
                                        std::stringstream ss;
                                        ss << gn->get_name() << "-" << j;
                                        PT( CollisionNode ) cnode = new CollisionNode( ss.str() );
                                        for ( int k = 0; k < geom->get_num_primitives(); k++ )
                                        {
                                                const GeomPrimitive *prim = geom->get_primitive( k );
                                                for ( int l = 0; l < prim->get_num_primitives(); l++ )
                                                {
                                                        int start = prim->get_primitive_start( l );
                                                        int end = prim->get_primitive_end( l );

                                                        pvector<LPoint3> verts;
                                                        for ( int m = start; m < end; m++ )
                                                        {
                                                                vreader.set_row( prim->get_vertex( m ) );
                                                                verts.push_back( vreader.get_data3f() );
                                                        }

                                                        PT( RADCollisionPolygon ) poly = nullptr;
                                                        if ( verts.size() == 3 )
                                                                poly = new RADCollisionPolygon( verts[0], verts[1], verts[2] );
                                                        else if ( verts.size() == 4 )
                                                                poly = new RADCollisionPolygon( verts[0], verts[1], verts[2], verts[3] );
                                                        if ( poly != nullptr )
                                                                sprop->polygons.push_back( poly );
                                                }
                                        }
                                }
                        }

                        g_static_props.push_back( sprop );

                        cout << "Successfully loaded static prop " << mdl_path << endl;
                }
                else
                {
                        cout << "Warning! Could not load static prop " << mdl_path << ", no shadows" << endl;
                }
        }
}

bool StaticPropIntersectionTest( const vec3_t start, const vec3_t stop, int leafnum )
{
        //LightMutexHolder holder( g_prop_lock );

        if ( start[0] == stop[0] && start[1] == stop[1] && start[2] == stop[2] )
        {
                return false;
        }

        int thread = GetCurrentThreadNumber();
        TestGroup *tg = &g_test_groups[thread];

        tg->seg->set_point_a( start[0], start[1], start[2] );
        tg->seg->set_point_b( stop[0], stop[1], stop[2] );

        for ( size_t i = 0; i < g_static_props.size(); i++ )
        {
                RADStaticProp *sprop = g_static_props[i];
                if ( !sprop->shadows || sprop->leafnum != leafnum )
                        continue;
                for ( size_t polynum = 0; polynum < sprop->polygons.size(); polynum++ )
                {
                        RADCollisionPolygon *poly = sprop->polygons[polynum];
                        if ( poly->segment_intersection_test( start, stop ) )
                        {
                                //cout << "Collided with polygon in leaf " << leafnum << endl;
                                return true;
                        }

                }
        }

        return false;
}

void BuildGeomNodes_r( PandaNode *node, pvector<PT( GeomNode )> &list )
{
        if ( node->is_of_type( GeomNode::get_class_type() ) )
        {
                list.push_back( DCAST( GeomNode, node ) );
        }

        for ( int i = 0; i < node->get_num_children(); i++ )
        {
                BuildGeomNodes_r( node->get_child( i ), list );
        }
}

pvector<PT( GeomNode )> BuildGeomNodes( const NodePath &root )
{
        pvector<PT( GeomNode )> list;

        BuildGeomNodes_r( root.node(), list );

        return list;
}

struct VDataDef
{
        CPT( GeomVertexData ) vdata;
        LMatrix4 mat_to_world;
};

void ComputeStaticPropLighting( int thread )
{
        int prop_idx = thread;//GetThreadWork();
        if ( prop_idx == -1 )
        {
                return;
        }
        RADStaticProp *prop = g_static_props[prop_idx];
        if ( prop == nullptr )
        {
                Warning( "ThreadComputeStaticPropLighting: prop is nullptr on thread %i\n", thread );
                return;
        }
        if ( ( g_bspdata->dstaticprops[prop->propnum].flags & STATICPROPFLAGS_STATICLIGHTING ) == 0 )
        {
                // baked lighting not wanted
                return;
        }

        pvector<VDataDef> vdatas;
        
        // transform all vertices to be in world space
        prop->mdl.clear_model_nodes();
        prop->mdl.flatten_light();

        // we're not using NodePath::find_all_matches() because we need a consistent order in the list
        pvector<PT( GeomNode )> geomnodes = BuildGeomNodes( prop->mdl );
        for ( size_t i = 0; i < geomnodes.size(); i++ )
        {
                PT( GeomNode ) gn = geomnodes[i];
                for ( int j = 0; j < gn->get_num_geoms(); j++ )
                {
                        CPT( Geom ) geom = gn->get_geom( j );
                        CPT( GeomVertexData ) vdata = geom->get_vertex_data();

                        VDataDef def;
                        def.vdata = vdata;
                        def.mat_to_world = NodePath( gn ).get_net_transform()->get_mat();
                        vdatas.push_back( def );
                }
        }

        dstaticprop_t *dprop = &g_bspdata->dstaticprops[prop->propnum];

        //ThreadLock();
        dprop->first_vertex_data = g_bspdata->dstaticpropvertexdatas.size();

        for ( size_t i = 0; i < vdatas.size(); i++ )
        {
                CPT( GeomVertexData ) vdata = vdatas[i].vdata;
                GeomVertexReader vtx_reader( vdata, InternalName::get_vertex() );
                GeomVertexReader norm_reader( vdata, InternalName::get_normal() );

                dstaticpropvertexdata_t dvdata;
                dvdata.first_lighting_sample = g_bspdata->staticproplighting.size();

                for ( int row = 0; row < vdata->get_num_rows(); row++ )
                {
                        vtx_reader.set_row( row );
                        norm_reader.set_row( row );

                        LVector3 world_pos( vtx_reader.get_data3f() );
                        LNormalf world_normal( norm_reader.get_data3f() );

                        LVector3 direct_col( 0 );
                        ComputeDirectLightingAtPoint( world_pos, world_normal, direct_col );

                        LVector3 indirect_col( 0 );
                        ComputeIndirectLightingAtPoint( world_pos, world_normal, indirect_col, false );

                        colorrgbexp32_t sample;
                        VectorToColorRGBExp32( direct_col + indirect_col, sample );
                        g_bspdata->staticproplighting.push_back( sample );
                }
                dvdata.num_lighting_samples = g_bspdata->staticproplighting.size() - dvdata.first_lighting_sample;

                g_bspdata->dstaticpropvertexdatas.push_back( dvdata );
        }
        dprop->num_vertex_datas = g_bspdata->dstaticpropvertexdatas.size() - dprop->first_vertex_data;

        //ThreadUnlock();
}

void DoComputeStaticPropLighting()
{
        //NamedRunThreadsOnIndividual( (int)g_static_props.size(), g_estimate, ComputeStaticPropLighting );

        Log( "Computing static prop lighting...\n" );
        for ( size_t i = 0; i < g_static_props.size(); i++ )
        {
                ComputeStaticPropLighting( i );
        }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

INLINE RADCollisionPolygon::
RADCollisionPolygon( const LVecBase3 &a, const LVecBase3 &b,
                     const LVecBase3 &c )
{
        LPoint3 array[3];
        array[0] = a;
        array[1] = b;
        array[2] = c;
        setup_points( array, array + 3 );
}

/**
*
*/
INLINE RADCollisionPolygon::
RADCollisionPolygon( const LVecBase3 &a, const LVecBase3 &b,
                     const LVecBase3 &c, const LVecBase3 &d )
{
        LPoint3 array[4];
        array[0] = a;
        array[1] = b;
        array[2] = c;
        array[3] = d;
        setup_points( array, array + 4 );
}

/**
* Returns the number of vertices of the CollisionPolygon.
*/
INLINE int RADCollisionPolygon::
get_num_points() const
{
        return _points.size();
}

/**
* Returns the nth vertex of the CollisionPolygon, expressed in 3-D space.
*/
INLINE LPoint3 RADCollisionPolygon::
get_point( int n ) const
{
        nassertr( n >= 0 && n < (int)_points.size(), LPoint3::zero() );
        LMatrix4 to_3d_mat;
        rederive_to_3d_mat( to_3d_mat );
        return to_3d( _points[n]._p, to_3d_mat );
}

/**
* Returns true if the 2-d v1 is to the right of v2.
*/
INLINE bool RADCollisionPolygon::
is_right( const LVector2 &v1, const LVector2 &v2 )
{
        return ( v1[0] * v2[1] - v1[1] * v2[0] ) > 1.0e-6f;
}

/**
* Returns the linear distance of p to the line defined by f and f+v, where v
* is a normalized vector.  The result is negative if p is left of the line,
* positive if it is right of the line.
*/
INLINE PN_stdfloat RADCollisionPolygon::
dist_to_line( const LPoint2 &p,
              const LPoint2 &f, const LVector2 &v )
{
        LVector2 v1 = ( p - f );
        return ( v1[0] * v[1] - v1[1] * v[0] );
}

/**
* Assuming the indicated point in 3-d space lies within the polygon's plane,
* returns the corresponding point in the polygon's 2-d definition space.
*/
INLINE LPoint2 RADCollisionPolygon::
to_2d( const LVecBase3 &point3d ) const
{
        LPoint3 point = LPoint3( point3d ) * _to_2d_mat;
        return LPoint2( point[0], point[2] );
}

/**
* Fills the indicated matrix with the appropriate rotation transform to move
* points from the 2-d plane into the 3-d (X, 0, Z) plane.
*/
INLINE void RADCollisionPolygon::
calc_to_3d_mat( LMatrix4 &to_3d_mat ) const
{
        // We have to be explicit about the coordinate system--we specifically mean
        // CS_zup_right, because that points the forward vector down the Y axis and
        // moves the coords in (X, 0, Z).  We want this effect regardless of the
        // user's coordinate system of choice.

        // The up vector, on the other hand, is completely arbitrary.

        look_at( to_3d_mat, -get_plane().get_normal(),
                 LVector3( 0.0f, 0.0f, 1.0f ), CS_zup_right );
        to_3d_mat.set_row( 3, get_plane().get_point() );
}

/**
* Fills the indicated matrix with the appropriate rotation transform to move
* points from the 2-d plane into the 3-d (X, 0, Z) plane.
*
* This is essentially similar to calc_to_3d_mat, except that the matrix is
* rederived from whatever is stored in _to_2d_mat, guaranteeing that it will
* match whatever algorithm produced that one, even if it was produced on a
* different machine with different numerical precision.
*/
INLINE void RADCollisionPolygon::
rederive_to_3d_mat( LMatrix4 &to_3d_mat ) const
{
        to_3d_mat.invert_from( _to_2d_mat );
}

/**
* Extrude the indicated point in the polygon's 2-d definition space back into
* 3-d coordinates.
*/
INLINE LPoint3 RADCollisionPolygon::
to_3d( const LVecBase2 &point2d, const LMatrix4 &to_3d_mat )
{
        return LPoint3( point2d[0], 0.0f, point2d[1] ) * to_3d_mat;
}

/**
*
*/
INLINE RADCollisionPolygon::PointDef::
PointDef( const LPoint2 &p, const LVector2 &v ) : _p( p ), _v( v )
{
}

/**
*
*/
INLINE RADCollisionPolygon::PointDef::
PointDef( PN_stdfloat x, PN_stdfloat y ) : _p( x, y ), _v( 0.0f, 0.0f )
{
}

/**
*
*/
INLINE RADCollisionPolygon::PointDef::
PointDef( const RADCollisionPolygon::PointDef &copy ) : _p( copy._p ), _v( copy._v )
{
}

/**
*
*/
INLINE void RADCollisionPolygon::PointDef::
operator = ( const RADCollisionPolygon::PointDef &copy )
{
        _p = copy._p;
        _v = copy._v;
}

/**
* Returns the linear distance of p to the line segment defined by f and t,
* where v = (t - f).normalize(). The result is negative if p is left of the
* line, positive if it is right of the line.  If the result is positive, it
* is constrained by endpoints of the line segment (i.e.  the result might be
* larger than it would be for a straight distance-to-line test).  If the
* result is negative, we don't bother.
*/
PN_stdfloat RADCollisionPolygon::
dist_to_line_segment( const LPoint2 &p,
                      const LPoint2 &f, const LPoint2 &t,
                      const LVector2 &v )
{
        LVector2 v1 = ( p - f );
        PN_stdfloat d = ( v1[0] * v[1] - v1[1] * v[0] );
        if ( d < 0.0f )
        {
                return d;
        }

        // Compute the nearest point on the line.
        LPoint2 q = p + LVector2( -v[1], v[0] ) * d;

        // Now constrain that point to the line segment.
        if ( v[0] > 0.0f )
        {
                // X+
                if ( v[1] > 0.0f )
                {
                        // Y+
                        if ( v[0] > v[1] )
                        {
                                // X-dominant.
                                if ( q[0] < f[0] )
                                {
                                        return ( p - f ).length();
                                } if ( q[0] > t[0] )
                                {
                                        return ( p - t ).length();
                                }
                                else
                                {
                                        return d;
                                }
                        }
                        else
                        {
                                // Y-dominant.
                                if ( q[1] < f[1] )
                                {
                                        return ( p - f ).length();
                                } if ( q[1] > t[1] )
                                {
                                        return ( p - t ).length();
                                }
                                else
                                {
                                        return d;
                                }
                        }
                }
                else
                {
                        // Y-
                        if ( v[0] > -v[1] )
                        {
                                // X-dominant.
                                if ( q[0] < f[0] )
                                {
                                        return ( p - f ).length();
                                } if ( q[0] > t[0] )
                                {
                                        return ( p - t ).length();
                                }
                                else
                                {
                                        return d;
                                }
                        }
                        else
                        {
                                // Y-dominant.
                                if ( q[1] > f[1] )
                                {
                                        return ( p - f ).length();
                                } if ( q[1] < t[1] )
                                {
                                        return ( p - t ).length();
                                }
                                else
                                {
                                        return d;
                                }
                        }
                }
        }
        else
        {
                // X-
                if ( v[1] > 0.0f )
                {
                        // Y+
                        if ( -v[0] > v[1] )
                        {
                                // X-dominant.
                                if ( q[0] > f[0] )
                                {
                                        return ( p - f ).length();
                                } if ( q[0] < t[0] )
                                {
                                        return ( p - t ).length();
                                }
                                else
                                {
                                        return d;
                                }
                        }
                        else
                        {
                                // Y-dominant.
                                if ( q[1] < f[1] )
                                {
                                        return ( p - f ).length();
                                } if ( q[1] > t[1] )
                                {
                                        return ( p - t ).length();
                                }
                                else
                                {
                                        return d;
                                }
                        }
                }
                else
                {
                        // Y-
                        if ( -v[0] > -v[1] )
                        {
                                // X-dominant.
                                if ( q[0] > f[0] )
                                {
                                        return ( p - f ).length();
                                } if ( q[0] < t[0] )
                                {
                                        return ( p - t ).length();
                                }
                                else
                                {
                                        return d;
                                }
                        }
                        else
                        {
                                // Y-dominant.
                                if ( q[1] > f[1] )
                                {
                                        return ( p - f ).length();
                                } if ( q[1] < t[1] )
                                {
                                        return ( p - t ).length();
                                }
                                else
                                {
                                        return d;
                                }
                        }
                }
        }
}

/**
* Now that the _p members of the given points array have been computed, go
* back and compute all of the _v members.
*/
void RADCollisionPolygon::
compute_vectors( Points &points )
{
        size_t num_points = points.size();
        for ( size_t i = 0; i < num_points; i++ )
        {
                points[i]._v = points[( i + 1 ) % num_points]._p - points[i]._p;
                points[i]._v.normalize();
        }
}

/**
* Returns true if the indicated point is within the polygon's 2-d space,
* false otherwise.
*/
bool RADCollisionPolygon::
point_is_inside( const LPoint2 &p, const RADCollisionPolygon::Points &points ) const
{
        // We insist that the polygon be convex.  This makes things a bit simpler.

        // In the case of a convex polygon, defined with points in counterclockwise
        // order, a point is interior to the polygon iff the point is not right of
        // each of the edges.
        for ( int i = 0; i < (int)points.size() - 1; i++ )
        {
                if ( is_right( p - points[i]._p, points[i + 1]._p - points[i]._p ) )
                {
                        return false;
                }
        }
        if ( is_right( p - points[points.size() - 1]._p,
                       points[0]._p - points[points.size() - 1]._p ) )
        {
                return false;
        }

        return true;
}

/**
* Returns the linear distance from the 2-d point to the nearest part of the
* polygon defined by the points vector.  The result is negative if the point
* is within the polygon.
*/
PN_stdfloat RADCollisionPolygon::
dist_to_polygon( const LPoint2 &p, const RADCollisionPolygon::Points &points ) const
{

        // We know that that the polygon is convex and is defined with the points in
        // counterclockwise order.  Therefore, we simply compare the signed distance
        // to each line segment; we ignore any negative values, and take the minimum
        // of all the positive values.

        // If all values are negative, the point is within the polygon; we therefore
        // return an arbitrary negative result.

        bool got_dist = false;
        PN_stdfloat best_dist = -1.0f;

        size_t num_points = points.size();
        for ( size_t i = 0; i < num_points - 1; ++i )
        {
                PN_stdfloat d = dist_to_line_segment( p, points[i]._p, points[i + 1]._p,
                                                      points[i]._v );
                if ( d >= 0.0f )
                {
                        if ( !got_dist || d < best_dist )
                        {
                                best_dist = d;
                                got_dist = true;
                        }
                }
        }

        PN_stdfloat d = dist_to_line_segment( p, points[num_points - 1]._p, points[0]._p,
                                              points[num_points - 1]._v );
        if ( d >= 0.0f )
        {
                if ( !got_dist || d < best_dist )
                {
                        best_dist = d;
                        got_dist = true;
                }
        }

        return best_dist;
}

/**
* Projects the polygon onto the given axis, returning the center on the line
* and the half extent.
*/
void RADCollisionPolygon::
project( const LVector3 &axis, PN_stdfloat &center, PN_stdfloat &extent ) const
{
        PN_stdfloat begin, end;

        Points::const_iterator pi;
        pi = _points.begin();

        const LPoint2 &point = ( *pi )._p;
        begin = point[0] * axis[0] + point[1] * axis[2];
        end = begin;

        for ( ; pi != _points.end(); ++pi )
        {
                const LPoint2 &point = ( *pi )._p;

                PN_stdfloat t = point[0] * axis[0] + point[1] * axis[2];
                begin = std::min( begin, t );
                end = std::max( end, t );
        }

        center = ( end + begin ) * 0.5f;
        extent = cabs( ( end - begin ) * 0.5f );
}

/**
*
*/
void RADCollisionPolygon::
setup_points( const LPoint3 *begin, const LPoint3 *end )
{
        int num_points = end - begin;
        nassertv( num_points >= 3 );

        _points.clear();

        // Tell the base CollisionPlane class what its plane will be.  To do this,
        // we must first compute the polygon normal.
        LVector3 normal = LVector3::zero();

        // Project the polygon into each of the three major planes and calculate the
        // area of each 2-d projection.  This becomes the polygon normal.  This
        // works because the ratio between these different areas corresponds to the
        // angle at which the polygon is tilted toward each plane.
        for ( int i = 0; i < num_points; i++ )
        {
                const LPoint3 &p0 = begin[i];
                const LPoint3 &p1 = begin[( i + 1 ) % num_points];
                normal[0] += p0[1] * p1[2] - p0[2] * p1[1];
                normal[1] += p0[2] * p1[0] - p0[0] * p1[2];
                normal[2] += p0[0] * p1[1] - p0[1] * p1[0];
        }

        if ( normal.length_squared() == 0.0f )
        {
                // The polygon has no area.
                return;
        }

        set_plane( LPlane( normal, begin[0] ) );

        // Construct a matrix that rotates the points from the (X,0,Z) plane into
        // the 3-d plane.
        LMatrix4 to_3d_mat;
        calc_to_3d_mat( to_3d_mat );

        // And the inverse matrix rotates points from 3-d space into the 2-d plane.
        _to_2d_mat.invert_from( to_3d_mat );

        // Now project all of the points onto the 2-d plane.

        const LPoint3 *pi;
        for ( pi = begin; pi != end; ++pi )
        {
                LPoint3 point = ( *pi ) * _to_2d_mat;
                _points.push_back( PointDef( point[0], point[2] ) );
        }

        nassertv( _points.size() >= 3 );

        compute_vectors( _points );
}

/**
* Converts the indicated point to 3-d space according to the way
* CollisionPolygons used to be stored in bam files prior to 4.9.
*/
LPoint3 RADCollisionPolygon::
legacy_to_3d( const LVecBase2 &point2d, int axis ) const
{
        nassertr( !point2d.is_nan(), LPoint3( 0.0f, 0.0f, 0.0f ) );

        LVector3 normal = get_normal();
        PN_stdfloat D = get_plane()[3];

        nassertr( !normal.is_nan(), LPoint3( 0.0f, 0.0f, 0.0f ) );
        nassertr( !cnan( D ), LPoint3( 0.0f, 0.0f, 0.0f ) );

        switch ( axis )
        {
        case 0:  // AT_x:
                return LPoint3( -( normal[1] * point2d[0] + normal[2] * point2d[1] + D ) / normal[0], point2d[0], point2d[1] );

        case 1:  // AT_y:
                return LPoint3( point2d[0],
                                -( normal[0] * point2d[0] + normal[2] * point2d[1] + D ) / normal[1], point2d[1] );

        case 2:  // AT_z:
                return LPoint3( point2d[0], point2d[1],
                                -( normal[0] * point2d[0] + normal[1] * point2d[1] + D ) / normal[2] );
        }

        nassertr( false, LPoint3( 0.0f, 0.0f, 0.0f ) );
        return LPoint3( 0.0f, 0.0f, 0.0f );
}

bool RADCollisionPolygon::segment_intersection_test( const vec3_t start, const vec3_t stop )
{
        if ( _points.size() < 3 )
        {
                return false;
        }

        LPoint3 from_a( start[0], start[1], start[2] );
        LPoint3 from_b( stop[0], stop[1], stop[2] );
        LPoint3 from_direction = from_b - from_a;

        PN_stdfloat t;
        if ( !get_plane().intersects_line( t, from_a, from_direction ) )
        {
                // No intersection.
                return false;
        }

        if ( t < 0.0f || t > 1.0f )
        {
                // The intersection point is before the start of the segment or after the
                // end of the segment.
                return false;
        }

        LPoint3 plane_point = from_a + t * from_direction;
        LPoint2 p = to_2d( plane_point );

        // No clip plane is in effect.  Do the default test.
        if ( !point_is_inside( p, _points ) )
        {
                return false;
        }

        // The segment intersects our polygon
        return true;
}