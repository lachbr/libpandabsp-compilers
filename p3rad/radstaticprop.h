#ifndef RADSTATICPROP_H
#define RADSTATICPROP_H

#include <nodePath.h>
#include <bitMask.h>
#include <pvector.h>
#include <collisionSegment.h>
#include <collisionPlane.h>
#include <lightMutex.h>

#include "mathtypes.h"

// A more efficient, general purpose collision polygon which can do direct testing with a CollisionSegment
// without the overhead of using a CollisionHandlerQueue and CollisionTraverser, which serves no purpose
// for RAD.
class RADCollisionPolygon : public CollisionPlane
{
public:
	INLINE RADCollisionPolygon( const LVecBase3 &a, const LVecBase3 &b,
				 const LVecBase3 &c );
	INLINE RADCollisionPolygon( const LVecBase3 &a, const LVecBase3 &b,
				 const LVecBase3 &c, const LVecBase3 &d );

	bool segment_intersection_test( const vec3_t start, const vec3_t end );

	INLINE int get_num_points() const;
	INLINE LPoint3 get_point( int n ) const;

	INLINE static bool is_right( const LVector2 &v1, const LVector2 &v2 );
	INLINE static PN_stdfloat dist_to_line( const LPoint2 &p,
						const LPoint2 &f, const LVector2 &v );
	static PN_stdfloat dist_to_line_segment( const LPoint2 &p,
						 const LPoint2 &f, const LPoint2 &t,
						 const LVector2 &v );

	class PointDef
	{
	public:
		INLINE PointDef( const LPoint2 &p, const LVector2 &v );
		INLINE PointDef( PN_stdfloat x, PN_stdfloat y );
		INLINE PointDef( const PointDef &copy );
		INLINE void operator = ( const PointDef &copy );

		LPoint2 _p;  // the point in 2-d space
		LVector2 _v; // the normalized vector to the next point
	};
	typedef pvector<PointDef> Points;

	static void compute_vectors( Points &points );

	bool point_is_inside( const LPoint2 &p, const Points &points ) const;
	PN_stdfloat dist_to_polygon( const LPoint2 &p, const Points &points ) const;
	void project( const LVector3 &axis, PN_stdfloat &center, PN_stdfloat &extent ) const;

	void setup_points( const LPoint3 *begin, const LPoint3 *end );
	INLINE LPoint2 to_2d( const LVecBase3 &point3d ) const;
	INLINE void calc_to_3d_mat( LMatrix4 &to_3d_mat ) const;
	INLINE void rederive_to_3d_mat( LMatrix4 &to_3d_mat ) const;
	INLINE static LPoint3 to_3d( const LVecBase2 &point2d, const LMatrix4 &to_3d_mat );
	LPoint3 legacy_to_3d( const LVecBase2 &point2d, int axis ) const;

	Points _points;
	LMatrix4 _to_2d_mat;

public:
	static TypeHandle get_class_type()
	{
		return _type_handle;
	}
	static void init_type()
	{
		CollisionPlane::init_type();
		register_type( _type_handle, "RADCollisionPolygon",
			       CollisionPlane::get_class_type() );
	}
	virtual TypeHandle get_type() const
	{
		return get_class_type();
	}
	virtual TypeHandle force_init_type()
	{
		init_type(); return get_class_type();
	}

private:
	static TypeHandle _type_handle;
};

extern bool g_collisions_loaded;

extern LightMutex g_prop_lock;

struct RADStaticProp
{
	int leafnum;
	pvector<PT( RADCollisionPolygon )> polygons;
};

struct TestGroup
{
	PT( CollisionSegment ) seg;
};
// One for each thread
extern pvector<TestGroup> g_test_groups;

extern pvector<RADStaticProp *> g_static_props;

extern void LoadStaticProps();
extern bool StaticPropIntersectionTest( const vec3_t start, const vec3_t stop, int leafnum );

#endif // RADSTATICPROP_H