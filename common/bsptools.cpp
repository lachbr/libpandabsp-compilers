#include "bsptools.h"

BaseBSPEnumerator::BaseBSPEnumerator( bspdata_t *dat ) :
        data( dat )
{
}

bool r_enumerate_nodes_along_ray( int node_id, const Ray &ray, float start,
                                  float end, BaseBSPEnumerator *surf, int context,
                                  float scale )
{
        float front, back;
        float start_dot_n, delta_dot_n;

        while ( node_id >= 0 )
        {
                dnode_t *node = &surf->data->dnodes[node_id];
                dplane_t *plane = &surf->data->dplanes[node->planenum];

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

                front = start_dot_n + start * delta_dot_n - (plane->dist / scale);
                back = start_dot_n + end * delta_dot_n - (plane->dist / scale);

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
                                split_frac = ( (plane->dist / scale) - start_dot_n ) / delta_dot_n;
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
                                                              split_frac, surf, context, scale );
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
                                                            end, surf, context, scale );
                }
        }

        // Visit the leaf...
        return surf->enumerate_leaf( ~node_id, ray, start, end, context );
}

bool enumerate_nodes_along_ray( const Ray &ray, BaseBSPEnumerator *surf, int context, float scale )
{
        return r_enumerate_nodes_along_ray( 0, ray, 0.0, 1.0, surf, context, scale );
}