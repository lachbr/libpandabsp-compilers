#include "qrad.h"
#include <virtualFileSystem.h>
#include <texturePool.h>

#include <viftokenizer.h>
#include <vifparser.h>

#ifdef WORDS_BIGENDIAN
#error "HLRAD_TEXTURE doesn't support WORDS_BIGENDIAN, because I have no big endian machine to test it"
#endif

int g_numtextures;
radtexture_t *g_textures;

void DefaultTexture( radtexture_t *tex, const char *name )
{
        int i;
        PNMImage *img = new PNMImage( DEFAULT_LIGHTMAP_SIZE, DEFAULT_LIGHTMAP_SIZE );
        img->fill( 1.0 );
        tex->image = img;
        strcpy( tex->name, name );
        tex->name[MAX_TEXTURE_NAME - 1] = '\0';
        tex->bump = nullptr;
}

void LoadTextures()
{
        if ( !g_notextures )
        {
                Log( "Load Textures:\n" );
        }

        VirtualFileSystem *vfs = VirtualFileSystem::get_global_ptr();

        g_numtextures = g_bspdata->numtexrefs;
        g_textures = (radtexture_t *)malloc( g_numtextures * sizeof( radtexture_t ) );
        hlassume( g_textures != NULL, assume_NoMemory );
        int i;
        for ( i = 0; i < g_numtextures; i++ )
        {
                radtexture_t *tex = &g_textures[i];
                tex->bump = nullptr;
                tex->image = nullptr;

                if ( g_notextures )
                {
                        DefaultTexture( tex, "DEFAULT" );
                }
                else
                {
                        texref_t *tref = &g_bspdata->dtexrefs[i];
                        string name = tref->name;
                        Filename fname = Filename::from_os_specific( name );
                        if ( !vfs->exists( fname ) )
                        {
                                Warning( "Material %s not found", fname.get_fullpath().c_str() );
                                continue;
                        }
                                
                        string mdata = vfs->read_file( fname, true );
                        TokenVec toks = tokenizer( mdata );
                        Parser p( toks );
                        Object obj = p._base_objects[0];
                        std::cout << obj.name << std::endl;
                        for ( size_t n = 0; n < obj.properties.size(); n++ )
                        {
                                std::cout << "\t" << obj.properties[n].name << "\t:\t" << obj.properties[n].value << std::endl;
                        }

                        if ( !p.has_property( obj, "$basetexture" ) )
                        {
                                Warning( "Material %s has no basetexture", fname.get_fullpath().c_str() );
                                continue;
                        }
                                

                        size_t ext_idx = name.find_last_of( "." );
                        string basename = name.substr( 0, ext_idx );
                        string ext = name.substr( ext_idx );

                        PNMImage *img = new PNMImage;
                        if ( img->read( Filename::from_os_specific(
                                p.get_property_value( obj, "$basetexture" ) ) ) )
                        {
                                tex->image = img;
                                strcpy( tex->name, name.c_str() );
                                tex->name[MAX_TEXTURE_NAME - 1] = '\0';

                                PNMImage *bump = new PNMImage;
                                if ( p.has_property( obj, "$bumpmap" ) )
                                {
                                        Log( "Loaded RAD texture from %s and corresponding bump map %s.\n", tref->name,
                                             p.get_property_value( obj, "$bumpmap" ).c_str() );
                                        tex->bump = bump;
                                }
                                else
                                {
                                        Log( "Loaded RAD texture from %s.\n", tref->name );
                                        delete bump;
                                }
                                
                        }
                        else
                        {
                                Warning( "Could not load texture %s!\n", tref->name );
                                delete img;
                        }
                }


                {
                        vec3_t total;
                        VectorClear( total );
                        int width = tex->image->get_x_size();
                        int height = tex->image->get_y_size();
                        PNMImage *img = tex->image;
                        for ( int row = 0; row < height; row++ )
                        {
                                for ( int col = 0; col < width; col++ )
                                {
                                        vec3_t reflectivity;
                                        VectorCopy( img->get_xel( col, row ), reflectivity );
                                        for ( int k = 0; k < 3; k++ )
                                        {
                                                reflectivity[k] = pow( reflectivity[k], g_texreflectgamma );
                                        }
                                        VectorScale( reflectivity, g_texreflectscale, reflectivity );
                                        VectorAdd( total, reflectivity, total );
                                }
                        }
                        VectorScale( total, 1.0 / (double)( width * height ), total );
                        VectorCopy( total, tex->reflectivity );
                        Developer( DEVELOPER_LEVEL_MESSAGE, "Texture '%s': reflectivity is (%f,%f,%f).\n",
                                   tex->name, tex->reflectivity[0], tex->reflectivity[1], tex->reflectivity[2] );
                        if ( VectorMaximum( tex->reflectivity ) > 1.0 + NORMAL_EPSILON )
                        {
                                Warning( "Texture '%s': reflectivity (%f,%f,%f) greater than 1.0.", tex->name, tex->reflectivity[0], tex->reflectivity[1], tex->reflectivity[2] );
                        }
                }
        }
        if ( !g_notextures )
        {
                Log( "%i textures referenced\n", g_numtextures );
        }
}

