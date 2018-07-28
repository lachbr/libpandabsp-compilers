#include "qrad.h"
#include <virtualFileSystem.h>
#include <texturePool.h>

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
}

/*

void LoadTexture (radtexture_t *tex, const texref_t *mt, int size)
{
int i, j;
const texref_t *header = mt;
const byte *data = (const byte *)mt;
//tex->width = header->width;
//tex->height = header->height;
tex->width = DEFAULT_LIGHTMAP_SIZE;
tex->height = DEFAULT_LIGHTMAP_SIZE;
strcpy (tex->name, header->name);
tex->name[MAX_TEXTURE_NAME - 1] = '\0';
if (tex->width <= 0 || tex->height <= 0) //||
//	tex->width % (2 * 1 << (MIPLEVELS - 1)) != 0 || tex->height % (2 * (1 << (MIPLEVELS - 1))) != 0)
{
//Error ("Texture '%s': dimension (%dx%d) is not multiple of %d.", tex->name, tex->width, tex->height, 2 * (1 << (MIPLEVELS - 1)));
Error("Texture '%s': invalid dimensions (%dx%d)", tex->name, tex->width, tex->height);
}

int mipsize;
for (mipsize = 0, i = 0; i < MIPLEVELS; i++)
{
if ((int)mt->offsets[i] != (int)sizeof (texref_t) + mipsize)
{
Error ("Texture '%s': unexpected miptex offset.", tex->name);
}
mipsize += (tex->width >> i) * (tex->height >> i);
}
if (size < (int)sizeof (texref_t) + mipsize + 2 + 256 * 3)
{
Error ("Texture '%s': no enough data.", tex->name);
}
if (*(unsigned short *)&data[sizeof (texref_t) + mipsize] != 256)
{
Error ("Texture '%s': palette size is not 256.", tex->name);
}

tex->canvas = (byte *)malloc (tex->width * tex->height);
hlassume (tex->canvas != NULL, assume_NoMemory);
for (i = 0; i < tex->height; i++)
{
for (j = 0; j < tex->width; j++)
{
tex->canvas[i * tex->width + j] = data[sizeof (texref_t) + i * tex->width + j];
}
}
for (i = 0; i < 256; i++)
{
for (j = 0; j < 3; j++)
{
tex->palette[i][j] = data[sizeof (texref_t) + mipsize + 2 + i * 3 + j];
}
}
}

void LoadTextureFromWad (radtexture_t *tex, const texref_t *header)
{
tex->width = DEFAULT_LIGHTMAP_SIZE;
tex->height = DEFAULT_LIGHTMAP_SIZE;
strcpy (tex->name, header->name);
tex->name[MAX_TEXTURE_NAME - 1] = '\0';
wadfile_t *wad;
for (wad = g_wadfiles; wad; wad = wad->next)
{
lumpinfo_t temp, *found;
strcpy (temp.name, tex->name);
found = (lumpinfo_t *)bsearch (&temp, wad->lumpinfos, wad->numlumps, sizeof (lumpinfo_t), lump_sorter_by_name);
if (found)
{
Developer (DEVELOPER_LEVEL_MESSAGE, "Texture '%s': found in '%s'.\n", tex->name, wad->path);
if (found->type != 67 || found->compression != 0)
continue;
if (found->disksize < (int)sizeof (texref_t) || found->filepos < 0 || found->filepos + found->disksize > wad->filesize)
{
Warning ("Texture '%s': invalid texture data in '%s'.", tex->name, wad->path);
continue;
}
texref_t *mt = (texref_t *)malloc (found->disksize);
hlassume (mt != NULL, assume_NoMemory);
if (fseek (wad->file, found->filepos, SEEK_SET))
Error ("File read failure");
SafeRead (wad->file, mt, found->disksize);
if (!TerminatedString(mt->name, 16))
{
Warning("Texture '%s': invalid texture data in '%s'.", tex->name, wad->path);
free (mt);
continue;
}
Developer (DEVELOPER_LEVEL_MESSAGE, "Texture '%s': name '%s', width %d, height %d.\n", tex->name, mt->name, mt->width, mt->height);
if (strcasecmp (mt->name, tex->name))
{
Warning("Texture '%s': texture name '%s' differs from its reference name '%s' in '%s'.", tex->name, mt->name, tex->name, wad->path);
}
LoadTexture (tex, mt, found->disksize);
free (mt);
break;
}
}
if (!wad)
{
Warning ("Texture '%s': texture is not found in wad files.", tex->name);
DefaultTexture (tex, tex->name);
return;
}
}
*/
void LoadTextures()
{
        if ( !g_notextures )
        {
                Log( "Load Textures:\n" );
        }

        VirtualFileSystem *vfs = VirtualFileSystem::get_global_ptr();
        for ( size_t i = 0; i < g_multifiles.size(); i++ )
        {
                Filename file = Filename::from_os_specific( g_multifiles[i] );
                cout << file.get_fullpath() << endl;
                if ( !vfs->mount( file, ".", VirtualFileSystem::MF_read_only ) )
                {
                        Warning( "Could not mount multifile from %s!\n", file.get_fullpath().c_str() );
                }
                else
                {
                        Log( "Mounted multifile %s.\n", file.get_fullpath().c_str() );
                }
        }

        g_numtextures = g_numtexrefs;
        g_textures = (radtexture_t *)malloc( g_numtextures * sizeof( radtexture_t ) );
        hlassume( g_textures != NULL, assume_NoMemory );
        int i;
        for ( i = 0; i < g_numtextures; i++ )
        {
                radtexture_t *tex = &g_textures[i];

                if ( g_notextures )
                {
                        DefaultTexture( tex, "DEFAULT" );
                }
                else
                {
                        texref_t *tref = &g_dtexrefs[i];
                        string name = tref->name;

                        PNMImage *img = new PNMImage;
                        if ( img->read( name ) )
                        {
                                tex->image = img;
                                strcpy( tex->name, name.c_str() );
                                tex->name[MAX_TEXTURE_NAME - 1] = '\0';
                                Log( "Loaded RAD texture from %s.\n", tref->name );
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
                                        if ( img->get_num_channels() == 1 && img->get_channel( col, row, 0 ) * 0xFF == 0xFF )
                                        {
                                                VectorFill( reflectivity, 0.0 );
                                        }
                                        else
                                        {
                                                VectorScale( img->get_xel( col, row ) * 0xFF, 1.0 / 255.0, reflectivity );
                                                for ( int k = 0; k < 3; k++ )
                                                {
                                                        reflectivity[k] = pow( reflectivity[k], g_texreflectgamma );
                                                }
                                                VectorScale( reflectivity, g_texreflectscale, reflectivity );
                                        }
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

