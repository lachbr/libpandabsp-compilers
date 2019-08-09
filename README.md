# libpandabsp-compilers
`libpandabsp-compilers` is a suite of tools for compiling BSP levels, for use with my [libpandabsp](https://github.com/lachbr/libpandabsp) library and the [Panda3D](https://github.com/panda3d/panda3d) game engine. It is a heavy modification and upgrade of [Zoner's Half Life Tools](http://zhlt.info/).

---

### How does it work?
The compilation of a complete Panda3D **.bsp** file is a four step process.
1) CSG (Constructive Solid Geometry)
2) BSP (Binary Space Partitioning)
3) VIS (Visibility)
4) RAD (Lighting/Radiosity)

---

### CSG (Constructive Solid Geometry)
This is the first and simplest step in the compilation process of a level.

The purpose of this step is to perform boolean operations on the brushes from your map file and output a list of convex polygons, which represent the final geometry of your level.

The two outputs from `p3csg.exe` are:
- A .bsp file containing only a list of planes, brushes, and brush sides
- A raw text file containing the list of polygons (from the CSG operations on the brushes) along with their cooresponding plane, brush, and brush side indices in the .bsp file.
### BSP (Binary Space Partitioning)
Using the .bsp file and list of polygons from the CSG step, the BSP step builds a binary space partition tree and outputs it to the .bsp file.
The polygons are used to determine the leaf nodes of the tree. If a polygon spans over multiple leaves, it will be split up until each polygon resides in a single leaf.
### VIS (Visibility)
Using the BSP tree emitted to the .bsp file, this step figures out for each leaf node, what other leaf nodes are potentially visible? A compressed array of leaf index->list of potentially visible leaves is emitted to the .bsp file.
### RAD (Lighting/Radiosity)
This is the most computationally expensive and lengthiest step. It will create a lightmap for each face in the .bsp file using the light entities put into the map. Along with direct lighting, it will also perform radiosity passes to simulate light bouncing.
This step will also bake per-vertex lighting onto any static props placed in the map, and bake environment probes all throughout the map (for ambient lighting on dynamic objects).

Once all four of these steps complete, the .bsp file is ready to be loaded into your Panda3D game.
