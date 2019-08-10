project(bsp_common)

file (GLOB SRCS "*.cpp")
file (GLOB HEADERS "*.h")
file (GLOB ML_SRCS "mathlib/*.cpp")
file (GLOB ML_HDRS "mathlib/*.h")
file (GLOB CL_HDRS "CL/*.h")

source_group("Header Files" FILES ${HEADERS})
source_group("Source Files" FILES ${SRCS})
source_group("Header Files\\mathlib" FILES ${ML_HDRS})
source_group("Source Files\\mathlib" FILES ${ML_SRCS})
source_group("Header Files\\CL" FILES ${CL_HDRS})

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type.")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /DEBUG")
string(REPLACE "/DNDEBUG" "" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
string(REPLACE "/DNDEBUG" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

if (DEFINED LINK_ALL_STATIC)
    set(P3BUILT built_x64_static)
else()
    set(P3BUILT built_x64)
endif ()

set(PANDA3D_DIR ..\\..\\..\\cio\\cio-panda3d\\${P3BUILT} CACHE PATH "Panda3D directory")

add_definitions(-DBUILDING_BSPCOMMON -DWITHIN_PANDA -DNOMINMAX -DSTDC_HEADERS)
if (DEFINED LINK_ALL_STATIC)
    set(LIB_TYPE STATIC)
    add_definitions(-DLINK_ALL_STATIC)
else ()
    set(LIB_TYPE SHARED)
endif ()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    # 64 bits
    add_definitions(-DVERSION_64BIT)
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # 32 bits
    add_definitions(-DVERSION_32BIT)
endif()

include_directories(./ ${PANDA3D_DIR}/include)
link_directories(${PANDA3D_DIR}/lib)

add_library(bsp_common ${LIB_TYPE} ${SRCS} ${HEADERS} ${ML_SRCS} ${ML_HDRS} ${CL_HDRS})

STRING(REGEX REPLACE "/" "\\\\" WIN_PANDA3D_DIR ${PANDA3D_DIR})

if (DEFINED LINK_ALL_STATIC)
    add_custom_command(TARGET bsp_common POST_BUILD COMMAND copy /Y ${CMAKE_BUILD_TYPE}\\bsp_common.lib ${WIN_PANDA3D_DIR}\\lib\\bsp_common.lib)
else()
    add_custom_command(TARGET bsp_common POST_BUILD COMMAND copy /Y ${CMAKE_BUILD_TYPE}\\bsp_common.dll ${WIN_PANDA3D_DIR}\\bin\\bsp_common.dll)
    add_custom_command(TARGET bsp_common POST_BUILD COMMAND copy /Y ${CMAKE_BUILD_TYPE}\\bsp_common.pdb ${WIN_PANDA3D_DIR}\\bin\\bsp_common.pdb)
endif()

target_link_libraries(bsp_common
					  libpanda.lib
					  libpandaexpress.lib
					  libp3dtool.lib
					  libp3dtoolconfig.lib)
