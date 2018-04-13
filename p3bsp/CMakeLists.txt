set(P3BSP_CPPFILES ${COMMON_CPPFILES}
                   brink.cpp
                   merge.cpp
                   outside.cpp
                   portals.cpp
                   qbsp.cpp
                   solidbsp.cpp
                   surfaces.cpp
                   tjunc.cpp
                   writebsp.cpp)

set(P3BSP_INCLUDEFILES ${COMMON_INCLUDEFILES}
                       bsp5.h)
                       
add_definitions(-DHLBSP -DDOUBLEVEC_T)
add_executable(p3bspcm ${P3BSP_CPPFILES} ${P3BSP_INCLUDEFILES})