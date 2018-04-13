cmake_minimum_required (VERSION 2.6)

project(panda_bsp_toolscm)

set(COMMON_CPPFILES "../common/blockmem.cpp"
                    "../common/bspfile.cpp"
                    "../common/cmdlib.cpp"
                    "../common/cmdlinecfg.cpp"
                    "../common/filelib.cpp"
                    "../common/log.cpp"
                    "../common/mathlib.cpp"
                    "../common/messages.cpp"
                    "../common/scriplib.cpp"
                    "../common/threads.cpp"
                    "../common/winding.cpp" [INTERNAL])
                    
set(COMMON_INCLUDEDIRS "../template ../common" [INTERNAL])

set(COMMON_INCLUDEFILES "../template/basictypes.h"
                        "../common/blockmem.h"
                        "../common/boundingbox.h"
                        "../common/bspfile.h"
                        "../common/cmdlib.h"
                        "../common/cmdlinecfg.h"
                        "../common/filelib.h"
                        "../common/hlassert.h"
                        "../common/log.h"
                        "../common/mathlib.h"
                        "../common/mathtypes.h"
                        "../common/messages.h"
                        "../common/scriplib.h"
                        "../common/threads.h"
                        "../common/win32fix.h"
                        "../common/winding.h" [INTERNAL])
                        
set(COMMON_DEFINITIONS "VERSION_LINUX"
                       "SYSTEM_POSIX"
                       "NDEBUG"
                       "STDC_HEADERS"
                       "HAVE_FCNTL_H"
                       "HAVE_PTHREAD_H"
                       "HAVE_SYS_RESOURCE_H"
                       "HAVE_SYS_STAT_H"
                       "HAVE_SYS_TIME_H"
                       "HAVE_UNISTD_H" [INTERNAL])
                       
set(COMMON_FLAGS "-Wall -O2 -fno-strict-aliasing -pthread -pipe" INTERNAL)

add_subdirectory(p3csg)
add_subdirectory(p3bsp)
add_subdirectory(p3vis)
add_subdirectory(p3rad)
add_subdirectory(ripent)
                    