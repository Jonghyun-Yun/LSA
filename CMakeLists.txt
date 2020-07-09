cmake_minimum_required(VERSION 3.0)
 project(lsjm-code)

 set(MATH /Users/yunj/stan-dev/math)
 set(PROJROOT /Users/yunj/Dropbox/research/lsjm-art/lsjm-code)
 set(Xcode_CLT /Library/Developer/CommandLineTools)

 set(CMAKE_CXX_FLAGS_DEBUG "-g")
 set(CMAKE_CXX_FLAGS_RELEASE "-O3")

 # set(CMAKE_CXX_STANDARD 11)

 # # STAN math lib
 include_directories(${MATH}/lib/tbb_2019_U8/include)
 include_directories(${MATH})
 include_directories(${MATH}/lib/eigen_3.3.7)
 include_directories(${MATH}/lib/boost_1.72.0)
 include_directories(${MATH}/lib/sundials_5.2.0/include)

 # xcode commandline toolchain
 include_directories(${Xcode_CLT}/usr/include/c++/v1)
 include_directories(${Xcode_CLT}/usr/lib/clang/11.0.0/include)
 include_directories(${Xcode_CLT}/usr/include)
 include_directories(${Xcode_CLT}/SDKs/MacOSX.sdk/usr/include)

 # set compiler
 # SET(CMAKE_CXX_COMPILER "/usr/bin/clang++")
 # SET(CMAKE_CC_COMPILER "/usr/bin/clnag")

 # c++ linker
 # SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -isystem/Library/Developer/CommandLineTools/usr/include/c++/v1")
 SET(CMAKE_CXX_FLAGS  "-std=c++1y -Wno-unknown-warning-option -Wno-tautological-compare -Wno-sign-compare -D_REENTRANT ${CMAKE_CXX_FLAGS}")
 SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DSTAN_THREADS")
 SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DBOOST_DISABLE_ASSERTS")
 SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wl,-L,${MATH}/lib/tbb -Wl,-rpath,${MATH}/lib/tbb")
 SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -isysroot ${Xcode_CLT}/SDKs/MacOSX.sdk/usr/include")
 SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk")
 # SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-unused-local-typedef -Wno-unused-function -Wno-ignored-qualifiers -Wno-unused-variable -Wno-sign-compare -Wno-infinite-recursion")

include_directories(${PROJROOT}/include)
file(GLOB SOURCES
    ${PROJROOT}/include/*.h
    ${PROJROOT}/include/*.hpp
    ${PROJROOT}/src/*.cpp
)

 add_executable(main main.cpp ${SOURCES})
 target_link_libraries(main ${CMAKE_CXX_FLAGS} ${MATH}/lib/sundials_5.2.0/lib/libsundials_nvecserial.a ${MATH}/lib/sundials_5.2.0/lib/libsundials_cvodes.a ${MATH}/lib/sundials_5.2.0/lib/libsundials_idas.a ${MATH}/lib/sundials_5.2.0/lib/libsundials_kinsol.a ${MATH}/lib/tbb/libtbb.dylib ${MATH}/lib/tbb/libtbbmalloc.dylib ${MATH}/lib/tbb/libtbbmalloc_proxy.dylib)
