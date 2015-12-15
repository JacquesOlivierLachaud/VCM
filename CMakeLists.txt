# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(VCM)

cmake_minimum_required (VERSION 2.6) 


FIND_PACKAGE(DGtal 0.9 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

 # -------------------------------------------------------------------------
 # This test is for instance used for ITK v3.x. ITK forces a limited
 # template depth. We remove this option.
 # --------------------------------------------------------------------------
 if (CMAKE_CXX_FLAGS MATCHES "-ftemplate-depth-[0-9]*")
    message( "Warning: some package has enabled a limited template depth for the C++ compiler." )
    message( "         Disabling option -ftemplate-depth-xx in CMAKE_CXX_FLAGS." )
    set( CMAKE_CXX_FLAGS_TMP ${CMAKE_CXX_FLAGS} )
    STRING( REGEX REPLACE "-ftemplate-depth-[0-9]*" ""
    CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_TMP}" )
     message ("         CMAKE_CXX_FLAGS=" ${CMAKE_CXX_FLAGS} )
  endif (CMAKE_CXX_FLAGS MATCHES "-ftemplate-depth-[0-9]*")




# -----------------------------------------------------------------------------
# Define variables and cmake parameters, and display some information
# -----------------------------------------------------------------------------
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
SET(VCM_VERSION_MAJOR 0)
SET(VCM_VERSION_MINOR 6)
SET(VCM_VERSION_PATCH 0)
SET(DGTALTOOLS_VERSION "${VCM_VERSION_MAJOR}.${VCM_VERSION_MINOR}.${VCM_VERSION_PATCH}")
SET(PROJECT_VERSION "${VCM_VERSION_MAJOR}.${VCM_VERSION_MINOR}.${VCM_VERSION_PATCH}")
SET(VERSION ${VCM_VERSION_MAJOR}.${VCM_VERSION_MINOR}.${VCM_VERSION_PATCH})

SET(Boost_USE_STATIC_LIBS   OFF)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
  SET(VCMLibDependencies ${VCMLibDependencies} 
     ${Boost_LIBRAIRIES}  
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )

#------Tools subdirectories ------

# Add source files not needing QGLViewer here.
SET(VCM_SRC
  dvcm-2d
  delta-distance
  delta-distance3D
  delta-vcm
  )
# Add source files needing QGLViewer here.
SET(VCM_QGL_SRC
  dvcm-3d
  dvcm-features-3d
  implicitShape3NormalEstimation
  noisyImplicitShape3NormalEstimation
  generic3dNormalEstimators
  angle-viewer
  area-integration
  vol2offAndNormals
  heightfield2offAndNormals
  )

#add_definitions( -DWITH_EIGEN=YES )

FOREACH(FILE ${VCM_SRC})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${VCMLibDependencies})
  SET_TARGET_PROPERTIES(${FILE} PROPERTIES LINKER_LANGUAGE CXX)
  install(TARGETS ${FILE} 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)


if ( WITH_VISU3D_QGLVIEWER )
  FOREACH(FILE ${VCM_QGL_SRC}) 
    add_executable(${FILE} ${FILE})
    target_link_libraries ( ${FILE} ${DGTAL_LIBRARIES} ${VCMLibDependencies})   
	SET_TARGET_PROPERTIES(${FILE} PROPERTIES LINKER_LANGUAGE CXX)
    install(TARGETS ${FILE} 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
  ENDFOREACH(FILE)
endif ( WITH_VISU3D_QGLVIEWER )


