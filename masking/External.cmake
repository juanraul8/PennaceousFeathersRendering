##################################################################################
# LIBPNG
##################################################################################
find_package(PNG)
if(PNG_FOUND)
	include_directories(${PNG_INCLUDE_DIR})
	list(APPEND cimg_libs ${PNG_LIBRARY})
	list(APPEND cimg_defs "cimg_use_png")
	list(APPEND svg_cpp_plot_libs ${PNG_LIBRARY})
	list(APPEND svg_cpp_plot_defs "USE_PNG")
endif(PNG_FOUND)

find_package(JPEG)
if(JPEG_FOUND)
	include_directories(${JPEG_INCLUDE_DIR})
	list(APPEND cimg_libs ${JPEG_LIBRARY})
	list(APPEND cimg_defs "cimg_use_jpeg")
endif(JPEG_FOUND)

find_package(TIFF)
if(TIFF_FOUND)
	include_directories(${TIFF_INCLUDE_DIR})
	list(APPEND cimg_libs ${TIFF_LIBRARY})
	list(APPEND image_defs "cimg_use_tiff")
endif(TIFF_FOUND)

if(FFTW3_FOUND)
	include_directories(${FFTW3_INCLUDE_DIR})
	list(APPEND cimg_libs ${FFTW3_LIBRARY})
	list(APPEND cimg_defs "cimg_use_fftw3")
endif(FFTW3_FOUND)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Windows specific code
    list(APPEND cimgdisplay_libs gdi32)
    list(APPEND cimgdisplay_defs cimg_display=2)
    list(APPEND cimg_defs cimg_OS=2)
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    list(APPEND cimg_libs pthread)
    list(APPEND cimgdisplay_libs X11)
    list(APPEND cimgdisplay_defs cimg_display=1)
    list(APPEND cimg_defs cimg_OS=1)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")


######################################################################
# EXTERNAL LIBRARIES
######################################################################
if (NOT EXTERNAL_INSTALL_LOCATION)
	set(EXTERNAL_INSTALL_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/ext)
endif()
if (NOT IS_DIRECTORY ${EXTERNAL_INSTALL_LOCATION})
	file(MAKE_DIRECTORY ${EXTERNAL_INSTALL_LOCATION})
endif()

include(ExternalProject)
# External include directory
add_custom_target(update)

ExternalProject_Add(svg.cc
  GIT_REPOSITORY https://github.com/adolfomunoz/svg.cc.git
  SOURCE_DIR ${EXTERNAL_INSTALL_LOCATION}/svg.cc
  GIT_TAG origin/main
  UPDATE_DISCONNECTED 1
  STEP_TARGETS update
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)
add_dependencies(update svg.cc-update)

include_directories(${EXTERNAL_INSTALL_LOCATION}/svg.cc/ext)
include_directories(${EXTERNAL_INSTALL_LOCATION}/svg.cc/ext/patterns/ext)
#This below does not seem to work so we add the above
file(GLOB dirs ${EXTERNAL_INSTALL_LOCATION}/*/ext)
foreach(dir ${dirs})
  message("Include ${dir}")
  include_directories(${dir})
  file(GLOB subdirs ${dir}/*/ext)
  foreach(subdir ${subdirs})
    message("Include ${subdir}")
    include_directories(${subdir})
  endforeach()
endforeach()

