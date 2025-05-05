project(masking)
cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

include(External)
include(Paths)
include(Compiler)

add_library(svgplot ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/svgplot.cc ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/plot.cc ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/scatter.cc ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/imshow.cc ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/bar.cc ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/barh.cc ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/hist.cc ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/axis-scale.cc ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/color.cc ${EXTERNAL_INSTALL_LOCATION}/svg.cc/plot/src/colormap.cc)


###########################################################################################
# TARGETS
###########################################################################################
add_executable(barbule-masking main/barbule-masking.cc)
add_dependencies(barbule-masking svg.cc)
add_executable(barb-masking main/barb-masking.cc)
add_dependencies(barb-masking svg.cc)

add_executable(plot-barbs main/plot-barbs.cc)
add_dependencies(plot-barbs svg.cc)
target_link_libraries(plot-barbs svgplot)

add_executable(plot-barbules main/plot-barbules.cc)
add_dependencies(plot-barbules svg.cc)
target_link_libraries(plot-barbules svgplot)

add_executable(htest main/htest.cc)
add_dependencies(htest svg.cc)
target_link_libraries(htest svgplot)

add_executable(test-projection main/test-projection.cc)

add_executable(plot-feathers main/plot-feathers.cc)
add_dependencies(plot-feathers svg.cc)
target_link_libraries(plot-feathers svgplot)

