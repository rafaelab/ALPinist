cmake_minimum_required(VERSION 2.8.12)

set(PLUGINNAME "alpinist")
project(${PLUGINNAME} C CXX)

find_package(PythonInterp)
find_package(PythonLibs)
find_package(SWIG REQUIRED)


# require C++11
# Note: CRPropa requires C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")


# SWIG options
option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" OFF) 
set(crpropaSwigFile "crpropa.i")
set(swigModeFlag "")
if(ENABLE_SWIG_BUILTIN)
  set(crpropaSwigFile "crpropa-builtin.i")
  set(swigModeFlag "-builtin")
endif(ENABLE_SWIG_BUILTIN)

# NumPy
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/checkNumpy.py" OUTPUT_VARIABLE numpyIncludePath)
  if(numpyIncludePath)
    MESSAGE(STATUS "Found NumPy headers in: " ${numpyIncludePath})
    SET(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
    list(APPEND CRPROPA_SWIG_DEFINES -DWITHNUMPY)
    include_directories(${numpyIncludePath})
  else(numpyIncludePath)
    MESSAGE(STATUS "NumPy not found.")
  endif(numpyIncludePath)

# get information on CRPropa
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" swig_interface OUTPUT_VARIABLE CRPROPA_SWIG_PATH)
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" install_prefix OUTPUT_VARIABLE CRPROPA_INSTALL_PREFIX)

find_path(CRPROPA_SWIG_INTERFACE_FILES ${crpropaSwigFile} HINTS crpropa/ share/crpropa python/ ${CRPROPA_SWIG_PATH})

find_path(CRPROPA_INCLUDES CRPropa.h HINTS crpropa/ include/crpropa ${CRPROPA_INSTALL_PREFIX}/include/)

find_library(CRPROPA_LIB NAMES crpropa libcrpropa HINTS crpropa/ lib/crpropa crpropa/lib ${CRPROPA_INSTALL_PREFIX}/lib/)

message(STATUS "CRPropa SWIG interface files: ${CRPROPA_SWIG_INTERFACE_FILES}")
message(STATUS "CRPropa includes: ${CRPROPA_INCLUDES}")
message(STATUS "CRPropa lib: ${CRPROPA_LIB}")


# kiss (provided by CRPropa)
include_directories(include ${crpropaInstallPrefix}/libs/kiss)

# Eigen library
include_directories(include libs/eigen3)



# SWIG plugin
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx PROPERTIES GENERATED true)
add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  COMMAND swig -builtin -c++ -python ${PY3} -I${CRPROPA_SWIG_INTERFACE_FILES} -I${CRPROPA_INCLUDES} -I${EIGEN_INCLUDES} -I${CMAKE_SOURCE_DIR}/include -I${PYTHON_INCLUDE_DIRS} ${CRPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/python/${PLUGINNAME}.i)

include_directories(include ${CMAKE_CURRENT_BINARY_DIR}/include ${CRPROPA_INCLUDES} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

# build interface and source in one library
add_library(${PLUGINNAME} MODULE 
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PlasmaDensity.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/SourceFeatures.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ALPPhotonMixing.cc
  )

target_link_libraries(${PLUGINNAME} ${CRPROPA_LIB} ${PYTHON_LIBRARIES})

set_target_properties(${PLUGINNAME} PROPERTIES PREFIX "")
set_target_properties(${PLUGINNAME} PROPERTIES OUTPUT_NAME "_${PLUGINNAME}")

execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

# install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}.py" DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")
# install(TARGETS ${PLUGINNAME} LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")

# ----------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------

option(ENABLE_TESTING "Build tests and enable test target" ON)
if(ENABLE_TESTING)
  set(GTEST_INCLUDE_PATH ${CRPROPA_INSTALL_PREFIX}/../libs/gtest/include)
  include_directories(${GTEST_INCLUDE_PATH})
  link_directories(${CRPROPA_INSTALL_PREFIX}/libs/gtest)
  message(STATUS "Google testing suite (from CRPropa): ${CRPROPA_INSTALL_PREFIX}/libs/gtest/include")
  if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
  endif(APPLE)


  enable_testing()

#   add_executable(testVector test/testVector.cpp)
#   target_link_libraries(testVector ${CRPROPA_LIB} gtest gtest_main pthread)
#   add_test(testVector testVector)

#   add_executable(testMatrix test/testMatrix.cpp)
#   target_link_libraries(testMatrix ${CRPROPA_LIB} gtest gtest_main pthread)
#   add_test(testMatrix testMatrix)

#   add_executable(testElectromagneticWave test/testElectromagneticWave.cpp src/ElectromagneticWave.cc)
#   target_link_libraries(testElectromagneticWave ${CRPROPA_LIB} gtest gtest_main pthread)
#   add_test(testElectromagneticWave testElectromagneticWave)

#   add_executable(testSourceFeatures test/testSourceFeatures.cpp src/SourceFeatures.cc src/ElectromagneticWave.cc)
#   target_link_libraries(testSourceFeatures  ${CRPROPA_LIB} gtest gtest_main pthread)
#   add_test(testSourceFeatures testSourceFeatures)

  #add_executable(testMixing test/testMixing.cpp src/SourceFeatures.cc src/ElectromagneticWave.cc src/AFields.cc)
  #target_link_libraries(testMixing  ${CRPROPA_LIB} gtest gtest_main pthread)
  #add_test(testMixing testMixing)

endif(ENABLE_TESTING)

# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
# install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}.py" DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")
# install(TARGETS ${PLUGINNAME} DESTINATION lib)
# install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
# install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}.py" DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")
install(TARGETS ${PLUGINNAME} LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")
