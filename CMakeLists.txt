cmake_minimum_required(VERSION 3.14)
project(alpinist CXX)

message(STATUS "Building ALPinist...")


# set default build-type to release to enable performance improvements
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# compilation flags
set(CMAKE_CXX_FLAGS "-std=c++17")

# workaround for LLVM > 15
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15.0)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION")
	endif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15.0)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

# find CRPropa
include(cmake/FindCRPropa.cmake)

# note that CRPropa's eigen3 library SHOULD NOT be used here
# assemble includes
include_directories(include 
  libs/eigen3
  ${CRPropa_kiss_INCLUDE_DIR}
  ${CRPropa_HepPID_INCLUDE_DIR}
  ${CRPropa_INCLUDE_DIR}  
  ${CMAKE_CURRENT_BINARY_DIR}/include 
)

# source files for compilation
set(SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Data.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PlasmaDensity.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/SourceFeatures.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PolarisationCorrections.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PairProduction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/InverseComptonScattering.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ALPPhotonMixing.cc
  )

# decide on whether to use Python or not
option(ENABLE_PYTHON "Create python library via SWIG" On)
message(STATUS "Build Python bindings via SWIG: ${ENABLE_PYTHON}")


# compile with Python support
if(ENABLE_PYTHON)
	# find Python
	find_package(Python COMPONENTS Interpreter Development)

	if(Python_FOUND)
		message(STATUS "Python version: ${Python_VERSION}")
		message(STATUS "Python executable: ${Python_EXECUTABLE}")
		message(STATUS "Python headers: ${Python_INCLUDE_DIRS}")
		message(STATUS "Python library: ${Python_LIBRARY}")
		message(STATUS "Python site-packges: ${Python_SITE_PACKAGES}")
	else(Python_FOUND)
		set(Python_HELPER_EXECUTABLE ". Python_EXECUTABLE (path to python executable)")
		set(Python_HELPER_INCLUDE_DIRS ". Python_INCLUDE_DIRS (directory where Python.h is located)")
		set(Python_HELPER_LIBRARY ". Python_LIBRARY (path to file libpython.so or libpython.dylib [in OSX]; sometimes contains version before extension)")
		message(STATUS 	"Python: NOT FOUND. Try setting the relevant flags manually: \n\t${Python_HELPER_EXECUTABLE} \n\t${Python_HELPER_LIBRARY} \n\t${Python_HELPER_INCLUDE_DIRS}")
	endif(Python_FOUND)

	# find SWIG
	find_package(SWIG REQUIRED)

	# only support versions of SWIG higher than 4.0
	if(SWIG_VERSION VERSION_LESS 4.0)
		message(STATUS "SWIG version too old (${SWIG_VERSION} < 4.0). Code will not be compile with Python support.")
		set(SWIG_FOUND False)
	endif()

	if(SWIG_FOUND)
		message(STATUS "SWIG installation directory: ${SWIG_DIR}")
		message(STATUS "SWIG executable: ${SWIG_EXECUTABLE}")
	else(SWIG_FOUND)
		message(STATUS "SWIG: NOT FOUND (cannot build with Python bindings)")
	endif(SWIG_FOUND)

	# add Python headers
	include_directories(include ${Python_INCLUDE_DIRS})


	# SWIG plugin
	set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/alpinist_wrap.cxx PROPERTIES GENERATED True)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/alpinist_wrap.cxx
		COMMAND swig ${SWIG_MODE_FLAG} -c++ -python ${PY3} 
		-I${CRPropa_SWIG_PATH} 
		-I${CRPropa_SWIG_INTERFACE_FILE} 
		-I${CRPropa_INCLUDE_DIR} ${CRPropa_SWIG_DEFINES} 
		-I${CMAKE_SOURCE_DIR}/include 
		-I${Python_INCLUDE_DIRS} 
		-dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/alpinist_wrap.cxx 
		-outdir ${CMAKE_CURRENT_BINARY_DIR} 
		${CMAKE_CURRENT_SOURCE_DIR}/python/alpinist.i
	)
	
	# build interface and source in one library
	add_library(alpinist MODULE ${CMAKE_CURRENT_BINARY_DIR}/alpinist_wrap.cxx ${SOURCE_FILES})

	target_link_libraries(alpinist ${CRPropa_LIBRARY} ${Python_LIBRARIES})

# compile code to work in C++ only, without Python support
elseif()
	# build interface and source in one library
	add_library(alpinist MODULE ${SOURCE_FILES})
	target_link_libraries(alpinist ${CRPropa_LIBRARY})
endif()

set_target_properties(alpinist PROPERTIES PREFIX "")
set_target_properties(alpinist PROPERTIES COMPILE_FLAGS "-w")
set_target_properties(alpinist PROPERTIES OUTPUT_NAME "_alpinist")

option(ENABLE_TESTING "Build tests and enable test target" On)
if(ENABLE_TESTING)
	set(GTEST_INCLUDE_PATH ${CRPropa_INSTALL_PREFIX}/../libs/gtest/include)
	include_directories(${GTEST_INCLUDE_PATH})
	link_directories(${CRPropa_INSTALL_PREFIX}/libs/gtest)
	message(STATUS "Google testing suite (from CRPropa): ${CRPropa_INSTALL_PREFIX}/libs/gtest/include")
	enable_testing()
endif(ENABLE_TESTING)


# Install
add_definitions(-DALPinist_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
add_definitions(-DALPinist_DATA_PATH="${CMAKE_INSTALL_PREFIX}/data")
install(TARGETS alpinist DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_INSTALL_PREFIX}/.)
install(DIRECTORY ${CMAKE_BINARY_DIR}/data/ DESTINATION share/alpinist PATTERN ".git" EXCLUDE)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/alpinist.i DESTINATION swig_interface)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_INSTALL_PREFIX}/.)