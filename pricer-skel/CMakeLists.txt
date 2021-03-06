cmake_minimum_required(VERSION 2.8)
project(mc-pricer CXX)

# Prefix
if (NOT PREFIX)
    set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()


# Detect Boost
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
message(STATUS "Boost include directory ${Boost_INCLUDE_DIR}")

# Detect PNL
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()

# Add math lib
set(LIBS ${LIBS} m)

# Testing the parser
add_executable(parser-test src/parser.cpp src/test_parser.cpp)
target_link_libraries(parser-test ${LIBS})
add_executable(pricer src/BlackScholesModel.cpp src/MonteCarlo.cpp src/BasketOption.cpp src/AsianOption.cpp src/PerformanceOption.cpp src/parser.cpp src/Market.cpp src/pricer.cpp)
target_link_libraries(pricer ${LIBS})

#Testing pricing different options at t=0
add_executable(testInitPricingOptions src/Market.cpp src/BlackScholesModel.cpp src/MonteCarlo.cpp src/BasketOption.cpp src/AsianOption.cpp src/PerformanceOption.cpp src/parser.cpp src/testInitPricingOptions.cpp)
target_link_libraries(testInitPricingOptions ${LIBS})

#Testing deltas calculation for different options at t=0
add_executable(testInitiDelta src/Market.cpp src/BlackScholesModel.cpp src/MonteCarlo.cpp src/BasketOption.cpp src/AsianOption.cpp src/PerformanceOption.cpp src/parser.cpp src/testInitiDelta.cpp)
target_link_libraries(testInitiDelta ${LIBS})

#Testing delta calculation at t >= 0
add_executable(testPastBasketDelta src/Market.cpp src/BlackScholesModel.cpp src/MonteCarlo.cpp src/BasketOption.cpp src/AsianOption.cpp src/PerformanceOption.cpp src/parser.cpp src/testPastBasketDelta.cpp)
target_link_libraries(testPastBasketDelta ${LIBS})

#Testing delta and ic calculation at t >= 0
add_executable(testIcDelta src/Market.cpp src/BlackScholesModel.cpp src/MonteCarlo.cpp src/BasketOption.cpp src/AsianOption.cpp src/PerformanceOption.cpp src/parser.cpp src/testIcDelta.cpp)
target_link_libraries(testIcDelta ${LIBS})

#Testing pricing at t > 0 using a vanilla
add_executable(testPastPricingVanilla src/Market.cpp src/BlackScholesModel.cpp src/MonteCarlo.cpp src/BasketOption.cpp src/AsianOption.cpp src/PerformanceOption.cpp src/parser.cpp src/testPastPricingVanilla.cpp)
target_link_libraries(testPastPricingVanilla ${LIBS})


# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_LINKER_EXE_FLAGS}" )
else()
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_LINKER_EXE_FLAGS}" )
endif()
