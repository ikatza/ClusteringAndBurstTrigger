# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(DAQWork)
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net XMLIO)
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# ROOT_GENERATE_DICTIONARY(G__OpticalHit      OpticalHit.hh)
# ROOT_GENERATE_DICTIONARY(G__WireHit         WireHit.hh)
# ROOT_GENERATE_DICTIONARY(G__OpticalCluster  OpticalCluster.hh)
# ROOT_GENERATE_DICTIONARY(G__WireCluster     WireCluster.hh)
# ROOT_GENERATE_DICTIONARY(G__ClusterSelector ClusterSelector.hh)
# ROOT_GENERATE_DICTIONARY(G__ClusterEngine   ClusterEngine.hh)
# ROOT_GENERATE_DICTIONARY(G__Trigger         Trigger.hh)

add_library(Clustering SHARED WireCluster.cc WireHit.cc OpticalCluster.cc OpticalHit.cc ClusterEngine.cc ClusterEnergyEstimator.cc Clustering.cc SmartTrigger.cc  SmartTriggerToyThrower.cc)
target_link_libraries(Clustering ${ROOT_LIBRARIES})

add_library(Display SHARED PhotonTimingDisplay.cc PhotonWirePositionDisplay.cc WirePositionTimingDisplay.cc)
target_link_libraries(Display ${ROOT_LIBRARIES})

add_library(SimpleOptHit SHARED HitSimplePlot.cc)
target_link_libraries(SimpleOptHit ${ROOT_LIBRARIES})

link_directories(build)

add_executable(RunDAQClustering RunDAQClustering.cxx)
target_link_libraries(RunDAQClustering Clustering)

add_executable(AnalyseOpticalCluster AnalyseOpticalCluster.cxx)
target_link_libraries(AnalyseOpticalCluster Clustering)

add_executable(GetEffBackRate GetEffBackRate.cxx)
target_link_libraries(GetEffBackRate Clustering)

add_executable(RunEventDiplay RunEventDisplay.cxx)
target_link_libraries(RunEventDiplay Display)

add_executable(AnalyseDAQClustering AnalyseDAQClustering.cxx)
target_link_libraries(AnalyseDAQClustering Clustering)

add_executable(RunSimpleHitInfo RunSimpleHitInfo.cxx)
target_link_libraries(RunSimpleHitInfo SimpleOptHit)

add_executable(PDSTiming PDSTiming.cxx)
target_link_libraries(PDSTiming ${ROOT_LIBRARIES})

add_executable(NeutronInfo NeutronInfo.cxx)
target_link_libraries(NeutronInfo ${ROOT_LIBRARIES})

add_executable(RunSmartTrigger RunSmartTrigger.cxx)
target_link_libraries(RunSmartTrigger Clustering)

add_executable(ThrowToySmartTrigger ThrowToySmartTrigger.cxx)
target_link_libraries(ThrowToySmartTrigger Clustering)

# add_executable(OpHitTrigger OpHitTrigger.C)
# target_link_libraries(OpHitTrigger ${ROOT_LIBRARIES})
