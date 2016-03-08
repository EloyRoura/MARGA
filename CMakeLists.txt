PROJECT(MARGA)

SET(CMAKE_BUILD_TYPE "Release")
cmake_minimum_required(VERSION 2.6)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Please Set ITK_DIR")
ENDIF(ITK_FOUND)

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_EXECUTABLE(Marga main.cpp brainio.cpp brainExtraction.cpp imageDefinitions.h brainExtraction.h brainio.h)
#TARGET_LINK_LIBRARIES(Marga ITKCommon ITKIO ITKStatistics ITKAlgorithms ITKNumerics ITKBasicFilters)
TARGET_LINK_LIBRARIES( Marga ${ITK_LIBRARIES})

ADD_EXECUTABLE(fillInHoles fillInHoles.cxx brainio.cpp brainio.h imageDefinitions.h)
#TARGET_LINK_LIBRARIES(fillInHoles ITKCommon ITKIO ITKStatistics ITKAlgorithms ITKNumerics ITKBasicFilters)
TARGET_LINK_LIBRARIES( fillInHoles ${ITK_LIBRARIES})
