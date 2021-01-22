cmake_minimum_required (VERSION 3.0)
project( capsule_contact_manifold_qr )

# https://eigen.tuxfamily.org/dox/TopicCMakeGuide.html
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

file(GLOB header ${CMAKE_CURRENT_SOURCE_DIR}/include/capsule_contact_manifold_qr.h)
install(FILES ${header} DESTINATION include)

add_library(capsule_contact_manifold_qr INTERFACE)
target_sources(capsule_contact_manifold_qr INTERFACE ${header})
target_include_directories(capsule_contact_manifold_qr INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(capsule_contact_manifold_qr INTERFACE Eigen3::Eigen)
target_compile_features(capsule_contact_manifold_qr INTERFACE cxx_std_17)

IF(WIN32)
	target_compile_definitions(capsule_contact_manifold_qr INTERFACE NOMINMAX _USE_MATH_DEFINES)
ENDIF(WIN32)

option(CAPSULE_BUILD_TEST "build test" ON)
if(CAPSULE_BUILD_TEST)
	add_executable(capsule_test main.cpp)
	target_link_libraries(capsule_test capsule_contact_manifold_qr)
endif()
