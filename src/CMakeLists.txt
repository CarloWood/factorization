include(AICxxProject)

add_executable(plot plot.cpp)
target_include_directories(plot PRIVATE ${Python_INCLUDE_DIRS} ${NUMPY_INCLUDE})
target_link_libraries(plot PRIVATE ${AICXX_OBJECTS_LIST} ${Python_LIBRARIES})
