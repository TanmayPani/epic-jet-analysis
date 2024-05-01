cmake_minimum_required(VERSION 3.5)
project(read_events)

set(ROOT_DIR, "/usr/local")
set(FASTJET_DIR, "/usr/local")

find_package(EDM4HEP 0.10.3 REQUIRED)
find_package(podio 0.17.3 REQUIRED)
find_package(ROOT 6.30 REQUIRED)

find_library(FASTJET_LIB libfastjet.so)

add_executable(${PROJECT_NAME} edm4eicDataReader.cpp edm4eicDataReader.hh)
separate_arguments(ROOT_EXE_LINKER_FLAGS)
target_link_libraries(${PROJECT_NAME} PUBLIC edm4hepRDF podio::podioRootIO  
                                             ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS} 
                                             ${FASTJET_LIB})
#mkdir build; cd build; cmake ..; make; ./read_events