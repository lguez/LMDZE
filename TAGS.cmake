#
# (C) Copyright 2018 Tillmann Heidsieck
#
# SPDX-License-Identifier: MIT
#

find_program(CTAGS ctags)

function(tags_add_dir _dir)
  get_property(_subdirlist DIRECTORY "${_dir}" PROPERTY SUBDIRECTORIES)
  
  foreach(_sd ${_subdirlist})
    tags_add_dir(${_sd})
  endforeach()

  get_property(_targets_list DIRECTORY "${_dir}" PROPERTY
    BUILDSYSTEM_TARGETS)

  foreach(_tgt ${_targets_list})
    get_property(_t_files  TARGET ${_tgt} PROPERTY SOURCES)
    list(FILTER _t_files EXCLUDE REGEX TARGET_OBJECTS)
    get_property(_t_defs   TARGET ${_tgt} PROPERTY
      COMPILE_DEFINITIONS)
    
    list(APPEND _all_source_files ${_t_files})
    string(REGEX REPLACE "${CMAKE_SOURCE_DIR}/" "" _t_files
      "${_t_files}")
    string(REGEX REPLACE ";" "\n" _t_files "${_t_files}")
    
    string(REGEX REPLACE ";" "\n" _t_defs "${_t_defs}")
    string(REGEX REPLACE "^([^=]+)$" "\\1=1" _t_defs "${_t_defs}")
  endforeach()

  set(_all_source_files ${_all_source_files} PARENT_SCOPE)
  set(_all_definitions ${_all_definitions} PARENT_SCOPE)
endfunction()

if(CTAGS)
  set(SOURCES_LIST "${CMAKE_BINARY_DIR}/sources_list")
  set(_all_source_files "")
  set(_all_definitions "")
  file(REMOVE "${SOURCES_LIST}")
  tags_add_dir(${CMAKE_SOURCE_DIR})
  list(SORT _all_source_files)
  list(REMOVE_DUPLICATES _all_source_files)

  string(REGEX REPLACE "${CMAKE_SOURCE_DIR}/" "" _asf
    "${_all_source_files}")

  string(REGEX REPLACE ";" "\n" _asf "${_asf}")
  file(APPEND ${SOURCES_LIST} "\n${_asf}")
  string(REGEX REPLACE ";" "\n" _all_definitions "${_all_definitions}")

  string(REGEX REPLACE "^([^=]+)$" "\\1=1" _all_definitions
    "${_all_definitions}")

  add_custom_target(TAGS DEPENDS ${CMAKE_SOURCE_DIR}/TAGS)

  add_custom_command(OUTPUT ${CMAKE_SOURCE_DIR}/TAGS
    COMMAND ctags -e --language-force=fortran -L ${SOURCES_LIST}
    DEPENDS ${_all_source_files}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM)
endif()
