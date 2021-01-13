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

  if(NOT _dir STREQUAL CMAKE_SOURCE_DIR)
    string(REPLACE "${CMAKE_SOURCE_DIR}/" "" rel_dir "${_dir}")
  endif()

  get_property(_targets_list DIRECTORY "${_dir}" PROPERTY
    BUILDSYSTEM_TARGETS)

  foreach(_tgt ${_targets_list})
    get_property(_t_files  TARGET ${_tgt} PROPERTY SOURCES)
    list(FILTER _t_files EXCLUDE REGEX TARGET_OBJECTS)
    
    foreach(my_file IN LISTS _t_files)
      if(my_file MATCHES "${CMAKE_SOURCE_DIR}/(.*)")
	list(APPEND _all_source_files ${CMAKE_MATCH_1})
      else()
	get_filename_component(my_dir ${my_file} DIRECTORY)
	
	if(my_dir OR _dir STREQUAL CMAKE_SOURCE_DIR)
	  list(APPEND _all_source_files ${my_file})
	else()
	  list(APPEND _all_source_files "${rel_dir}/${my_file}")
	endif()
      endif()
    endforeach()
  endforeach()

  set(_all_source_files ${_all_source_files} PARENT_SCOPE)
endfunction()

if(CTAGS)
  ##variable_watch(_all_source_files)
  set(SOURCES_LIST "${CMAKE_BINARY_DIR}/sources_list")
  set(_all_source_files "")
  file(REMOVE "${SOURCES_LIST}")
  tags_add_dir(${CMAKE_SOURCE_DIR})
  list(SORT _all_source_files)
  list(REMOVE_DUPLICATES _all_source_files)
  string(REGEX REPLACE ";" "\n" _asf "${_all_source_files}")
  file(APPEND ${SOURCES_LIST} "${_asf}\n")
  add_custom_target(TAGS DEPENDS ${CMAKE_SOURCE_DIR}/TAGS)

  add_custom_command(OUTPUT ${CMAKE_SOURCE_DIR}/TAGS
    COMMAND ctags -e --language-force=fortran -L ${SOURCES_LIST}
    DEPENDS ${_all_source_files}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM)
endif()
