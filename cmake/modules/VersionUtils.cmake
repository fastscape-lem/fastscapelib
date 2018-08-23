# Some utilities for handling versions
# (mostly inspired by Python versioneer).

include(GetGitRevisionDescription)

function(get_git_version_pieces _version_pieces)
  get_git_head_revision(_ _hash_full)
  git_describe(GIT_REV_DESCRIPTION --tags --always --dirty --long)

  if(${GIT_REV_DESCRIPTION} MATCHES
      "^[v]*([0-9\\.]+)\\-([0-9]+)\\-g([A-Za-z0-9]+)\\-*([dirty]*).*")
    # found tag + additional commits
    set(_version_tag ${CMAKE_MATCH_1})
    set(_commit_count ${CMAKE_MATCH_2})
    set(_hash_short ${CMAKE_MATCH_3})
    set(_local_changes ${CMAKE_MATCH_4})

  elseif(${GIT_REV_DESCRIPTION} MATCHES
      "^[v]*([0-9]+\\.[0-9]+\\.[0-9]+)\\-*([dirty]*).*")
    # HEAD is a tag (normally not needed if --long is used for git describe)
    set(_version_tag ${CMAKE_MATCH_1})
    set(_commit_count "")
    set(_hash_short "")
    set(_local_changes ${CMAKE_MATCH_2})

  elseif(${GIT_REV_DESCRIPTION} MATCHES
      "^([A-Za-z0-9]+)\\-*([dirty]*).*")
    # no tag found
    set(_version_tag "0")
    git_get_commit_count(_commit_count)
    set(_hash_short ${CMAKE_MATCH_1})
    set(_local_changes ${CMAKE_MATCH_2})

  else()
    message(FATAL_ERROR
      "Could not get version info correctly. Output of git describe: "
      ${GIT_REV_DESCRIPTION})
  endif()

  set(${_version_pieces}
    "${_version_tag};${_commit_count};${_hash_full};${_hash_short};${_local_changes}"
    PARENT_SCOPE)
endfunction()


function(get_source_dir_version_pieces _version_pieces)
  get_filename_component(_dirname ${CMAKE_CURRENT_SOURCE_DIR} NAME)

  if(${_dirname} MATCHES "^.*\\-[v]*([0-9\\.]+)")
    set(_version_tag ${CMAKE_MATCH_1})
  else()
    set(_version_tag "VERSION-NOTFOUND")
  endif()

  set(${_version_pieces} "${_version_tag};0;NOHASH;;" PARENT_SCOPE)
endfunction()


function(get_version_pieces _version_pieces _version_tag)
  if(_version_tag)
    set(_pieces "${_version_tag};0;NOHASH;;")

  else()
    # TODO: proper version/error handling when git is not installed
    get_git_head_revision(_ _hash_full)

    if(_hash_full STREQUAL "GITDIR-NOTFOUND")
      get_source_dir_version_pieces(_pieces)
    else()
      get_git_version_pieces(_pieces)
    endif()

  endif()

  set(${_version_pieces} "${_pieces}" PARENT_SCOPE)
endfunction()


function(format_version_pep440 _version_str _version_pieces)
  list(GET ${_version_pieces} 0 _version_tag)
  list(GET ${_version_pieces} 1 _commit_count)
  list(GET ${_version_pieces} 3 _hash_short)
  list(GET ${_version_pieces} 4 _local_changes)

  set(_local_version_items "")

  if(_version_tag EQUAL "0")
    list(APPEND _local_version_items "untagged")
  endif()
  if(_commit_count)
    list(APPEND _local_version_items ${_commit_count})
    if(_hash_short)
      list(APPEND _local_version_items "g${_hash_short}")
    endif()
  endif()
  if(_local_changes)
    list(APPEND _local_version_items ${_local_changes})
  endif()

  if(_local_version_items)
    string(REPLACE ";" "." _items_str "${_local_version_items}")
    set(_local_version_label "+${_items_str}")
  else()
    set(_local_version_label "")
  endif()

  set(${_version_str} ${_version_tag}${_local_version_label} PARENT_SCOPE)
endfunction()


function(get_version_numbers
    _version_major _version_minor _version_patch
    _version_pieces)
  list(GET ${_version_pieces} 0 _version_tag)

  if(${_version_tag} MATCHES "([0-9]+)\\.([0-9]+)\\.([0-9]+)")
    set(${_version_major} ${CMAKE_MATCH_1} PARENT_SCOPE)
    set(${_version_minor} ${CMAKE_MATCH_2} PARENT_SCOPE)
    set(${_version_patch} ${CMAKE_MATCH_3} PARENT_SCOPE)
  else()
    set(${_version_major} "0" PARENT_SCOPE)
    set(${_version_minor} "0" PARENT_SCOPE)
    set(${_version_patch} "0" PARENT_SCOPE)
  endif()
endfunction()
