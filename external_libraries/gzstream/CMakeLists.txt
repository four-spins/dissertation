add_library(gzstream gzstream.C)
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( gzstream ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

