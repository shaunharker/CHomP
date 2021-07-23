
set ( LIBS ${LIBS}
           ${CMAKE_THREAD_LIBS_INIT}
           ${Boost_LIBRARIES} )

set ( TARGETS chomp-cubical
              chomp-simplicial
	            chomp-simplicial-z2
	            chomp-simplicial-z3)

foreach ( TARGET ${TARGETS} )
  add_executable ( ${TARGET} ${TARGET}.cpp )
  target_link_libraries ( ${TARGET} ${LIBS} )
endforeach ( TARGET )

install( TARGETS ${TARGETS}
         RUNTIME DESTINATION bin )
