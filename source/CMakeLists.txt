
set ( LIBS ${LIBS}
           ${X11_LIBRARIES}
           ${CMAKE_THREAD_LIBS_INIT}
           ${Boost_LIBRARIES} )

set ( TARGETS chomp-cubical 
              chomp-matrix 
              chomp-simplicial 
              chomp-greyscale-to-cubical)

foreach ( TARGET ${TARGETS} ) 
  add_executable ( ${TARGET} ${TARGET}.cpp )
  target_link_libraries ( ${TARGET} ${LIBS} ) 
endforeach ( TARGET ) 

install( TARGETS ${TARGETS}
         RUNTIME DESTINATION bin )
