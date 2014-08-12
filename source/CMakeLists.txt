# add the executable

set ( LIBS      ${LIBS}
                boost_system
                boost_thread
                boost_serialization )

set ( TARGETS chomp-cubical 
              chomp-matrix 
              chomp-simplicial 
              chomp-greyscale-to-cubical)

foreach ( TARGET ${TARGETS} ) 
  add_executable ( ${TARGET} ${TARGET}.cpp ${HEADERS} )
  target_link_libraries ( ${TARGET} ${LIBS} ) 
endforeach ( TARGET ) 

install( TARGETS ${TARGETS}
         RUNTIME DESTINATION bin )
