// GraphComplex.h
#ifndef CHOMP_GRAPHCOMPLEX_H
#define CHOMP_GRAPHCOMPLEX_H

#include "chomp/Complex.h"
#include <boost/foreach.hpp>

namespace chomp {
  /**************************
   *      GRAPH CELL        *
   **************************/
  class GraphCell {
  public:
    Index domain_index;
    Index codomain_index;
    int domain_dim;
    int codomain_dim;
    GraphCell ( void ) {}
    GraphCell ( Index dc, int dd, Index cc, int cd ) :
    domain_index ( dc ), codomain_index (cc), domain_dim (dd), codomain_dim (cd) {}
  };
  
  bool operator == (const GraphCell & lhs,
                    const GraphCell & rhs)
  {
    if ( lhs . domain_index != rhs . domain_index ) return false;
    if ( lhs . codomain_index != rhs . codomain_index ) return false;
    if ( lhs . domain_dim != rhs . domain_dim ) return false;
    if ( lhs . codomain_dim != rhs . codomain_dim ) return false;
    return true;
  }
  
  std::size_t hash_value ( const GraphCell & gc )
  {
    std::size_t seed = 0;
    boost::hash_combine(seed, gc . domain_index );
    boost::hash_combine(seed, gc . codomain_index );
    boost::hash_combine(seed, gc . domain_dim );
    boost::hash_combine(seed, gc . codomain_dim );
    return seed;
  }
  
  /*******************************
   *      GRAPH COMPLEX        *
   *******************************/
  class GraphComplex : public Complex {
  public:
    // Cell type: GraphCell
    CHOMP_COMPLEX(GraphCell)
    /*******************************
     *      COMPLEX INTERFACE      *
     *******************************/
    /// boundary. See Complex.h 
    virtual void boundary ( Chain * output, const Index input, int dim ) const;
    /// coboundary. See Complex.h
    virtual void coboundary ( Chain * output, const Index input, int dim ) const;
    
    /*****************************
     * SPECIFIC TO GRAPH COMPLEX *
     *****************************/
    GraphComplex (Complex * domain,
                  Complex * codomain);
    
    const Complex & domain ( void ) const;
    
    const Complex & codomain ( void ) const;
    
    void insert ( const GraphCell & gc );
    
    void insert ( Index dc, int dd, Index cc, int cd );

    void projectToDomain (Chain * output, 
                          const Chain & input ) const;
    
    void projectToCodomain (Chain * output, 
                            const Chain & input ) const;
    
    Chain projectToDomain ( const Chain & input ) const {
      Chain output; projectToDomain ( &output, input ); return output; 
    }
    
    Chain projectToCodomain ( const Chain & input ) const {
      Chain output; projectToCodomain ( &output, input ); return output; 
    }
  private:
    Complex * domain_;
    Complex * codomain_;
    
  };
  
  // DEFINITIONS
  
  GraphComplex::GraphComplex (Complex * domain,
                              Complex * codomain) : domain_ (domain), codomain_ (codomain) {
    
  }
  
  const Complex & GraphComplex::domain ( void ) const {
    return * domain_;
  }
  
  const Complex & GraphComplex::codomain ( void ) const {
    return * codomain_;
  }
  
  void GraphComplex::boundary (Chain * output, 
                                       const Index input, 
                                       int dim ) const {
    output -> dimension () = dim - 1;
    Cell c = indexToCell ( input, dim );
    Chain dom_bd = domain () . boundary ( c . domain_index, c . domain_dim );
    Chain cod_bd = codomain () . boundary ( c . codomain_index, c . codomain_dim );
    Cell bd_cell;
    // domain bd x codomain original terms
    bd_cell . codomain_index = c . codomain_index;
    bd_cell . codomain_dim = c . codomain_dim;
    bd_cell . domain_dim = c . domain_dim - 1; 
    BOOST_FOREACH ( const Term & t, dom_bd () ) {
      bd_cell . domain_index = t . index ();
      Index i = cellToIndex ( bd_cell, dim - 1 );
      if ( i != size ( dim - 1 ) ) { /* i == size(dim-1) signals invalid cell */
        *output += Term ( i, t . coef () );
      }
    }
    // Determine sign of domain original x codomain bd terms
    Ring sign;
    if ( c . domain_dim % 2 == 0 ) {
      sign = Ring ( 1 ); 
    } else { 
      sign = Ring ( -1 );
    }
    // domain original x codomain bd terms
    bd_cell . domain_index = c . domain_index;
    bd_cell . domain_dim = c . domain_dim;
    bd_cell . codomain_dim = c . codomain_dim - 1; 
    BOOST_FOREACH ( const Term & t, cod_bd () ) {
      bd_cell . codomain_index = t . index ();
      Index i = cellToIndex ( bd_cell, dim - 1 );
      if ( i != size ( dim - 1 ) ) { /* i == size(dim-1) signals invalid cell */
        *output += Term ( i, sign * t . coef () );
      }
    }
    
  }

  void GraphComplex::coboundary (Chain * output, 
                                         const Index input, 
                                         int dim ) const {
    output -> dimension () = dim + 1;
    Cell c = indexToCell ( input, dim );
    Chain dom_cbd = domain () . coboundary ( c . domain_index, c . domain_dim );
    Chain cod_cbd = codomain () . coboundary ( c . codomain_index, c . codomain_dim );
    Cell cbd_cell;
    // domain bd x codomain original terms
    cbd_cell . codomain_index = c . codomain_index;
    cbd_cell . codomain_dim = c . codomain_dim;
    cbd_cell . domain_dim = c . domain_dim + 1; 
    BOOST_FOREACH ( const Term & t, dom_cbd () ) {
      cbd_cell . domain_index = t . index ();
      Index i = cellToIndex ( cbd_cell, dim + 1 );
      if ( i != size ( dim + 1 ) ) { /* i == size(dim+1) signals invalid cell */
        *output += Term ( i, t . coef () );
      }
    }
    // Determine sign of domain original x codomain bd terms
    Ring sign;
    if ( c . domain_dim % 2 == 0 ) {
      sign = Ring ( 1 ); 
    } else { 
      sign = Ring ( -1 );
    }
    // domain original x codomain bd terms
    cbd_cell . domain_index = c . domain_index;
    cbd_cell . domain_dim = c . domain_dim;
    cbd_cell . codomain_dim = c . codomain_dim + 1; 
    BOOST_FOREACH ( const Term & t, cod_cbd () ) {
      cbd_cell . codomain_index = t . index ();
      Index i = cellToIndex ( cbd_cell, dim + 1 );
      if ( i != size ( dim + 1 ) ) { /* i == size(dim+1) signals invalid cell */
        * output += Term ( i, sign * t . coef () );
      }
    }
    
    
  }
  
  void GraphComplex::insert ( const GraphCell & gc ) {
    int dim = gc . domain_dim + gc . codomain_dim;
    insertCell ( gc, dim );
  }
  
  inline void GraphComplex::insert ( Index dc, int dd, Index cc, int cd ) {
    GraphCell gc ( dc, dd, cc, cd );
    insert ( gc );
  }
  
  inline void GraphComplex::projectToDomain (Chain * output, 
                                      const Chain & input ) const {
    int D = output -> dimension () = input . dimension ();
    BOOST_FOREACH ( const Term & t, input () ) {
      Cell c = indexToCell ( t . index (), input . dimension () );
      if ( c . domain_dim == D ) * output += Term ( c . domain_index, t . coef () ); 
    }
    *output = simplify ( *output ); 
  }
  

  
  inline void GraphComplex::projectToCodomain (Chain * output, 
                                        const Chain & input ) const {
    int D = output -> dimension () = input . dimension ();
    BOOST_FOREACH ( const Term & t, input () ) {
      Cell c = indexToCell ( t . index (), input . dimension () );
      if ( c . codomain_dim == D ) * output += Term ( c . codomain_index, t . coef () ); 
    }
    *output = simplify ( *output );    
  }
  
} // namespace chomp

#endif
