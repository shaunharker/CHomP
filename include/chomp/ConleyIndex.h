/// ConleyIndex.h
/// Shaun Harker
/// 9/22/11

#ifndef CHOMP_CONLEYINDEX_H
#define CHOMP_CONLEYINDEX_H

#include "chomp/SparseMatrix.h"
#include "chomp/RelativeMapHomology.h"

/* Conley Index */
// Later separate this into a distinct file
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/nvp.hpp>

namespace chomp { 
  
class ConleyIndex_t {
private:
  std::vector < SparseMatrix < Ring > >  data_;
  bool undefined_;  
public:
  ConleyIndex_t ( void );
  std::vector < SparseMatrix < Ring > > & data ( void );
  const std::vector < SparseMatrix < Ring > > & data ( void ) const;
  bool & undefined ( void );
  const bool & undefined ( void ) const;
  
  /// The serialization method.
  friend class boost::serialization::access;
  template < class Archive >
  void serialize ( Archive & ar , const unsigned int version ) {
    ar & boost::serialization::make_nvp("data",data_);
    ar & boost::serialization::make_nvp("undefined",undefined_);
  }
};

template < class Grid, class Subset, class Map > void
ConleyIndex ( ConleyIndex_t * output,
		           const Grid & grid,
		           const Subset & S,
		      /* const */ Map & F ) {
  /* Method.
  1. Let F be the combinatorial map on "grid" that "f" induces
  2. Generate a combinatorial index pair (X, A)
     A combinatorial index pair satisfies the following:
      a) F(X\A) \subset X
      b) F(A) \cap X \subset A
     Since S is a Path-SCC (Combinatorial Morse set) in "grid" with respect to F,
     the following is a combinatorial index pair:
     X = S \cup F(S), A = F(S) \ S.
     proof: We show (a). Observe that X \ A = S, so F(X\A) = F(S) \subset X trivially. Hence (a). 
            We show (b). First we make the following observation from the definition of Path-SCC:
            Since S is a combinatorial Morse set, (t \in F(S) and s \in F(t) \cap S) implies t \in S.
            Suppose t \in A, s \in F(t), and s \in X. We show s \in A. Assume, to the contrary, that s \notin A.
            Then s \in X \ A = S. We then have t \in A \subset F(S), and s \in F(t) \cap S. Thus t \in S, by
            the observation. Yet t \in A \cap S = \emptyset is impossible! Hence (b).
   
  3. Generate the restricted combinatorial map G : (X, A) -> (X, A) via the formula
     G = F\vert_A \cap A
     That is, for each x \in X, define G(x) := F(x) \cap A.
     A theorem gives us F_* = G_*.
  4. Return the Relative Homology G_* of G
  */
  
  typedef typename Grid::GridElement Cell;
  typedef std::map < Cell, Subset > Combinatorial_Map;
  
  int depth = grid . getDepth ( S );
  
  //std::cout << "Conley Index. S . size () " << S . size () << "\n";
  //std::cout << "Conley Index. getDepth(S) returned " << depth << "\n";
  
  //clock_t start, start0, stop;
  //std::cout << "Conley Index. Preparing computation...\n";
  //start0 = start = clock ();
  
  /* Construct G on S, and also X and S in hash set forms called X_cells and S_cells */
  // note: All cells in X may be obtained as images of cells in S.
  typedef boost::unordered_set<Cell> CellDictionary;
  CellDictionary X_cells;
  CellDictionary S_cells;
  Combinatorial_Map G;
  std::insert_iterator<CellDictionary> Xc_ii ( X_cells, X_cells . begin () );
  BOOST_FOREACH ( Cell cell, S ) {
    S_cells . insert ( cell );
    std::vector < typename Grid::GridElement > image = 
      grid . cover ( F ( grid . geometry ( cell ) ) );
    std::copy ( image . begin (), image . end (), Xc_ii );
  } /* boost_foreach */
  /* Construct X and A */
  // note: The cells in A are those in X not in S.
  Subset X, A;
  std::insert_iterator<Subset> X_ii ( X, X . begin () );
  std::insert_iterator<Subset> A_ii ( A, A . begin () );
  BOOST_FOREACH ( Cell cell, X_cells ) {
    * X_ii ++ = cell;
    if ( S_cells . count ( cell ) == 0 ) {
      * A_ii ++ = cell;
      
    }
  } /* boost_foreach */ 
  
  // BEGIN DEBUG
  //std::cout << "Conley Index. getDepth(A) returned " << grid . getDepth ( A ) << "\n";
  // END DEBUG

  /* Compute G for domain cells in A */
#if 0
  // note: we restrict the ranges to A 
  BOOST_FOREACH ( Cell cell, A ) {
    Subset image;
    std::insert_iterator<Subset> image_ii ( image, image . end () );
    grid . cover ( image_ii, F ( grid . geometry ( grid . find ( cell ) ) ) );
    G [ cell ] = Subset ();
    std::insert_iterator<Subset> G_ii ( G [ cell ], G [ cell ] . begin () );
    BOOST_FOREACH ( Cell image_cell, image ) {
      if ( X_cells . count ( image_cell ) != 0 ) * G_ii ++ = image_cell; // better to look in A_cells, but i didn't make it
    }
    
  } /* boost_foreach */  
#endif
  
  //stop = clock ();
  //std::cout << "Conley Index computation prepared as relative map homology problem.\n";
  //std::cout << "Elapsed time = " << (float) ( stop - start ) / (float) CLOCKS_PER_SEC << "\n";
  //start = clock ();
  
  std::cout << "ConleyIndex: calling RelativeMapHomology.\n";
  int error_code = RelativeMapHomology ( &(output -> data ()), grid, X, A, grid, X, A, F, depth );
  if ( error_code == 1 ) {
    std::cout << "Problem computing conley index. Returning undefined result.\n";
    output -> undefined () = true;
    return;
  }
  
  //stop = clock ();  
  //std::cout << "Conley Index computed. Total time = " << (float) ( stop - start0 ) / (float) CLOCKS_PER_SEC << "\n";
  //std::cout << "Number of cells in X = " << X . size () << "\n Cells per second = " << (float) X . size () * (float) CLOCKS_PER_SEC / (float) ( stop - start0 )  << "\n";
  
  return;
} /* void ConleyIndex(...) */

  
  template < class Grid, class Subset, class Map > void
  ConleyIndexODE ( ConleyIndex_t * output,
                  const Grid & grid,
                  const Subset & S,
                  /* const */ Map & F ) {
    
    
    typedef typename Grid::GridElement Cell;
    
    int depth = grid . getDepth ( S );
    
    clock_t total_time_start = clock ();
    
    
    // note: All cells in X may be obtained as images of cells in S.
    typedef boost::unordered_set<Cell> CellDictionary;
    CellDictionary X_cells; // The cells in the total index pair (X, A)
    CellDictionary S_cells; // The cells just in the morse set S
    
    PRINT "ConleyIndexODE: Computing forward image of Morse Set S to find exit set A \n";

    std::insert_iterator<CellDictionary> Xc_ii ( X_cells, X_cells . begin () );
    BOOST_FOREACH ( Cell cell, S ) {
      S_cells . insert ( cell );
      Subset image = Subset ();
      std::insert_iterator<Subset> image_ii ( image, image . end () );
      grid . cover ( image_ii, F ( grid . geometry ( grid . find ( cell ) ) ) );
      std::copy ( image . begin (), image . end (), Xc_ii );
    } /* boost_foreach */
    
    PRINT "ConleyIndexODE: Constructing Index Pair (X, A) (note: X = S \\cup A) \n";

    /* Construct X and A */
    // note: The cells in A are those in X not in S.
    Subset X, A;
    std::insert_iterator<Subset> X_ii ( X, X . begin () );
    std::insert_iterator<Subset> A_ii ( A, A . begin () );
    BOOST_FOREACH ( Cell cell, X_cells ) {
      * X_ii ++ = cell;
      if ( S_cells . count ( cell ) == 0 ) {
        * A_ii ++ = cell;
        
      }
    } /* boost_foreach */
    
    
    /// NOW WE HAVE grid, X, A, F, depth
    
    
    
    
    // Create the Relative Complexes
    RelativePair domain_pair;
    // Determine depth.
    // TODO: determine it, don't ask for it as input
    grid . relativeComplex ( &domain_pair, X, A, depth );
    
    PRINT "ConleyIndexODE: Constructing Relative Complex Of Morse Set and Exit Set (X, A) \n";
    PRINT "ConleyIndexODE: Note: X . size () = " << X . size () << " and A . size () = " << A . size () << "\n";
    
    CubicalComplex & full_domain = domain_pair . base ();
    BitmapSubcomplex & domain = domain_pair . pair ();
    BitmapSubcomplex & exit = domain_pair . relative ();

    int D = full_domain . dimension ();
    
    PRINT "ConleyIndexODE: Size of Relative Complex = " << full_domain . size () << "\n";
    PRINT "ConleyIndexODE: Dimension of Relative Complex = " << D << "\n";
    
    int cutoff_dimension = full_domain . dimension ();
#ifdef CONLEYINDEXCUTOFF
    cutoff_dimension = CONLEYINDEXCUTOFF;
#endif
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Compute the homology generators (X, A)
    PRINT "ConleyIndexODE: Computing Relative Homology H_*(X, A)\n";
    Generators_t domain_gen = MorseGenerators ( domain, cutoff_dimension );
    for ( int d = 0; d <= domain . dimension (); ++ d ) {
#ifdef CONLEYINDEXCUTOFF
      if ( d == CONLEYINDEXCUTOFF ) break;
#endif
      int num_gen = domain_gen [ d ] . size ();
      Matrix MapHom ( num_gen, num_gen );
      for ( int i = 0; i < num_gen; ++ i ) MapHom . write ( i, i , 1 );
      std::cout << "ConleyIndexODE: Betti Number at Dimension " << d << ": " << num_gen << "\n";
      output -> data () . push_back ( MapHom );
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Compute the homology generators of exit set
    PRINT "ConleyIndexODE: Computing H_*(A)\n";
    Generators_t exit_gen = MorseGenerators ( exit, cutoff_dimension );
    for ( int d = 0; d <= exit . dimension (); ++ d ) {
#ifdef CONLEYINDEXCUTOFF
      if ( d == CONLEYINDEXCUTOFF ) break;
#endif
      int num_gen = exit_gen [ d ] . size ();
      std::cout << "ConleyIndexODE: Betti Number at Dimension " << d << ": " << num_gen << "\n";
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Compute the homology generators
    PRINT "ConleyIndexODE: Computing H_*(X) \n";
    Generators_t full_gen = MorseGenerators ( full_domain, cutoff_dimension );
    for ( int d = 0; d <= full_domain . dimension (); ++ d ) {
#ifdef CONLEYINDEXCUTOFF
      if ( d == CONLEYINDEXCUTOFF ) break;
#endif
      int num_gen = full_gen [ d ] . size ();
      std::cout << "ConleyIndexODE: Betti Number at Dimension " << d << ": " << num_gen << "\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////////

    double total_time = ((double)(clock()-total_time_start)/(double)CLOCKS_PER_SEC);
    PRINT "ConleyIndexODE: Time Elapsed (total) = " << total_time << "\n";

  }

/* Conley Index */
inline ConleyIndex_t::ConleyIndex_t ( void ) {
  undefined_ = false;
} /* ConleyIndex_t::ConleyIndex_t */

inline std::vector < SparseMatrix < Ring > > & ConleyIndex_t::data ( void ) {
  return data_;
}

inline const std::vector < SparseMatrix < Ring > > & ConleyIndex_t::data ( void ) const {
  return data_;
}

inline bool & ConleyIndex_t::undefined ( void ) {
  return undefined_;
} /* ConleyIndex_t::undefined */

inline const bool & ConleyIndex_t::undefined ( void ) const {
  return undefined_;
} /* ConleyIndex_t::undefined */

} // namespace chomp

#endif
