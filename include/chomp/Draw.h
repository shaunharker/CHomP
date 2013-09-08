/*
 *  Draw.h 
 *  Shaun Harker
 *  9/30/11
 *
 */

#ifndef CHOMP_DRAW_H
#define CHOMP_DRAW_H 
#include <iostream>
#include <vector>
/* include the X library headers */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
// X11 has some foolish macros that conflict with us
#undef index
#undef Complex

#ifdef None
#undef None
#endif

#include "chomp/Chain.h"

namespace chomp {
  
class GraphicsWindow {
public:
  GraphicsWindow( const char * title );
  ~GraphicsWindow ( void );

  void clear ( void );
  char wait ( void ); // return a keystroke
  void line ( int color, int x0, int y0, int x1, int y1 );
  void rect (int color, int x, int y, int w, int h);
  std::vector < int > palette;

protected:
  Display *dis;
  int screen;
  Window win;
  GC gc;
};
  
inline GraphicsWindow::GraphicsWindow( const char * title ) {  
  unsigned long black,white;
  dis=XOpenDisplay(NULL);
  screen=DefaultScreen(dis);
  black=BlackPixel(dis,screen),
  white=WhitePixel(dis, screen);
  win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0, 
                          512, 512, 5,black, white);
  XSetStandardProperties(dis,win, title,"Hi", 0L/* None */,NULL,0,NULL);
  XSelectInput(dis, win, KeyPressMask);
  // get Graphics Context
  gc=XCreateGC(dis, win, 0,0);        
  XSetBackground(dis,gc,white);
  XSetForeground(dis,gc,black);
  XClearWindow(dis, win);
  XMapRaised(dis, win);
  // SRH BEGIN
  /* first, find the default visual for our screen. */
  Visual* default_visual = DefaultVisual(dis, DefaultScreen(dis));
  /* this creates a new color map. the number of color entries in this map */
  /* is determined by the number of colors supported on the given screen.  */
  Colormap my_colormap = XCreateColormap(dis,
                                         win,
                                         default_visual,
                                         AllocNone);
  palette . resize ( 256 );
  for ( int i = 0; i <  128 ; ++ i ) {
    XColor rgb_color;
    
    rgb_color.red = 0; 
    rgb_color.green = 0;
    rgb_color.blue = i * 512;
    //rgb_color.red = rand () % 65536; 
    //rgb_color.green = rand () % 65536; 
    //rgb_color.blue = rand () % 65536;
    /*Status rc = */XAllocColor(dis,
                                my_colormap,
                                &rgb_color);
    palette [ i ] = rgb_color . pixel;
  }
  
  for ( int i = 0; i < 128; ++ i ) {
    XColor rgb_color;
    rgb_color.red = i * 512; 
    rgb_color.green = 0;
    rgb_color.blue = 0;
    XAllocColor(dis, my_colormap, &rgb_color);
    palette [ i + 128 ] = rgb_color . pixel;
  } 
  
}

inline GraphicsWindow::~GraphicsWindow ( void ) {
  XFreeGC(dis, gc);
  XDestroyWindow(dis,win);
  XCloseDisplay(dis);     
}

inline void GraphicsWindow::clear ( void ) {
  XClearWindow(dis, win);
}

inline char GraphicsWindow::wait ( void ) {
  XFlush ( dis );
  XEvent event;
  while(1){
		XNextEvent(dis, &event);
    //std::cout << "event occurred.\n";
		switch(event.type){
      case KeyPress:
        
				char string[25];
				KeySym keysym;
				XLookupString(&event.xkey, string, 25, &keysym, NULL);
				return string [ 0 ];
        break;
        
		}
	}
}

inline void GraphicsWindow::line (int color, int x0, int y0, int x1, int y1) {
  XSetForeground(dis, gc, palette[color]);
  XDrawLine(dis, win, gc, x0, y0, x1, y1);
}

inline void GraphicsWindow::rect (int color, int x, int y, int w, int h) {
  XSetForeground(dis, gc, palette[color]);
  //std::cout << x << " " << y << " " << w << " " << h << "\n";
  XFillRectangle(dis, win, gc, x, y, w + 3, h + 3);
}

/// ComplexVisualization. A class for visualizing 2D FiberComplexes.

class ComplexVisualization : public GraphicsWindow {
public:
  typedef std::pair < Rect, int > Object;
  ComplexVisualization ( const char * title );
  void drawObject ( const Object & obj );
  void drawNow ( void );
  void drawRect ( const Rect & p, int color );
  template < class C >
  void drawCell ( const C & fiber, Index i, int d, int color );
  template < class C >
  void drawSpecialCell ( const C & fiber, Index i, int d, int color );
  template < class C >
  void drawChain ( const C & fiber, const Chain &, int color );
  template < class C >
  void drawComplex ( const C & fiber, int color );
  template < class P >
  inline void drawRelativeComplex ( const P & pair, int color1, int color2 );
  void explore ( void );
private:
  bool rescale_cell ( Rect & bounds );
  std::vector < std::pair < Rect, int > > cells_;
  std::vector < std::pair < Rect, int > > special_cells_;

  Rect bounds_;
};


/* ComplexVisualization */

inline ComplexVisualization
::ComplexVisualization ( const char * title ) : GraphicsWindow ( title ), bounds_ ( 2 ) {
  bounds_ . lower_bounds [ 0 ] = 0.0f;
  bounds_ . lower_bounds [ 1 ] = 0.0f;
  bounds_ . upper_bounds [ 0 ] = 1.0f;
  bounds_ . upper_bounds [ 1 ] = 1.0f;
}


inline void ComplexVisualization::drawObject ( const Object & object ) {
  float x = object . first . lower_bounds [ 0 ];
  float y = object . first . lower_bounds [ 1 ];
  float width = object . first . upper_bounds [ 0 ] - object . first . lower_bounds [ 0 ];
  float height = object . first . upper_bounds [ 1 ] - object . first . lower_bounds [ 1 ];
  
  float dx = 480.0f * (x - bounds_ . lower_bounds [ 0 ]) / ( bounds_ . upper_bounds [ 0 ] 
                                                  - bounds_ . lower_bounds [ 0 ] );
  float dy = 480.0f * (y - bounds_ . lower_bounds [ 1 ] ) / ( bounds_ . upper_bounds [ 1 ] 
                                                            - bounds_ . lower_bounds [ 1 ] );
  float dw = width * 480.0f / ( bounds_ . upper_bounds [ 0 ] 
                              - bounds_ . lower_bounds [ 0 ] );
  float dh = height * 480.0f / ( bounds_ . upper_bounds [ 1 ] 
                               - bounds_ . lower_bounds [ 1 ] );
  
  if ( dx + dw < 0 ) return;
  if ( dy + dh < 0 ) return;
  if ( dx > 480 ) return;
  if ( dy > 480 ) return;
  
  if ( object . second > 100 ) {
  GraphicsWindow::rect ( object . second, 
                        (int) (dx + dw / 4.0f),  
                        (int) (dy + dh / 4.0f), 
                        (int) (dw * 2.0f / 4.0f), 
                        (int) (dh * 2.0f / 4.0f));
  } else {
    GraphicsWindow::rect ( object . second, 
                          (int) (dx ),  
                          (int) (dy ), 
                          (int) (dw ), 
                          (int) (dh ) );  
  }
  //std::cout << object.first << "\n";
  //std::cout << bounds_ << "\n";
  //std::cout << "(dx,dy,dw,dh) " << dx << " , " << dy << ", " << (int)dw <<","<< (int) dh<<"\n";
}

inline void ComplexVisualization::drawNow ( void ) {
  clear ();
  BOOST_FOREACH ( const Object & object, cells_ ) drawObject ( object );
  BOOST_FOREACH ( const Object & object, special_cells_ ) drawObject ( object );
  XFlush ( dis );
}


inline void ComplexVisualization::drawRect ( const Rect & p, int color ) {
  special_cells_  . push_back ( std::make_pair ( p, color ) );
}

template < class C >
inline void ComplexVisualization::drawCell ( const C & fiber, Index i, int d, int color ) {
  //std::cout << "drawCell (" << i << ", " << d << ")\n";
  Rect bounds = fiber . geometry ( i, d );
  cells_ . push_back ( std::make_pair ( bounds, color ) );
}


template < class C >
inline void ComplexVisualization::drawSpecialCell ( const C & fiber, Index i, int d, int color ) {
  Rect bounds = fiber . geometry ( i, d );
  special_cells_ . push_back ( std::make_pair ( bounds, color ) );
}


template < class C >
inline void ComplexVisualization::drawChain ( const C & fiber, const Chain & chain, int color ) {
  BOOST_FOREACH ( const Term & term, chain () ) {
    drawSpecialCell ( fiber, term . index (), chain . dimension (), color );
  }
}


template < class C >
inline void ComplexVisualization::drawComplex ( const C & fiber, int color ) {
  // TODO
  for ( int d = 0; d <= fiber . dimension (); ++ d ) {
    for ( Index i = 0; i < fiber . size ( d ); ++ i ) {
      drawCell ( fiber, i, d, color );
    }
  }

}

template < class P >
inline void ComplexVisualization::drawRelativeComplex ( const P & pair, int color1, int color2 ) {
  for ( int d = 0; d <= pair . pair () . dimension (); ++ d ) {
    for ( Index i = 0; i < pair . pair () . size ( d ); ++ i ) {
      Index j = pair . pair () . indexToCell ( i, d );
      drawCell ( pair . base (), j, d, color1 );
    }
  }
  for ( int d = 0; d <= pair . relative () . dimension (); ++ d ) {
    for ( Index i = 0; i < pair . relative () . size ( d ); ++ i ) {
      Index j = pair . relative () . indexToCell ( i, d );
      drawCell ( pair . base (), j, d, color2 );
    }
  }
  
}

/*
inline void ComplexVisualization::drawRelativeFiberComplex ( const FiberComplex & FiberComplex, 
                                                 const SubFiberComplex & subFiberComplex, 
                                                 int color, int color2 ) {
  for ( int d = 0; d <= FiberComplex . dimension (); ++ d ) {
    for ( Index i = 0; i <= FiberComplex . size ( d ); ++ i ) {
      drawCell ( FiberComplex, i, d, color );
    }
  }
  for ( int d = 0; d <= subFiberComplex . dimension (); ++ d ) {
    for ( Index i = 0; i <= subFiberComplex . size ( d ); ++ i ) {
      drawCell ( FiberComplex, i, d, color2 );
    }
  }

}
*/

inline void ComplexVisualization::explore ( void ) {

  if ( not special_cells_ . empty () ) {
    bounds_ . lower_bounds [ 0 ] = 1000.0f;
    bounds_ . lower_bounds [ 1 ] = 1000.0f;
    bounds_ . upper_bounds [ 0 ] = 0.0f;
    bounds_ . upper_bounds [ 1 ] = 0.0f;
    BOOST_FOREACH ( const Object & object, cells_ ) {
      Rect bounds = object . first;
      bounds_ . lower_bounds [ 0 ] = std::min ( bounds_ . lower_bounds [ 0 ], 
                                                bounds . lower_bounds [ 0 ]);
      bounds_ . lower_bounds [ 1 ] = std::min ( bounds_ . lower_bounds [ 1 ], 
                                               bounds . lower_bounds [ 1 ]);      
      bounds_ . upper_bounds [ 0 ] = std::max ( bounds_ . upper_bounds [ 0 ], 
                                                bounds . upper_bounds [ 0 ]);      
      bounds_ . upper_bounds [ 1 ] = std::max ( bounds_ . upper_bounds [ 1 ],
                                                bounds . upper_bounds [ 1 ]);
    }

  }
  
  drawNow ();
  while ( 1 ) {
    char c = wait ();
    float xscale = bounds_ . upper_bounds[ 0 ] - bounds_ . lower_bounds [ 0 ];
    float yscale = bounds_ . upper_bounds[ 1 ] - bounds_ . lower_bounds [ 1 ];
    
    if ( c == 'q' ) return;
    if ( c == 'z' ) drawNow ();

    if ( c == '-' ) {
      bounds_ . lower_bounds [ 0 ] -= xscale / 8.0f;
      bounds_ . lower_bounds [ 1 ] -= yscale / 8.0f;
      bounds_ . upper_bounds [ 0 ] += xscale / 8.0f;
      bounds_ . upper_bounds [ 1 ] += yscale / 8.0f;
      drawNow ();
    }
    if ( c == '=' ) {

      bounds_ . lower_bounds [ 0 ] += xscale / 8.0f;
      bounds_ . lower_bounds [ 1 ] += yscale / 8.0f;
      bounds_ . upper_bounds [ 0 ] -= xscale / 8.0f;
      bounds_ . upper_bounds [ 1 ] -= yscale / 8.0f;
      drawNow ();
    }
    if ( c == 'w' ) {
      bounds_ . lower_bounds [ 1 ] += yscale / 8.0f;
      bounds_ . upper_bounds [ 1 ] += yscale / 8.0f;
      drawNow ();
    }
    if ( c == 'a' ) {
      bounds_ . lower_bounds [ 0 ] -= xscale / 8.0f;
      bounds_ . upper_bounds [ 0 ] -= xscale / 8.0f;
      drawNow ();
    }
    if ( c == 's' ) {
      bounds_ . lower_bounds [ 1 ] -= yscale / 8.0f;
      bounds_ . upper_bounds [ 1 ] -= yscale / 8.0f;
      drawNow ();
    }   
    if ( c == 'd' ) {
      bounds_ . lower_bounds [ 0 ] += xscale / 8.0f;
      bounds_ . upper_bounds [ 0 ] += xscale / 8.0f;
      drawNow ();
    }

  }
}
  
} // namespace chomp

#endif
