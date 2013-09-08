#include <iostream>
#include <fstream>
// commandline arguments: inputfile threshold outputfile
#include "CImg.h"
using namespace cimg_library;

int main ( int argc, char * argv [] ) {
  if ( argc != 4 ) {
    std::cout << "command line args: inputfile thresholdvalue outfile\n";
    return 0;
  }
  CImg<double> image(argv[1]);
  CImgDisplay main_disp(image,"Image",0);
  CImg<double> thresimage(image.width(),image.height(),1,1,0);
  
  int threshold = atoi(argv[2]);
  std::ofstream outfile ( argv[3] );
  for (int i=0;i<image.width();i++){
    for (int j=0;j<image.height();j++) {
      if ( image(i,j,0,0) > threshold ) {
	outfile << "(" << i << ", " << j << ")\n";
	thresimage(i,j,0,0)=255;
      }
    }
  }
  outfile.close ();
  //for (int k=0;k<3;k++) darkimage(i,j,0,k)=image(i,j,0,k)/2;
  CImgDisplay dark_disp (thresimage,"Thresholded Image",0); 
  while (!main_disp.is_closed())
    main_disp.wait();
  return 0;
}
