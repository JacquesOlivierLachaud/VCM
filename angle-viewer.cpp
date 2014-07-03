/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file noisyImplicitShape3NormalEstimation.cpp
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/05/02
 *
 * Estimates the normal vector field of an implicitly defined shape
 * for several estimators. The implicit shape is perturbated by a
 * Kanungo noise.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <QtGui/qapplication.h>

#include "DGtal/base/Common.h"
#include "DGtal/base/CountedPtr.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/TrueDigitalSurfaceLocalEstimator.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/ShapeGeometricFunctors.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/io/Display3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
  QApplication application(argc,argv);

  // parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("polynomial,p", po::value<string>(), "the implicit polynomial whose zero-level defines the shape of interest." )
    ("reach,R",  po::value<double>(), "the reach of the shape." )
    ("minAABB,a",  po::value<double>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
    ("maxAABB,A",  po::value<double>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
    ("gridstep,g", po::value< double >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " );

  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  if( !parseOK || vm.count("help"))
    {
      cout << "Usage: " << argv[0] << " -p \"90-x^2-y^2-z^2\" -h 1 -R 5 -r 6 -t 2\n"
		<< "Computes the angle of surface normal and displays potential problematic zones." 
                << endl
		<< general_opt << "\n";
      cout << "Example:\n"
           << "./angle-viewer -p \"90-3*x^2-2*y^2-z^2\"" << endl
           << " - ellipse  : 90-3*x^2-2*y^2-z^2 " << endl
           << " - torus    : -1*(x^2+y^2+z^2+6*6-2*2)^2+4*6*6*(x^2+y^2) " << endl
           << " - rcube    : 6561-x^4-y^4-z^4" << endl
           << " - goursat  : 8-0.03*x^4-0.03*y^4-0.03*z^4+2*x^2+2*y^2+2*z^2" << endl
           << " - distel   : 10000-(x^2+y^2+z^2+1000*(x^2+y^2)*(x^2+z^2)*(y^2+z^2))" << endl
           << " - leopold  : 100-(x^2*y^2*z^2+4*x^2+4*y^2+3*z^2)" << endl
           << " - diabolo  : x^2-(y^2+z^2)^2" << endl
           << " - heart    : -1*(x^2+2.25*y^2+z^2-1)^3+x^2*z^3+0.1125*y^2*z^3" << endl
           << " - crixxi   : -0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3" << endl
           << " - goursat_dodecahedron: z^6-5*(x^2+y^2)*z^4+5*(x^2+y^2)^2*z^2-2*(x^4-10*x^2*y^2+5*y^4)*x*z+1*(x^2+y^2+z^2)^3+(-1)*(5)^2*(x^2+y^2+z^2)^2+1*(5)^4*(x^2+y^2+z^2)+(-1)*(5)^6" << endl
           << " - goursat_icosahedron : z^6-5*(x^2+y^2)*z^4+5*(x^2+y^2)^2*z^2-2*(x^4-10*x^2*y^2+5*y^4)*x*z+(-1)*(x^2+y^2+z^2)^3+(0)*(5)^2*(x^2+y^2+z^2)^2+(-1)*(5)^4*(x^2+y^2+z^2)+(1)*(5)^6" << endl
           << " - goursat_60lines : z^6-5*(x^2+y^2)*z^4+5*(x^2+y^2)^2*z^2-2*(x^4-10*x^2*y^2+5*y^4)*x*z+(0)*(x^2+y^2+z^2)^3+(5)*(5)^2*(x^2+y^2+z^2)^2+(-45)*(5)^4*(x^2+y^2+z^2)+(71)*(5)^6" << endl;
      return 0;
    }
  if ( ! vm.count( "polynomial" ) ) 
    {
      cerr << "Need parameter --polynomial" << endl;
      return 1;
    }
  
  trace.beginBlock( "Make implicit shape..." );
  typedef Z3i::Space Space;
  typedef double Scalar;
  typedef MPolynomial< 3, Scalar > Polynomial3;
  typedef MPolynomialReader<3, Scalar> Polynomial3Reader;
  typedef ImplicitPolynomial3Shape<Space> ImplicitShape;
  string poly_str = vm[ "polynomial" ].as<string>();
  Polynomial3 poly;
  Polynomial3Reader reader;
  string::const_iterator iter = reader.read( poly, poly_str.begin(), poly_str.end() );
  if ( iter != poly_str.end() )
    {
      trace.error() << "ERROR reading polynomial: I read only <" << poly_str.substr( 0, iter - poly_str.begin() )
                    << ">, and I built P=" << poly << std::endl;
      return 2;
    }
  CountedPtr<ImplicitShape> shape( new ImplicitShape( poly ) ); // smart pointer
  trace.endBlock();

  trace.beginBlock( "Make implicit digital shape..." );
  typedef Z3i::KSpace KSpace;
  typedef KSpace::Point Point;
  typedef Space::RealPoint RealPoint;
  typedef GaussDigitizer< Space, ImplicitShape > ImplicitDigitalShape;
  typedef ImplicitDigitalShape::Domain Domain;
  Scalar min_x = vm[ "minAABB" ].as<double>();
  Scalar max_x = vm[ "maxAABB" ].as<double>();
  Scalar h = vm[ "gridstep" ].as<double>();
  RealPoint p1( min_x, min_x, min_x );
  RealPoint p2( max_x, max_x, max_x );
  CountedPtr<ImplicitDigitalShape> dshape( new ImplicitDigitalShape() );
  dshape->attach( *shape );
  dshape->init( p1, p2, h );
  Domain domain = dshape->getDomain();
  KSpace K;
  K.init( domain.lowerBound(), domain.upperBound(), true );
  trace.info() << "- domain is " << domain << std::endl;
  trace.endBlock();

  trace.beginBlock( "Make digital surface..." );
  typedef LightImplicitDigitalSurface<KSpace,ImplicitDigitalShape> SurfaceContainer;
  typedef DigitalSurface< SurfaceContainer > Surface;
  typedef Surface::ConstIterator ConstIterator;
  typedef typename Surface::Surfel Surfel;
  SurfelAdjacency< KSpace::dimension > surfAdj( true );
  Surfel bel;
  try {
    bel = Surfaces<KSpace>::findABel( K, *dshape, 10000 );
  } catch (DGtal::InputException e) {
    trace.error() << "ERROR Unable to find bel." << std::endl;
    return 3;
  }
  SurfaceContainer* surfaceContainer = new SurfaceContainer( K, *dshape, surfAdj, bel );
  CountedPtr<Surface> ptrSurface( new Surface( surfaceContainer ) ); // acquired
  trace.info() << "- surface component has " << ptrSurface->size() << " surfels." << std::endl; 
  trace.endBlock();

  trace.beginBlock( "Create estimator." );
  typedef ShapeGeometricFunctors::ShapeNormalVectorFunctor<ImplicitShape> NormalFunctor;
  typedef TrueDigitalSurfaceLocalEstimator<KSpace, ImplicitShape, NormalFunctor> NormalEstimator;
  typedef NormalEstimator::Quantity NQuantity;
  NormalEstimator normal_estimator;
  normal_estimator.attach( shape );
  normal_estimator.setParams( K, NormalFunctor(), 20, 0.1, 0.01 );
  normal_estimator.init( h, ptrSurface->begin(), ptrSurface->end() );

  typedef ShapeGeometricFunctors::ShapeMeanCurvatureFunctor<ImplicitShape> CurvatureFunctor;
  typedef TrueDigitalSurfaceLocalEstimator<KSpace, ImplicitShape, CurvatureFunctor> CurvatureEstimator;
  typedef CurvatureEstimator::Quantity CQuantity;
  CurvatureEstimator curv_estimator;
  curv_estimator.attach( shape );
  curv_estimator.setParams( K, CurvatureFunctor(), 20, 0.1, 0.01 );
  curv_estimator.init( h, ptrSurface->begin(), ptrSurface->end() );
  trace.endBlock();

  CQuantity mc = 0;
  for ( ConstIterator it = ptrSurface->begin(), itE= ptrSurface->end(); it != itE; ++it )
    {
      CQuantity c = curv_estimator.eval( it );
      mc = std::max( abs( c ), mc ); 
    }
  double R = 1.0 / mc;
  if ( vm.count( "reach" ) ) R = vm[ "reach" ].as<double>();
  typedef Viewer3D<Space,KSpace> MyViewer3D;
  typedef Display3DFactory<Space,KSpace> MyDisplay3DFactory;
  MyViewer3D viewer( K );
  viewer.show(); 
  viewer << SetMode3D( bel.className(), "Basic" );
  double s_manifold = h / (0.794*R);
  double s_homeo_xi = h*2.0*sqrt(3)/R;
  trace.info() << "Estimated reach    = " << (1.0/mc) << std::endl;
  trace.info() << "Chosen reach       = " << R << std::endl;
  trace.info() << "Threshold manifold = " << s_manifold << std::endl;
  trace.info() << "Threshold homeo xi = " << s_homeo_xi << std::endl;

  double areas[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
  for ( ConstIterator it = ptrSurface->begin(), itE= ptrSurface->end(); it != itE; ++it )
    {
      viewer.setFillColor( Color::White );
      Dimension k    = K.sOrthDir( *it );
      NQuantity n    = normal_estimator.eval( it );
      double area    = abs( n[ k ] );
      double angle_x = acos( abs( n[0] ) );
      double angle_y = acos( abs( n[1] ) );
      double angle_z = acos( abs( n[2] ) );
      bool non_manifold = ( ( angle_x < s_manifold ) 
                            || ( angle_y < s_manifold ) 
                            || ( angle_z < s_manifold ) );
      bool non_homeo_xi = ! ( ( abs( n[ 0 ] ) > s_homeo_xi )
                              && ( abs( n[ 1 ] ) > s_homeo_xi )
                              && ( abs( n[ 2 ] ) > s_homeo_xi ) );
      int color = ( non_manifold ? 2 : 0 ) + ( non_homeo_xi ? 1 : 0 );
      if ( color == 0 )      viewer.setFillColor( Color::White );
      else if ( color == 1 ) viewer.setFillColor( Color( 160, 160, 180 ) );
      else if ( color == 2 ) viewer.setFillColor( Color::Yellow );
      else if ( color == 3 ) viewer.setFillColor( Color( 60, 60, 80 ) );
      areas[ 0 ] += area;
      areas[ 1 ] += color == 1 ? area : 0.0;
      areas[ 2 ] += color == 2 ? area : 0.0;
      areas[ 3 ] += color == 3 ? area : 0.0;
      // for ( int i = 0; i <= color; ++i ) areas[ i ] += area;
      MyDisplay3DFactory::drawOrientedSurfelWithNormal( viewer, *it, n, false );
    }
  for ( int i = 0; i < 4; ++i ) areas[ i ] *= h*h;
  std::cout << setprecision( 12 );
  trace.info() << "Area              = " << setprecision( 12 ) << areas[ 0 ] << " " << "100.0 %" << std::endl;
  trace.info() << "Area Non-Homeo-Xi = " << setprecision( 12 ) << areas[ 1 ] << " " << (100.0*areas[1]/areas[0]) << " %" << std::endl;
  trace.info() << "Area Non-Manifold = " << setprecision( 12 ) << areas[ 2 ] << " " << (100.0*areas[2]/areas[0]) << " %" << std::endl;
  trace.info() << "Area NM + NHXi    = " << setprecision( 12 ) << areas[ 3 ] << " " << (100.0*areas[3]/areas[0]) << " %" << std::endl;
  viewer << MyViewer3D::updateDisplay;
  application.exec();
 
  return 0;
}
