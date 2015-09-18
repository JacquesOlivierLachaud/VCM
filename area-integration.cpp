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
#include "DGtal/io/viewers/Viewer3D.h"

#include "DGtal/base/Common.h"
#include "DGtal/base/CountedPtr.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/CNormalVectorEstimator.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/VCMDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/TrueDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"
#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/ShapeGeometricFunctors.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

/**
   Computes the normal estimations. Outputs statistics or export cell geometry.
 */
template <typename KSpace,
          typename ImplicitShape,
          typename Surface,
          typename Estimator>
void computeEstimation
( const po::variables_map& vm,     //< command-line parameters
  const KSpace& K,                 //< cellular grid space
  const ImplicitShape& shape,      //< implicit shape "ground truth"
  const Surface& surface,          //< digital surface approximating shape
  Estimator& estimator )           //< an initialized estimator
{
  typedef typename Surface::ConstIterator ConstIterator;
  typedef typename Surface::Surfel Surfel;
  typedef typename Estimator::Quantity Quantity;
  typedef double Scalar;
  typedef DepthFirstVisitor< Surface > Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef typename VisitorRange::ConstIterator VisitorConstIterator;
  
  std::string fname = vm[ "output" ].as<std::string>();
  string nameEstimator = vm[ "estimator" ].as<string>();

  trace.beginBlock( "Computing " + nameEstimator + " estimations." );
  CountedPtr<VisitorRange> range( new VisitorRange( new Visitor( surface, *(surface.begin()) )) );
  std::vector<Quantity> n_estimations;
  estimator.eval( range->begin(), range->end(), std::back_inserter( n_estimations ) );
  trace.info() << "- nb estimations  = " << n_estimations.size() << std::endl;
  trace.endBlock();

  trace.beginBlock( "Computing areas." );
  range = CountedPtr<VisitorRange>( new VisitorRange( new Visitor( surface, *(surface.begin()) )) );
  double area_est   = 0.0; // normal integration with absolute value.
  unsigned int i = 0;
  for ( typename VisitorRange::ConstIterator it = range->begin(), itE = range->end(); 
        it != itE; ++it, ++i )
    {
      Surfel s = *it;
      Dimension k = K.sOrthDir( s );
      area_est  += abs( n_estimations[ i ][ k ] );
    }  
  double h = vm["gridstep"].as<double>();
  trace.info() << setprecision(10) << "- Area_est     " << ( area_est * h * h )   << std::endl;
  std::ostringstream area_sstr;
  area_sstr << fname << "-" << nameEstimator << "-area-" << h << ".txt"; 
  std::ofstream area_output( area_sstr.str().c_str() );
  area_output << "# Area estimation by digital surface integration." << std::endl;
  area_output << "# X: " << nameEstimator << std::endl;
  area_output << "# h Area[X] nb_surf" << std::endl;
  area_output << setprecision(10) << h
              << " " << ( area_est * h * h )
              << " " << i << std::endl;
  area_output.close();
  trace.endBlock();
}

template <typename KSpace,
          typename ImplicitShape,
          typename Surface,
          typename KernelFunction,
          typename PointPredicate>
void chooseEstimator
( const po::variables_map& vm,     //< command-line parameters
  const KSpace& K,                 //< cellular grid space
  const ImplicitShape& shape, //< implicit shape "ground truth"
  const Surface& surface,     //< digital surface approximating shape
  const KernelFunction& chi,  //< the kernel function
  const PointPredicate& ptPred )   //< analysed implicit digital shape as a PointPredicate
{
  using namespace DGtal::functors;
  string nameEstimator = vm[ "estimator" ].as<string>();
  double h = vm["gridstep"].as<double>();
  typedef ShapeGeometricFunctors::ShapeNormalVectorFunctor<ImplicitShape> NormalFunctor;
  typedef TrueDigitalSurfaceLocalEstimator<KSpace, ImplicitShape, NormalFunctor> TrueEstimator;
  TrueEstimator true_estimator;
  true_estimator.attach( shape );
  true_estimator.setParams( K, NormalFunctor(), 20, 0.1, 0.01 );
  true_estimator.init( h, surface.begin(), surface.end() );
  if ( nameEstimator == "True" )
    {
      trace.beginBlock( "Chosen estimator is: True." );
      typedef TrueDigitalSurfaceLocalEstimator<KSpace, ImplicitShape, NormalFunctor> Estimator;
      int maxIter     = vm["maxiter"].as<int>();
      double accuracy = vm["accuracy"].as<double>();
      double gamma    = vm["gamma"].as<double>();
      Estimator estimator;
      estimator.attach( shape );
      estimator.setParams( K, NormalFunctor(), maxIter, accuracy, gamma );
      estimator.init( h, surface.begin(), surface.end() );
      trace.endBlock();
      computeEstimation( vm, K, shape, surface, estimator );
    }
  else if ( nameEstimator == "VCM" )
    {
      trace.beginBlock( "Chosen estimator is: VCM." );
      typedef typename KSpace::Space Space;
      typedef typename Surface::DigitalSurfaceContainer SurfaceContainer;
      typedef ExactPredicateLpSeparableMetric<Space,2> Metric;
      typedef VoronoiCovarianceMeasureOnDigitalSurface<SurfaceContainer,Metric,
                                                       KernelFunction> VCMOnSurface;
      typedef VCMNormalVectorFunctor<VCMOnSurface> NormalFunctor;
      typedef VCMDigitalSurfaceLocalEstimator<SurfaceContainer,Metric,
                                              KernelFunction, NormalFunctor> VCMNormalEstimator;
      int embedding = vm["embedding"].as<int>();
      Surfel2PointEmbedding embType = embedding == 0 ? Pointels :
                                      embedding == 1 ? InnerSpel : OuterSpel;     
      double R = vm["R-radius"].as<double>();
      double r = vm["r-radius"].as<double>();
      double t = vm["trivial-radius"].as<double>();
      double alpha = vm["alpha"].as<double>();
      if ( alpha != 0.0 ) R *= pow( h, alpha-1.0 );
      if ( alpha != 0.0 ) r *= pow( h, alpha-1.0 );
      trace.info() << "- R=" << R << " r=" << r << " t=" << t << std::endl;
      VCMNormalEstimator estimator;
      estimator.attach( surface );
      estimator.setParams( embType, R, r, chi, t, Metric(), true );
      estimator.init( h, surface.begin(), surface.end() );
      trace.endBlock();
      computeEstimation( vm, K, shape, surface, estimator );
    }
  else if ( nameEstimator == "II" )
    {
      trace.beginBlock( "Chosen estimator is: II." );
      typedef typename KSpace::Space Space;
      typedef HyperRectDomain<Space> Domain;
      typedef ImageContainerBySTLVector< Domain, bool> Image;
      typedef typename Domain::ConstIterator DomainConstIterator;
      typedef SimpleThresholdForegroundPredicate<Image> ThresholdedImage;
      typedef IINormalDirectionFunctor<Space> IINormalFunctor;
      typedef IntegralInvariantCovarianceEstimator<KSpace, ThresholdedImage, IINormalFunctor> IINormalEstimator;
      double r = vm["r-radius"].as<double>();
      double alpha = vm["alpha"].as<double>();
      if ( alpha != 0.0 ) r *= pow( h, alpha-1.0 );
      trace.info() << " r=" << r << std::endl;
      trace.beginBlock( "Preparing characteristic set." );
      Domain domain( K.lowerBound(), K.upperBound() );
      Image image( domain );
      for ( DomainConstIterator it = domain.begin(), itE = domain.end(); it != itE; ++it )
        {
          image.setValue( *it, ptPred( *it ) );
        }
      trace.endBlock();
      trace.beginBlock( "Initialize II estimator." );
      ThresholdedImage thresholdedImage( image, false );
      IINormalEstimator ii_estimator( K, thresholdedImage );
      ii_estimator.setParams( r );
      ii_estimator.init( h, surface.begin(), surface.end() );
      trace.endBlock();
      trace.endBlock();
      computeEstimation( vm, K, shape, surface, ii_estimator );
   }
  else if ( nameEstimator == "Trivial" )
    {
      trace.beginBlock( "Chosen estimator is: Trivial." );
      typedef HatFunction<double> Functor;
      typedef typename KSpace::Space Space;
      typedef typename KSpace::Surfel Surfel;
      typedef typename Surface::DigitalSurfaceContainer SurfaceContainer;
      typedef ExactPredicateLpSeparableMetric<Space,2> Metric;
      typedef ElementaryConvolutionNormalVectorEstimator< Surfel, CanonicSCellEmbedder<KSpace> > 
        SurfelFunctor;
      typedef LocalEstimatorFromSurfelFunctorAdapter< SurfaceContainer, Metric, SurfelFunctor, Functor>
        NormalEstimator;

      double t = vm["trivial-radius"].as<double>();
      Functor fct( 1.0, t );
      CanonicSCellEmbedder<KSpace> canonic_embedder( K );
      SurfelFunctor surfelFct( canonic_embedder, 1.0 );
      NormalEstimator estimator;
      estimator.attach( surface );
      estimator.setParams( Metric(), surfelFct, fct, t );
      estimator.init( 1.0, surface.begin(), surface.end() );
      trace.endBlock();
      computeEstimation( vm, K, shape, surface, estimator );
    }
}


///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
  QApplication application(argc,argv);

  // parse command line ----------------------------------------------
  using namespace DGtal::functors;
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("polynomial,p", po::value<string>(), "the implicit polynomial whose zero-level defines the shape of interest." )
    ("minAABB,a",  po::value<double>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
    ("maxAABB,A",  po::value<double>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
    ("gridstep,g", po::value< double >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " )
    ("estimator,e", po::value<string>()->default_value( "True" ), "the chosen normal estimator: True | VCM | II | Trivial" )
    ("R-radius,R", po::value<double>()->default_value( 5 ), "the constant for parameter R in R(h)=R h^alpha (VCM)." )
    ("r-radius,r", po::value<double>()->default_value( 3 ), "the constant for parameter r in r(h)=r h^alpha (VCM,II,Trivial)." )
    ("kernel,k", po::value<string>()->default_value( "hat" ), "the function chi_r, either hat or ball." )
    ("alpha", po::value<double>()->default_value( 0.0 ), "the parameter alpha in r(h)=r h^alpha (VCM)." )
    ("trivial-radius,t", po::value<double>()->default_value( 3 ), "the parameter t defining the radius for the Trivial estimator. Also used for reorienting the VCM." )
    ("embedding,E", po::value<int>()->default_value( 0 ), "the surfel -> point embedding for VCM estimator: 0: Pointels, 1: InnerSpel, 2: OuterSpel." )
    ("maxiter", po::value<int>()->default_value( 20 ), "the maximal number of iterations for True estimator (default is 20).")
    ("accuracy", po::value<double>()->default_value( 0.1 ), "the maximal accuracy for True estimator (default is 0.1).")
    ("gamma", po::value<double>()->default_value( 0.01 ), "the maximal gamma step for True estimator (default is 0.01).")
    ("output,o", po::value<string>()->default_value( "output" ), "the output basename. All generated files will have the form <arg>-*, for instance <arg>-area-<gridstep>.txt." )
    ;

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
		<< "Computes the area of a digital surface, defined as an implicit polynomial surface, by integration of normal estimation." 
                << endl
		<< general_opt << "\n";
      cout << "Example:\n"
           << "./area-integration -p \"81-x^2-y^2-z^2\" -e VCM -R 3 -r 3 -g 0.5  # aire de la sphere de rayon 9, discrétisé au pas 0.5" << endl
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

  string kernel = vm[ "kernel" ].as<string>();
  double r = vm["r-radius"].as<double>();
  double alpha = vm["alpha"].as<double>();
  if ( alpha != 0.0 ) r *= pow( h, alpha-1.0 );
  if ( kernel == "hat" ) {
    typedef typename KSpace::Point Point;
    typedef HatPointFunction<Point,double> KernelFunction;
    KernelFunction chi_r( 1.0, r );
    trace.info() << "- kernel hat r = " << r << std::endl; 
    chooseEstimator( vm, K, *shape, *ptrSurface, chi_r, *dshape );
  } else if ( kernel == "ball" ) {
    typedef typename KSpace::Point Point;
    typedef BallConstantPointFunction<Point,double> KernelFunction;
    KernelFunction chi_r( 1.0, r );
    trace.info() << "- kernel ball r = " << r << std::endl; 
    chooseEstimator( vm, K, *shape, *ptrSurface, chi_r, *dshape );
  }
  return 0;
}
