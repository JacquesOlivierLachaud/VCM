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

#include "DGtal/base/Common.h"
#include "DGtal/base/CountedPtr.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/CNormalVectorEstimator.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/VCMDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/TrueDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/ShapeGeometricFunctors.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/io/readers/MPolynomialReader.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

template <typename KSpace,
          typename ImplicitShape,
          typename Surface,
          typename Estimator>
void computeEstimation
( const po::variables_map& vm,     //< command-line parameters
  const KSpace& K,                 //< cellular grid space
  const ImplicitShape& shape, //< implicit shape "ground truth"
  const Surface& surface,     //< digital surface approximating shape
  Estimator& estimator )           //< an initialized estimator
{
}

template <typename KSpace,
          typename ImplicitShape,
          typename Surface,
          typename KernelFunction>
void chooseEstimator
( const po::variables_map& vm,     //< command-line parameters
  const KSpace& K,                 //< cellular grid space
  const ImplicitShape& shape, //< implicit shape "ground truth"
  const Surface& surface,     //< digital surface approximating shape
  const KernelFunction& chi ) //< the kernel function
{
  string nameEstimator = vm[ "estimator" ].as<string>();
  double h = vm["gridstep"].as<double>();
  if ( nameEstimator == "True" )
    {
      trace.beginBlock( "Setting up True estimator." );
      typedef ShapeGeometricFunctors::ShapeNormalVectorFunctor<ImplicitShape> NormalFunctor;
      typedef TrueDigitalSurfaceLocalEstimator<KSpace, ImplicitShape, NormalFunctor> TrueEstimator;
      TrueEstimator true_estimator;
      true_estimator.attach( shape );
      true_estimator.setParams( K, NormalFunctor(), 20, 0.1, 0.01 );
      true_estimator.init( h, surface.begin(), surface.end() );
      trace.endBlock();
    }
  else if ( nameEstimator == "VCM" )
    {
      trace.beginBlock( "Setting up VCM estimator." );
      typedef typename KSpace::Space Space;
      typedef typename Surface::DigitalSurfaceContainer SurfaceContainer;
      typedef ExactPredicateLpSeparableMetric<Space,2> Metric;
      typedef VoronoiCovarianceMeasureOnDigitalSurface<SurfaceContainer,Metric,
                                                       KernelFunction> VCMOnSurface;
      typedef VCMGeometricFunctors::VCMNormalVectorFunctor<VCMOnSurface> NormalFunctor;
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
      VCMNormalEstimator estimator;
      estimator.attach( surface );
      estimator.setParams( embType, R, r, chi, t, Metric(), true );
      estimator.init( h, surface.begin(), surface.end() );
      trace.endBlock();
    }

}

template <typename KSpace,
          typename ImplicitShape,
          typename Surface>
void chooseKernel
( const po::variables_map& vm,     //< command-line parameters
  const KSpace& K,                 //< cellular grid space
  const ImplicitShape& shape, //< implicit shape "ground truth"
  const Surface& surface )    //< digital surface approximating shape
{
  string kernel = vm[ "kernel" ].as<string>();
  double h = vm["gridstep"].as<double>();
  double r = vm["r-radius"].as<double>();
  double alpha = vm["alpha"].as<double>();
  if ( alpha != 0.0 ) r *= pow( h, alpha-1.0 );
  if ( kernel == "hat" ) {
    typedef typename KSpace::Point Point;
    typedef HatPointFunction<Point,double> KernelFunction;
    KernelFunction chi_r( 1.0, r );
    chooseEstimator( vm, K, shape, surface, chi_r );
  } else if ( kernel == "ball" ) {
    typedef typename KSpace::Point Point;
    typedef BallConstantPointFunction<Point,double> KernelFunction;
    KernelFunction chi_r( 1.0, r );
    chooseEstimator( vm, K, shape, surface, chi_r );
  }
}

          

///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("polynomial,p", po::value<string>(), "the implicit polynomial whose zero-level defines the shape of interest." )
    ("noise,N", po::value<double>()->default_value( 0.0 ), "the Kanungo noise level l=arg, with l^d the probability that a point at distance d is flipped inside/outside." )
    ("minAABB,a",  po::value<double>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
    ("maxAABB,A",  po::value<double>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
    ("gridstep,g", po::value< double >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " )
    ("estimator,e", po::value<string>()->default_value( "True" ), "the chosen normal estimator: True | VCM | Trivial" )
    ("R-radius,R", po::value<double>()->default_value( 5 ), "the constant for parameter R in R(h)=R h^alpha (VCM)." )
    ("r-radius,r", po::value<double>()->default_value( 3 ), "the constant for parameter r in r(h)=r h^alpha (VCM,Trivial)." )
    ("kernel,k", po::value<string>()->default_value( "hat" ), "the function chi_r, either hat or ball." )
    ("alpha", po::value<double>()->default_value( 0.0 ), "the parameter alpha in r(h)=r h^alpha (VCM)." )
    ("trivial-radius,t", po::value<double>()->default_value( 3 ), "the parameter t for reorienting the VCM." )
    ("embedding,E", po::value<int>()->default_value( 0 ), "the surfel -> point embedding: 0: Pointels, 1: InnerSpel, 2: OuterSpel." )
    ("export,x", po::value<string>(), "exports surfel normals which can be viewed with viewSetOfSurfels." )
    ("output,o", po::value<string>(), "the output basename." )
    ("noff,n", po::value<string>(), "exports the digital surface with normals as NOFF file <arg>" )
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
		<< "Computes a 3D normal vector field for several estimators, output it and compute some statistics." 
                << endl
		<< general_opt << "\n";
      cout << "Example:\n"
                << "./implicitShape3NormalEstimation -p \"90-3*x^2-2*y^2-z^2\" -o VCM-ellipse -a -10 -A 10 -e VCM -R 3 -r 3 -t 2 -E 0 -x vcm" << endl;
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
  typedef KSpace::Surfel Surfel;
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

  // Selecting a model of surface depending on noise / not noise.
  Scalar noiseLevel = vm[ "noise" ].as<double>();
  if ( noiseLevel == 0.0 )
    { // no noise
      trace.beginBlock( "Make digital surface..." );
      typedef LightImplicitDigitalSurface<KSpace,ImplicitDigitalShape> SurfaceContainer;
      typedef DigitalSurface< SurfaceContainer > Surface;
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
      chooseKernel( vm, K, *shape, *ptrSurface );
    }
  else
    { // noise
      trace.beginBlock( "Make digital surface..." );
      typedef KanungoNoise< ImplicitDigitalShape, Domain > KanungoPredicate;
      typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > SurfaceContainer;
      typedef DigitalSurface< SurfaceContainer > Surface;
      SurfelAdjacency< KSpace::dimension > surfAdj( true );
      Surfel bel;
      KanungoPredicate* noisified_dshape = new KanungoPredicate( *dshape, dshape->getDomain(), noiseLevel );
      // We have to search for a big connected component.
      CountedPtr<Surface> ptrSurface;
      double minsize = dshape->getUpperBound()[0] - dshape->getLowerBound()[0];
      unsigned int nb_surfels = 0;
      unsigned int tries = 0;
      do {
        try { // Search initial bel
          bel = Surfaces<KSpace>::findABel( K, *noisified_dshape, 10000 );
        } catch (DGtal::InputException e) {
          trace.error() << "ERROR Unable to find bel." << std::endl;
          return 3;
        }
        SurfaceContainer* surfaceContainer = new SurfaceContainer( K, *noisified_dshape, surfAdj, bel );
        ptrSurface = CountedPtr<Surface>( new Surface( surfaceContainer ) ); // acquired
        nb_surfels = ptrSurface->size();
      } while ( ( nb_surfels < 2 * minsize ) && ( tries++ < 150 ) );
      if( tries >= 150 )
        {
          trace.error() << "ERROR cannot find a proper bel in a big enough component." << std::endl;
          return 4;
        }
      trace.info() << "- surface component has " << nb_surfels << " surfels." << std::endl; 
      chooseKernel( vm, K, *shape, *ptrSurface );
      trace.endBlock();
    }
  return 0;
}

