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
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/CNormalVectorEstimator.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/VCMDigitalSurfaceNormalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/TrueDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/ShapeGeometricFunctors.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/io/readers/MPolynomialReader.h"

using namespace DGtal;

template <typename SCell, typename RealVector>
struct GradientMapAdapter {
  typedef std::map<SCell,RealVector> SCell2RealVectorMap;
  typedef SCell                                 Argument;
  typedef RealVector                               Value;
  GradientMapAdapter( ConstAlias<SCell2RealVectorMap> map )
    : myMap( map ) {}
  RealVector operator()( const Argument& arg ) const
  {
    typename SCell2RealVectorMap::const_iterator it = myMap->find( arg );
    if ( it != myMap->end() ) return it->second;
    else return RealVector();
  }
  CountedConstPtrOrConstPtr<SCell2RealVectorMap> myMap;
};

template <typename SCellEmbedder>
struct SCellEmbedderWithNormal : public SCellEmbedder
{
  using SCellEmbedder::space;
  using SCellEmbedder::operator();
  typedef typename SCellEmbedder::KSpace          KSpace;
  typedef typename SCellEmbedder::SCell            SCell;
  typedef typename SCellEmbedder::RealPoint    RealPoint;
  typedef SCell                                 Argument;
  typedef RealPoint                                Value;
  typedef typename KSpace::Space::RealVector  RealVector;
  typedef std::map<SCell,RealVector> SCell2RealVectorMap;
  typedef GradientMapAdapter<SCell,RealVector> GradientMap;

  SCellEmbedderWithNormal( ConstAlias<SCellEmbedder> embedder, 
                           ConstAlias<SCell2RealVectorMap> map )
    : SCellEmbedder( embedder ), myMap( map )
  {}
  
  GradientMap gradientMap() const
  {
    return GradientMap( myMap );
  }

  CountedConstPtrOrConstPtr<SCell2RealVectorMap> myMap;
};

template <typename DigitalSurface, 
          typename Estimator>
void exportNOFFSurface( const DigitalSurface& surface,
                        const Estimator& estimator,
                        std::ostream& output )
{
  typedef typename DigitalSurface::KSpace KSpace;
  typedef typename DigitalSurface::ConstIterator ConstIterator;
  typedef typename DigitalSurface::Surfel Surfel;
  typedef typename KSpace::SCell SCell;
  typedef typename Estimator::Quantity Quantity;
  const KSpace& ks = surface.container().space();
  std::map<Surfel,Quantity> normals;
  for ( ConstIterator it = surface.begin(), itE = surface.end(); it != itE; ++it )
    {
      Quantity n_est = estimator.eval( it );
      normals[ *it ] = n_est;
    }
  CanonicSCellEmbedder<KSpace> surfelEmbedder( ks );
  typedef SCellEmbedderWithNormal< CanonicSCellEmbedder<KSpace> > Embedder;
  Embedder embedder( surfelEmbedder, normals );
  surface.exportAs3DNOFF( output, embedder );
}


template <typename DigitalSurface, 
          typename Estimator,
          typename CompEstimator>
void computeNormalEstimation( const DigitalSurface& surface,
                              const Estimator& estimator,
                              const CompEstimator& compEstimator,
                              std::ostream& output_error,
                              std::ostream& output_export,
                              bool toExport )
{
  typedef typename DigitalSurface::KSpace KSpace;
  typedef typename DigitalSurface::ConstIterator ConstIterator;
  typedef typename DigitalSurface::Surfel Surfel;
  typedef typename Estimator::Quantity Quantity;
  typedef typename CompEstimator::Scalar Scalar;
  const KSpace& ks = surface.container().space();
  DGtal::Statistic<Scalar> stat_angle_error;
  for ( ConstIterator it = surface.begin(), itE = surface.end(); it != itE; ++it )
    {
      Quantity n_est = estimator.eval( it );
      Quantity n_comp_est = compEstimator.eval( it );
      Scalar angle_error = acos( n_est.dot( n_comp_est ) );
      stat_angle_error.addValue( angle_error );
      if ( toExport ) 
        {
          Surfel s = *it;
          output_export
            << "CellN " 
            << (256+ks.sKCoord( s, 0 )) << " " << (256+ks.sKCoord( s, 1 ))
            << " " << (256+ks.sKCoord( s, 2 )) << " " << ks.sSign( s ) << " 0.5 0.5 1.0" 
            << " " << n_est[ 0 ] << " " << n_est[ 1 ] << " " << n_est[ 2 ] << std::endl;
        }
    }
  stat_angle_error.terminate();
  output_error << "# Average error X of the angle between two vector estimations." << std::endl;
  output_error << "# h E[X] Var[X] Min[X] Max[X] Nb[X]." << std::endl;
  output_error << compEstimator.h() 
               << " " << stat_angle_error.mean()
               << " " << stat_angle_error.unbiasedVariance()
               << " " << stat_angle_error.min()
               << " " << stat_angle_error.max()
               << " " << stat_angle_error.samples()
               << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
  using namespace DGtal;

  // parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("polynomial,p", po::value<std::string>(), "the implicit polynomial whose zero-level defines the shape of interest." )
    ("output,o", po::value<std::string>(), "the output basename." )
    ("gridstep,g", po::value< double >()->default_value( 1.0 ), "the gridstep that the defines the digitization. " )
    ("estimator,e", po::value<std::string>()->default_value( "True" ), "the chosen normal estimator: True | VCM | Trivial" )
    ("minAABB,a",  po::value<double>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
    ("maxAABB,A",  po::value<double>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
    ("R-radius,R", po::value<double>()->default_value( 5 ), "the parameter R in the VCM." )
    ("r-radius,r", po::value<double>()->default_value( 3 ), "the parameter r in the VCM." )
    ("kernel,k", po::value<std::string>()->default_value( "hat" ), "the function chi_r, either hat or ball." )
    ("alpha", po::value<double>()->default_value( 0.25 ), "the parameter alpha in r(h)=r h^alpha in the VCM." )
    ("trivial-radius,t", po::value<double>()->default_value( 3 ), "the parameter r for the trivial normal estimator." )
    ("embedding,E", po::value<int>()->default_value( 0 ), "the surfel -> point embedding: 0: Pointels, 1: InnerSpel, 2: OuterSpel." )
    ("export,x", po::value<std::string>(), "exports surfel normals which can be viewed with viewSetOfSurfels." )
    ("noff,n", po::value<std::string>(), "exports the digital surface with normals as NOFF file <arg>" )
    ("noise,N", po::value<double>()->default_value( 0.0 ), "the Kanungo noise level l=arg, with l^d the probability that a point at distance d is flipped inside/outside." )
    ;  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify(vm);    
  if( !parseOK || vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " -p \"90-x^2-y^2-z^2\" -h 1 -R 5 -r 6 -t 2\n"
		<< "Computes a normal vector field, output it and compute some statistics." 
                << std::endl
		<< general_opt << "\n";
      std::cout << "Example:\n"
                << "./implicitShape3NormalEstimation -p \"90-3*x^2-2*y^2-z^2\" -o VCM/ellipse -a -10 -A 10 -e VCM -R 3 -r 3 -t 2 -E 0 -x vcm" << std::endl;
      return 0;
    }
  if ( ! vm.count( "polynomial" ) ) 
    {
      std::cerr << "Need parameter --polynomial" << std::endl;
      return 1;
    }
  trace.beginBlock( "Make shape..." );
  typedef Z3i::Space Space;
  typedef Z3i::KSpace KSpace;
  typedef double Scalar;
  typedef KSpace::Point Point;
  typedef KSpace::Surfel Surfel;
  typedef Space::RealPoint RealPoint;
  typedef MPolynomial< 3, Scalar > Polynomial3;
  typedef MPolynomialReader<3, Scalar> Polynomial3Reader;
  typedef ImplicitPolynomial3Shape<Space> ImplicitShape;
  typedef GaussDigitizer< Space, ImplicitShape > ImplicitDigitalShape;
  typedef ImplicitDigitalShape::Domain Domain;
  // typedef LightImplicitDigitalSurface<KSpace,ImplicitDigitalShape> SurfaceContainer;
  // typedef DigitalSurface< SurfaceContainer > Surface;
  typedef KanungoNoise< ImplicitDigitalShape, Domain > KanungoPredicate;
  typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > SurfaceContainer;
  typedef DigitalSurface< SurfaceContainer > Surface;

  std::string poly_str = vm[ "polynomial" ].as<std::string>();
  Polynomial3 poly;
  Polynomial3Reader reader;
  std::string::const_iterator iter = reader.read( poly, poly_str.begin(), poly_str.end() );
  if ( iter != poly_str.end() )
    {
      std::cerr << "ERROR: I read only <" << poly_str.substr( 0, iter - poly_str.begin() )
                << ">, and I built P=" << poly << std::endl;
      return 2;
    }
  CountedPtr<ImplicitShape> shape( new ImplicitShape( poly ) );
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

  // Extracts shape boundary
  Scalar noiseLevel = vm[ "noise" ].as<double>();
  KanungoPredicate* noisifiedObject = new KanungoPredicate( *dshape, dshape->getDomain(), 
                                                            noiseLevel );
  Surfel bel;
  try {
    bel = Surfaces<KSpace>::findABel( K, *noisifiedObject, 10000 );
  } catch (DGtal::InputException e) {
    std::cerr << "Unable to find bel." << std::endl;
    return 3;
  }
  SurfelAdjacency< KSpace::dimension > surfAdj( true );
  SurfaceContainer* surfaceContainer = 
    new SurfaceContainer( K, *noisifiedObject, surfAdj, bel );
  CountedPtr<Surface> ptrSurface( new Surface( surfaceContainer ) ); // acquired
  double minsize = dshape->getUpperBound()[0] - dshape->getLowerBound()[0];
  unsigned int tries = 0;
  while( ptrSurface->size() < 2 * minsize && ( tries++ < 150 ) )
    {
      delete surfaceContainer;
      bel = Surfaces< KSpace >::findABel( K, *noisifiedObject, 10000 );
      surfaceContainer = 
        new SurfaceContainer( K, *noisifiedObject, surfAdj, bel );
      ptrSurface = CountedPtr<Surface>( new Surface( surfaceContainer ) );
    }
  if( tries >= 150 )
    {
      std::cerr << "Can't found a proper bel. So .... I ... just ... kill myself." << std::endl;
      return false;
    }
  trace.info() << "- domain is " << domain << std::endl;
  trace.info() << "- surface has " << ptrSurface->size() << " surfels." << std::endl;
  trace.endBlock();

  std::string nameEstimator = vm[ "estimator" ].as<std::string>();
  std::string fname = vm[ "output" ].as<std::string>();
  std::ostringstream error_sstr;
  error_sstr << fname << "-" << nameEstimator << "-error-" << h << ".txt"; 
  std::ostringstream export_sstr;
  std::ostringstream noff_sstr;
  bool toExport = vm.count( "export" );
  std::string str = toExport ? vm[ "export" ].as<std::string>() : "tmp";
  export_sstr << fname << "-" << nameEstimator << "-" << str << "-" << h << ".txt"; 
  bool toNoff = vm.count( "noff" );
  str = toNoff ? vm[ "noff" ].as<std::string>() : "tmp";
  noff_sstr << fname << "-" << nameEstimator << "-" << str << "-" << h << ".off"; 
  std::ofstream error_output( error_sstr.str().c_str() );
  std::ofstream export_output( export_sstr.str().c_str() );
  std::ofstream noff_output( noff_sstr.str().c_str() );

  trace.beginBlock( "Setting up estimators." );
  typedef ShapeGeometricFunctors::ShapeNormalVectorFunctor<ImplicitShape> NormalFunctor;
  typedef TrueDigitalSurfaceLocalEstimator<KSpace, ImplicitShape, NormalFunctor> TrueEstimator;
  TrueEstimator true_estimator;
  true_estimator.setParams( K, NormalFunctor() );
  true_estimator.attach( shape );
  true_estimator.init( h, ptrSurface->begin(), ptrSurface->end() );

  if ( nameEstimator == "VCM" )
    {
      typedef ExactPredicateLpSeparableMetric<Space,2> Metric;
      int embedding = vm["embedding"].as<int>();
      Surfel2PointEmbedding embType = embedding == 0 ? Pointels :
                                      embedding == 1 ? InnerSpel : OuterSpel;     
      Scalar R = vm["R-radius"].as<double>();
      Scalar r = vm["r-radius"].as<double>();
      Scalar t = vm["trivial-radius"].as<double>();
      Scalar alpha = vm["alpha"].as<double>();
      R = R * pow( h, alpha-1.0 );
      r = r * pow( h, alpha-1.0 );
      std::string kernel = vm[ "kernel" ].as<std::string>();
      if ( kernel == "hat" ) {
        typedef HatPointFunction<Point,double> KernelFunction;
        typedef VoronoiCovarianceMeasureOnDigitalSurface<SurfaceContainer,Metric,
                                                         KernelFunction> VCMOnSurface;
        typedef VCMDigitalSurfaceNormalEstimator<SurfaceContainer,Metric,
                                                 KernelFunction> VCMNormalEstimator;
        trace.beginBlock("Computing VCM on surface." );
        KernelFunction chi_r( 1.0, r );
        CountedPtr<VCMOnSurface> vcm_surface( new VCMOnSurface( ptrSurface, embType,
                                                                R, r, chi_r,
                                                                t, Metric(), true ) );
        VCMNormalEstimator estimator( vcm_surface );
        trace.info() << "# VCM estimation: h=" << h << " R=" << R 
                     << " r=" << r << " hat" << std::endl;
        error_output << "# VCM estimation: h=" << h << " R=" << R 
                     << " r=" << r << " hat" << std::endl;
        if ( toExport )
          export_output << "# VCM estimation: h=" << h << " R=" << R 
                        << " r=" << r << " hat" << std::endl;
        trace.endBlock();
        trace.beginBlock("Statistics and export." );
        computeNormalEstimation( *ptrSurface, estimator, true_estimator,
                                 error_output, export_output, toExport );
        if ( vm.count( "noff" ) )
          exportNOFFSurface( *ptrSurface, estimator, noff_output );
        trace.endBlock();
      } else if ( kernel == "ball" ) {
        typedef BallConstantPointFunction<Point,double> KernelFunction;
        typedef VoronoiCovarianceMeasureOnDigitalSurface<SurfaceContainer,Metric,
                                                         KernelFunction> VCMOnSurface;
        typedef VCMDigitalSurfaceNormalEstimator<SurfaceContainer,Metric,
                                                 KernelFunction> VCMNormalEstimator;
        trace.beginBlock("Computing VCM on surface." );
        KernelFunction chi_r( 1.0, r );
        CountedPtr<VCMOnSurface> vcm_surface( new VCMOnSurface( ptrSurface, embType,
                                                                R, r, chi_r,
                                                                t, Metric(), true ) );
        VCMNormalEstimator estimator( vcm_surface );
        trace.info() << "# VCM estimation: h=" << h << " R=" << R 
                     << " r=" << r << " ball" << std::endl;
        error_output << "# VCM estimation: h=" << h << " R=" << R 
                     << " r=" << r << " ball" << std::endl;
        if ( toExport )
          export_output << "# VCM estimation: h=" << h << " R=" << R 
                        << " r=" << r << " ball" << std::endl;
        trace.endBlock();
        trace.beginBlock("Statistics and export." );
        computeNormalEstimation( *ptrSurface, estimator, true_estimator,
                                 error_output, export_output, toExport );
        if ( vm.count( "noff" ) )
          exportNOFFSurface( *ptrSurface, estimator, noff_output );
        trace.endBlock();
      }
    }
  error_output.close();
  export_output.close();
  noff_output.close();
  trace.endBlock();
  return 0;
}
