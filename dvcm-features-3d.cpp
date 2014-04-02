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
 * @file dvcm-3d.cpp
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/05/02
 *
 * Computes the 3d voronoi map of a list of digital points.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <QtGui/qapplication.h>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/SimpleMatrix.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/math/EigenValues3D.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/IntervalForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/geometry/volumes/estimation/VoronoiCovarianceMeasure.h"
#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/ElementaryConvolutionNormalVectorEstimator.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/colormaps/GradientColorMap.h"


///////////////////////////////////////////////////////////////////////////////

using namespace std;

typedef DGtal::Z3i::Space Space;
typedef DGtal::Z3i::KSpace KSpace;
typedef DGtal::Z3i::Vector Vector;
typedef DGtal::Z3i::Point Point;
typedef DGtal::Z3i::RealPoint RealPoint;
typedef DGtal::Z3i::RealVector RealVector;
typedef DGtal::HyperRectDomain<Space> Domain;
typedef DGtal::ImageContainerBySTLVector<Domain,bool> CharacteristicSet;
typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> Metric; // L2-metric
typedef DGtal::EigenValues3D<double> LinearAlgebraTool;
typedef LinearAlgebraTool::Matrix33 Matrix33;
typedef LinearAlgebraTool::Vector3 Vector3;
typedef KSpace::Surfel Surfel;
typedef KSpace::SCell SCell;


// Model of CPointPredicate
struct CharacteristicSetPredicate {
  typedef CharacteristicSetPredicate Self;
  typedef DGtal::Z3i::Point Point;
  CharacteristicSetPredicate() : ptrSet( 0 ) {}
  CharacteristicSetPredicate( const CharacteristicSet& aSet) : ptrSet( &aSet ) {}
  CharacteristicSetPredicate( const Self& other ) : ptrSet( other.ptrSet ) {}
  Self& operator=( const Self& other ) { ptrSet = other.ptrSet; return *this; }
  bool operator()( const Point& p ) const
  { 
    ASSERT( ptrSet != 0 );
    return (*ptrSet)( p );
  }
private:
  const CharacteristicSet* ptrSet;
};

struct EigenVCM {
  Vector3 values;   //< eigenvalues
  Matrix33 vectors; //< eigenvectors
  unsigned char sigma[ 3 ];   //< permutation
  EigenVCM() {
    sigma[ 0 ] = 0;
    sigma[ 1 ] = 1;
    sigma[ 2 ] = 2;
  }
};

using namespace DGtal;

template <typename Domain>
struct RegularSubdivision {
  typedef typename Domain::Space Space;
  typedef typename Domain::Point Point;
  typedef typename Point::Coordinate Coordinate;
  typedef typename Domain::Vector Vector;
  typedef std::vector<Point> Storage;
  typedef ImageContainerBySTLVector<Domain,Storage*> StorageArray;

  RegularSubdivision( Clone<Domain> aDomain, Coordinate size )
    : myDomain( aDomain ), mySize( size ), myArray( aDomain )
  {
    Point dimensions = myDomain.upperBound() - myDomain.lowerBound();
    dimensions /= mySize;
    myReducedDomain = Domain( Point::zero, dimensions );
    myArray = StorageArray( myReducedDomain );
    // for ( typename Domain::ConstIterator it = myReducedDomain.begin(), itE = myReducedDomain.end();
    //       it != itE; ++it )
    //   myArray.setValue( *it, 0 );
  }

  ~RegularSubdivision()
  {
    for ( typename Domain::ConstIterator it = myReducedDomain.begin(), itE = myReducedDomain.end();
          it != itE; ++it )
      {
        Storage* ptr = myArray( *it );
        if ( ptr ) delete ptr;
      }
  }

  Point bucket( Point p ) const
  {
    p -= myDomain.lowerBound();
    return p / mySize;
  }

  /// Lowest point in bucket.
  Point inf( Point bucket ) const
  {
    bucket *= mySize;
    return bucket + myDomain.lowerBound();
  }
  /// Uppermost point in bucket.
  Point sup( Point bucket ) const
  {
    bucket *= mySize;
    return bucket + myDomain.lowerBound() + Point::diagonal(mySize-1);
  }

  void store( Point p ) 
  {
    Point b = bucket( p );
    Storage* pts = myArray( b );
    if ( pts == 0 ) 
      {
        pts = new Storage;
        myArray.setValue( b, pts );
      }
    pts->push_back( p );
  }

  template <typename PointConstIterator>
  void store( PointConstIterator it, PointConstIterator itE )
  {
    for ( ; it != itE; ++it )
      store( *it );
  }

  /**
     @tparam the type of a point predicate
     @param pred a predicate on point that must be convex.
  */
  template <typename PointPredicate>
  void getPoints( std::vector<Point> & pts, 
                  Point bucket_min, Point bucket_max, PointPredicate pred )
  {
    Domain local( bucket_min.sup( myReducedDomain.lowerBound() ),
                  bucket_max.inf( myReducedDomain.upperBound() ) );
    for ( typename Domain::ConstIterator it = local.begin(), itE = local.end(); it != itE; ++it )
      { // First check that at least one vertex satisfies the predicate
        Point lo = inf( *it );
        Point hi = sup( *it );
        if ( pred( lo ) || pred( hi ) 
             || pred( Point( hi[ 0 ], lo[ 1 ], lo[ 2 ] ) )
             || pred( Point( lo[ 0 ], hi[ 1 ], lo[ 2 ] ) )
             || pred( Point( hi[ 0 ], hi[ 1 ], lo[ 2 ] ) )
             || pred( Point( lo[ 0 ], lo[ 1 ], hi[ 2 ] ) )
             || pred( Point( lo[ 0 ], hi[ 1 ], hi[ 2 ] ) )
             || pred( Point( hi[ 0 ], lo[ 1 ], hi[ 2 ] ) ) ) {
          // std::cout << "Predicate ok for bucket " << *it << std::endl;
          const Storage* storage = myArray( *it );
          if ( storage )
            for ( typename Storage::const_iterator its = storage->begin(), itsE = storage->end(); its != itsE; ++its )
              if ( pred( *its ) ) pts.push_back( *its );
        }
      }
  }

  Domain myDomain;
  Coordinate mySize;
  Domain myReducedDomain;
  StorageArray myArray;
};

template <typename Point>
struct PointInBallPredicate {
  const Point center;
  const double radius2;
  PointInBallPredicate( Point aCenter, double aRadius )
    : center( aCenter ), radius2( sqr( aRadius ) )
  {}

  static inline double sqr( double x ) { return x*x; }
  bool operator()( const Point & p ) const
  {
    double d = sqr( center[ 0 ] - p[ 0 ] )
      + sqr( center[ 1 ] - p[ 1 ] )
      + sqr( center[ 2 ] - p[ 2 ] );
    return d <= radius2;
  }
  
};
// template <typename Metric, typename Point>
// struct PointInBallPredicate {
//   Metric metric;
//   Point center;
//   double radius;
//   PointInBallPredicate( Clone<Metric> aMetric, Point aCenter, double aRadius )
//     : metric( aMetric ), center( aCenter ), radius( aRadius )
//   {}
  
//   bool operator()( const Point & p ) const
//   {
//     double d = metric( center, p );
//     return d <= radius;
//   }
  
// };

template <typename OutputPointIterator, typename KSpace, typename Surfel>
void getPointsFromSurfel( OutputPointIterator outIt,
                          const KSpace & ks,
                          Surfel cell )
{
  Dimension k = ks.sOrthDir( cell );
  Dimension i = (k+1)%3;
  Dimension j = (i+1)%3;
  SCell l1 = ks.sIncident( cell, i, true );
  SCell l2 = ks.sIncident( cell, i, false );
  *outIt++ = ks.sCoords( ks.sIncident( l1, j, true ) );
  *outIt++ = ks.sCoords( ks.sIncident( l1, j, false ) );
  *outIt++ = ks.sCoords( ks.sIncident( l2, j, true ) );
  *outIt++ = ks.sCoords( ks.sIncident( l2, j, false ) );
}


using namespace DGtal;

template <typename SCell, typename RealVector>
struct GradientMapAdapter {
  typedef std::map<SCell,RealVector> SCell2RealVectorMap;
  typedef SCell                                 Argument;
  typedef RealVector                               Value;
  GradientMapAdapter( DGtal::ConstAlias<SCell2RealVectorMap> map )
    : myMap( map ) {}
  RealVector operator()( const Argument& arg ) const
  {
    typename SCell2RealVectorMap::const_iterator it = myMap->find( arg );
    if ( it != myMap->end() ) return it->second;
    else return RealVector();
  }
  DGtal::CountedConstPtrOrConstPtr<SCell2RealVectorMap> myMap;
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

  SCellEmbedderWithNormal( DGtal::ConstAlias<SCellEmbedder> embedder, 
                           DGtal::ConstAlias<SCell2RealVectorMap> map )
    : SCellEmbedder( embedder ), myMap( map )
  {}
  
  GradientMap gradientMap() const
  {
    return GradientMap( myMap );
  }

  DGtal::CountedConstPtrOrConstPtr<SCell2RealVectorMap> myMap;
};



template <typename Viewer, typename Surface, typename KernelFunction>
void computeSurfaceVCMFeatures( Viewer& viewer,
                                ostream & output_vcm,
                                const Surface & surface, 
                                double R, double r,
                                KernelFunction chi,
                                double trivial_r, int embedding,
                                double T )
{
  typedef typename Surface::DigitalSurfaceContainer DigitalSurfaceContainer;
  BOOST_CONCEPT_ASSERT(( CDigitalSurfaceContainer<DigitalSurfaceContainer> ));
  typedef typename Surface::KSpace KSpace; 
  typedef typename Surface::Surfel Surfel;
  typedef typename Surface::Cell Cell;
  typedef typename KSpace::Space Space;
  typedef typename Space::Point Point;
  typedef ExactPredicateLpSeparableMetric<Space, 2> Metric; // L2-metric
  typedef VoronoiCovarianceMeasure< Space, Metric > VCM;
  typedef typename VCM::Domain Domain;
  typedef typename Space::RealPoint RealPoint;
  typedef typename Space::RealVector RealVector;
  typedef typename Surface::ConstIterator SurfelConstIterator;
  // typedef HatPointFunction<Point,double> KernelFunction;
  typedef VoronoiCovarianceMeasureOnDigitalSurface< DigitalSurfaceContainer, Metric,
                                                    KernelFunction > VCMOnSurface;
  typedef typename VCMOnSurface::VectorN VectorN;
  typedef typename VCMOnSurface::Surfel2Normals::const_iterator S2NConstIterator;

  const KSpace & ks = surface.container().space();
  Surfel2PointEmbedding embType = 
    embedding == 0 ? Pointels :
    embedding == 1 ? InnerSpel :
    OuterSpel;
  // KernelFunction chi( 1.0, r );
  VCMOnSurface vcm_surface( surface, embType, R, r, chi, trivial_r, Metric(), true );

  trace.beginBlock ( "Export des normales." );
  VectorN lambda; // eigenvalues of chi-vcm
  DGtal::GradientColorMap<double> grad( 0, T );
  grad.addColor( DGtal::Color( 128, 128, 255 ) );
  grad.addColor( DGtal::Color( 128, 128, 255 ) );
  grad.addColor( DGtal::Color( 128, 128, 255 ) );
  grad.addColor( DGtal::Color( 128, 255, 255 ) );
  grad.addColor( DGtal::Color( 255, 255, 0 ) );
  grad.addColor( DGtal::Color( 255, 0, 0 ) );
  Cell dummy;
  viewer << SetMode3D( dummy.className(), "Basic" );

  for ( S2NConstIterator it = vcm_surface.mapSurfel2Normals().begin(), 
          itE = vcm_surface.mapSurfel2Normals().end(); it != itE; ++it )
    {
      Surfel s = it->first;
      const VectorN & n = it->second.vcmNormal;
      const VectorN & t = it->second.trivialNormal;
      vcm_surface.getChiVCMEigenvalues( lambda, s );
      double ratio = lambda[ 1 ] / ( lambda[ 0 ] + lambda[ 1 ] + lambda[ 2 ] ); // 3D !!!
      //std::cout << " " << ratio;
      if ( ratio > T ) ratio = T;
      Color c = grad( ratio );
      // bool sharp = ratio > T;
      output_vcm 
        << "CellN " << ks.sKCoord( s, 0 ) << " " << ks.sKCoord( s, 1 ) << " " << ks.sKCoord( s, 2 )
        // << " " << ks.sSign( s ) << (sharp ? " 1.0 0.0 0.0" : " 0.5 0.5 1.0" )
        << " " << ks.sSign( s ) 
        << " " << (((float)c.red())/255.0)
        << " " << (((float)c.green())/255.0)
        << " " << (((float)c.blue())/255.0)
        << " " << n[ 0 ] << " " << n[ 1 ] << " " << n[ 2 ] << std::endl;
      // viewer.setFillColor( sharp ? DGtal::Color( 255, 0.0, 0.0 ) : DGtal::Color( 128, 128, 255 ) );
      viewer.setFillColor( c );
      viewer << ks.unsigns( s );
    }
  trace.endBlock();
  
}

///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
  QApplication application(argc,argv);

  // parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("sdp,s", po::value<std::string>(), "name of the file containing 3d discrete points (.sdp) " )
    ("vol,v", po::value<std::string>(), "the volume file (.vol)" )
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min (excluded) to define binary shape" ) 
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max (included) to define binary shape" )
    ("R-radius,R", po::value<double>()->default_value( 5 ), "the parameter R in the VCM." )
    ("r-radius,r", po::value<double>()->default_value( 3 ), "the parameter r in the VCM." )
    ("kernel,k", po::value<std::string>()->default_value( "hat" ), "the function chi_r, either hat or ball." )
    ("trivial-radius,t", po::value<double>()->default_value( 3 ), "the parameter r for the trivial normal estimator." )
    ("angle-threshold,T", po::value<double>()->default_value( 1 ), "the angle threshold for feature detection.")
    ("embedding,E", po::value<int>()->default_value( 0 ), "the surfel -> point embedding: 0: Pointels, 1: InnerSpel, 2: OuterSpel." )
    ;  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  if( !parseOK || vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " -s [file.pts] -R 5\n"
		<< "Reads a set of points and computes a Voronoi map.";
      std::cout << "Usage: " << argv[0] << " -v [file.vol] -R 5\n"
		<< "Reads a vol file, extract a surface boundary and computes a Voronoi map."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "dvcm-3d -s ../helix-10.pts \n";
      return 0;
    }
  KSpace ks;
  std::vector<Surfel> surfels;   // Contains the surfels if data comes from volume.
  std::vector<Point> vectPoints; // Contains the set of discrete points.
  if ( vm.count("sdp") )
    {
      // Get points as a list of points in a file.
      vector<unsigned int> vPos;
      vPos.push_back(0);
      vPos.push_back(1);
      vPos.push_back(2);
      std::string inputSDP = vm["sdp"].as<std::string>();
      trace.info() << "Reading input 3d discrete points file: " << inputSDP; 
      vectPoints=  PointListReader<Point>::getPointsFromFile(inputSDP, vPos); 
      trace.info() << " [done] " << std::endl ; 
    }
  else if ( vm.count( "vol" ) )
    {
      trace.beginBlock( "Loading image into memory." );
      typedef DGtal::HyperRectDomain<Space> Domain;
      typedef DGtal::ImageSelector<Domain, unsigned char>::Type Image;
      string inputFilename = vm["vol"].as<std::string>();
      int thresholdMin = vm["thresholdMin"].as<int>();
      int thresholdMax = vm["thresholdMax"].as<int>();
      Image image = GenericReader<Image>::import (inputFilename );
      Domain domain = image.domain();
      typedef IntervalForegroundPredicate<Image> ThresholdedImage;
      ThresholdedImage thresholdedImage( image, thresholdMin, thresholdMax );
      Point dsize = domain.upperBound() - domain.lowerBound();
      trace.info() << "Image size = " << dsize[ 0 ]
                   << "x" << dsize[ 1 ]
                   << "x" << dsize[ 2 ] << std::endl;
      trace.endBlock();
      trace.beginBlock( "Extracting boundary by scanning the space. " );
      bool space_ok = ks.init( image.domain().lowerBound(),
                               image.domain().upperBound(), true );
      if (!space_ok)
        {
          trace.error() << "Error in the Khamisky space construction."<<std::endl;
          return 2;
        }
      typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
      typedef KSpace::Surfel Surfel;
      typedef ImplicitDigitalSurface< KSpace, ThresholdedImage > MySurfaceContainer;
      typedef DigitalSurface< MySurfaceContainer > MyDigitalSurface;
      MySurfelAdjacency surfAdj( true ); // interior in all directions.
      Surfel bel = Surfaces<KSpace>::findABel( ks, thresholdedImage, 10000 );
      MySurfaceContainer* container = 
        new MySurfaceContainer( ks, thresholdedImage, surfAdj, bel, false  );
      MyDigitalSurface surface( container ); //acquired
      trace.info() << "Digital surface has " << surface.size() << " surfels."
                   << std::endl;
      trace.endBlock();

      Viewer3D<> viewer( ks );
      viewer.setWindowTitle("Voronoi 3D viewer");
      viewer.show();

      std::ofstream output1( "titi-vcm-features.txt" );
      std::string kernel = vm[ "kernel" ].as<std::string>();
      if ( kernel == "hat" ) {
        typedef HatPointFunction<Point,double> KernelFunction;
        computeSurfaceVCMFeatures( viewer, output1, 
                                   surface, 
                                   vm["R-radius"].as<double>(), vm["r-radius"].as<double>(),
                                   KernelFunction( 1.0, vm["r-radius"].as<double>() ),
                                   vm["trivial-radius"].as<double>(),
                                   vm["embedding"].as<int>(),
                                   vm["angle-threshold"].as<double>() );
      } else if ( kernel == "ball" ) {
        typedef BallConstantPointFunction<Point,double> KernelFunction;
        computeSurfaceVCMFeatures( viewer, output1, 
                                   surface, 
                                   vm["R-radius"].as<double>(), vm["r-radius"].as<double>(),
                                   KernelFunction( 1.0, vm["r-radius"].as<double>() ),
                                   vm["trivial-radius"].as<double>(),
                                   vm["embedding"].as<int>(),
                                   vm["angle-threshold"].as<double>() );
      }
      output1.close();
      viewer << Viewer3D<>::updateDisplay;
      return application.exec();
    }
  else
    {
      trace.error() << " Input filename is required." << endl;      
    }
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
