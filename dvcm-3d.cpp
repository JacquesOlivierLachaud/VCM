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
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/geometry/volumes/estimation/VoronoiCovarianceMeasure.h"
#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/ElementaryConvolutionNormalVectorEstimator.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/GenericReader.h"


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
typedef DGtal::EigenDecomposition<3,double> LinearAlgebraTool;
typedef LinearAlgebraTool::Matrix Matrix33;
typedef LinearAlgebraTool::Vector Vector3;
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



template <typename Surface>
void computeSurfaceVCM( std::ostream& output_vcm, 
                        std::ostream& output_trivial, 
                        const Surface & surface, 
                        double R, double r,
                        double trivial_r, int embedding )
{
  typedef typename Surface::DigitalSurfaceContainer DigitalSurfaceContainer;
  BOOST_CONCEPT_ASSERT(( concepts::CDigitalSurfaceContainer<DigitalSurfaceContainer> ));
  typedef typename Surface::KSpace KSpace; 
  typedef typename Surface::Surfel Surfel;
  typedef typename KSpace::Space Space;
  typedef typename Space::Point Point;
  typedef ExactPredicateLpSeparableMetric<Space, 2> Metric; // L2-metric
  typedef VoronoiCovarianceMeasure< Space, Metric > VCM;
  typedef typename VCM::Domain Domain;
  typedef typename Space::RealPoint RealPoint;
  typedef typename Space::RealVector RealVector;
  typedef typename Surface::ConstIterator SurfelConstIterator;
  typedef functors::HatPointFunction<Point,double> KernelFunction;
  typedef VoronoiCovarianceMeasureOnDigitalSurface< DigitalSurfaceContainer, Metric,
                                                    KernelFunction > VCMOnSurface;
  typedef typename VCMOnSurface::VectorN VectorN;
  typedef typename VCMOnSurface::Surfel2Normals::const_iterator S2NConstIterator;

  const KSpace & ks = surface.container().space();
  Surfel2PointEmbedding embType = 
    embedding == 0 ? Pointels :
    embedding == 1 ? InnerSpel :
    OuterSpel;
  KernelFunction chi( 1.0, r );
  VCMOnSurface vcm_surface( surface, embType, R, r, chi, trivial_r, Metric(), true );

  trace.beginBlock ( "Export des normales." );
  std::map<Surfel,RealVector> surfel2normals;
  for ( S2NConstIterator it = vcm_surface.mapSurfel2Normals().begin(), 
          itE = vcm_surface.mapSurfel2Normals().end(); it != itE; ++it )
    {
      Surfel s = it->first;
      const VectorN & n = it->second.vcmNormal;
      const VectorN & t = it->second.trivialNormal;
      output_vcm 
        << "CellN " << ks.sKCoord( s, 0 ) << " " << ks.sKCoord( s, 1 ) << " " << ks.sKCoord( s, 2 )
        << " " << ks.sSign( s ) << " 0.5 0.5 1.0" 
        << " " << n[ 0 ] << " " << n[ 1 ] << " " << n[ 2 ] << std::endl;
      output_trivial
        << "CellN " << ks.sKCoord( s, 0 ) << " " << ks.sKCoord( s, 1 ) << " " << ks.sKCoord( s, 2 )
        << " " << ks.sSign( s ) << " 0.5 0.5 1.0" 
        << " " << t[ 0 ] << " " << t[ 1 ] << " " << t[ 2 ] << std::endl;
      surfel2normals[ s ] = n;
    }
  trace.endBlock();

  trace.beginBlock ( "Export OFF surface." );
  CanonicSCellEmbedder<KSpace> surfelEmbedder( ks );
  typedef SCellEmbedderWithNormal< CanonicSCellEmbedder<KSpace> > Embedder;
  Embedder embedder( surfelEmbedder, surfel2normals );
  std::ofstream output_off( "surface.off" );
  surface.exportAs3DNOFF( output_off, embedder );
  output_off.close();
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
    ("trivial-radius,t", po::value<double>()->default_value( 3 ), "the parameter r for the trivial normal estimator." )
    ("embedding,E", po::value<int>()->default_value( 0 ), "the surfel -> point embedding: 0: Pointels, 1: InnerSpel, 2: OuterSpel." )
    ("view,V", "view the set of digital points." )
    ("normal-dir,n", po::value<double>()->default_value( 2.0 ), "display the normal direction with the given size." )
    ("first-principal-dir,p", po::value<double>()->default_value( 0.0 ), "display the first principal direction with the given size." )
    ("second-principal-dir,q", po::value<double>()->default_value( 0.0 ), "display the second principal direction with the given size." )
    ("export,e", po::value<std::string>(), "exports surfel normals which can be viewed with viewSetOfSurfels." )
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
      typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
      ThresholdedImage thresholdedImage( image, thresholdMin, thresholdMax );
      trace.endBlock();
      trace.beginBlock( "Extracting boundary by scanning the space. " );
      bool space_ok = ks.init( image.domain().lowerBound(),
                               image.domain().upperBound(), true );
      if (!space_ok)
        {
          trace.error() << "Error in the Khamisky space construction."<<std::endl;
          return 2;
        }
      // std::back_insert_iterator< std::vector<Surfel> > outIt = std::back_inserter( surfels );
      // Surfaces<KSpace>::sWriteBoundary( outIt,
      //                                   ks, thresholdedImage,
      //                                   domain.lowerBound(),
      //                                   domain.upperBound() );
      // typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
      // typedef KSpace::SurfelSet SurfelSet;
      // typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
      // typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;
      // MySurfelAdjacency surfAdj( true ); // interior in all directions.
      // MySetOfSurfels* theSetOfSurfels = new MySetOfSurfels( ks, surfAdj );
      // for ( std::vector<Surfel>::const_iterator it = surfels.begin(), itE = surfels.end(); it != itE; ++it )
      //   theSetOfSurfels->surfelSet().insert( *it );
      // MyDigitalSurface surface( theSetOfSurfels ); //acquired
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

      std::ofstream output1( "titi-vcm.txt" );
      std::ofstream output2( "titi-triv.txt" );
      computeSurfaceVCM( output1, output2, surface, 
                         vm["R-radius"].as<double>(), vm["r-radius"].as<double>(),
                         vm["trivial-radius"].as<double>(),
                         vm["embedding"].as<int>() );
      output1.close();
      output2.close();
      return 0;
    }
  else
    {
      trace.error() << " Input filename is required." << endl;      
      return 0;
    }
  
  int R = (int) vm["R-radius"].as<double>();
  int r = (int) vm["r-radius"].as<double>();

  // Calcul le VCM.
  typedef VoronoiCovarianceMeasure< Space, Metric > VCM;
  Metric l2;
  VCM vcm( R, r, l2, true );
  vcm.init( vectPoints.begin(), vectPoints.end() );
  Domain domain = vcm.domain();  

  // // Calcule le domaine englobant.
  // Point lower = *vectPoints.begin();
  // Point upper = *vectPoints.begin();
  // for ( std::vector<Point>::const_iterator it = vectPoints.begin()+1, itE = vectPoints.end();
  //       it != itE; ++it )
  //   {
  //     lower = lower.inf( *it );
  //     upper = upper.sup( *it );
  //   }
  // lower -= Point::diagonal( R );
  // upper += Point::diagonal( R );
  // Domain domain( lower, upper );

  // // Définit les points de X.
  // CharacteristicSet charSet( domain );
  // for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
  //       it != itE; ++it )
  //   charSet.setValue( *it, true );

  // Display input set.
  Viewer3D<> viewer( ks );
  viewer.setWindowTitle("Voronoi 3D viewer");
  viewer.show();

  if ( vm.count( "view" ) )
    for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
          it != itE; ++it )
      viewer << *it;

  // // Le diagramme de Voronoi est calculé sur le complément de X.
  // CharacteristicSetPredicate inCharSet( charSet );
  // typedef NotPointPredicate<CharacteristicSetPredicate> NotPredicate;
  // NotPredicate notSetPred( inCharSet);

  // trace.beginBlock ( "Calcul du diagramme de Voronoi 3D" );
  // typedef VoronoiMap<Z3i::Space, NotPredicate, Metric > Voronoi3D;
  // Metric l2;
  // Voronoi3D voronoimap(domain,notSetPred,l2);
  // trace.endBlock();

  // // Calcul du VCM
  // trace.beginBlock ( "Calcul du VCM" );
  // RegularSubdivision<Domain> subdivision( domain, r ); // on subdivise en cellules de taille R.
  // typedef std::map<Point,Matrix33> Pt2VCM;
  // Pt2VCM pt2vcm;   // mapping point -> VCM
  // Matrix33 m; // zero matrix
  // for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
  //       it != itE; ++it )
  //   {
  //     pt2vcm[ *it ] = m;      // on initialise le VCM.
  //     subdivision.store( *it ); // on stocke les points dans la structure de subdivision.
  //   }
  // // On parcourt le domaine pour calculer le VCM.
  // int domain_size = voronoimap.domain().size();
  // int i = 0;
  // for(Voronoi3D::Domain::ConstIterator it = voronoimap.domain().begin(),
  //       itend = voronoimap.domain().end(); it != itend; ++it)
  //   {
  //     trace.progressBar(++i,domain_size);
  //     Point p = *it;
  //     Voronoi3D::Value q = voronoimap( p );   // site le plus proche de p
  //     if ( q != p )
  //       {
  //         double d = l2( q, p );
  //         if ( d <= (double)R ) // on se limite à l'offset de taille R
  //           { 
  //             RealVector v( p[ 0 ] - q[ 0 ], p[ 1 ] - q[ 1 ],p[ 2 ] - q[ 2 ] );
  //             // calcul le produit tensoriel V^t x V
  //             for ( Dimension i = 0; i < 3; ++i ) 
  //               for ( Dimension j = 0; j < 3; ++j )
  //                 m.setComponent( i, j, v[ i ] * v[ j ] ); 
  //             pt2vcm[ q ] += m;
  //           }
  //       }
  //   }
  // trace.endBlock();

  // trace.beginBlock ( "Intégration de VCM_r" );
  // // typedef PointInBallPredicate<Metric,Point> BallPredicate;
  // typedef PointInBallPredicate<Point> BallPredicate;
  // typedef std::map<Point,EigenVCM> Pt2EigenVCM;
  // Pt2EigenVCM pt2eigen_vcm;
  // std::vector<Point> neighbors;
  // int pts_size = vectPoints.size();
  // i = 0;
  // for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
  //       it != itE; ++it )
  //   {
  //     trace.progressBar(++i,pts_size);
  //     Point p = *it;
  //     BallPredicate pred( p, r ); // ( l2, p, r )
  //     Point bucket = subdivision.bucket( p ); 
  //     subdivision.getPoints( neighbors, 
  //                            bucket - Point::diagonal(1),
  //                            bucket + Point::diagonal(1),
  //                            pred );
  //     Matrix33 vcm;
  //     // std::cout << *it << " has " << neighbors.size() << " neighbors." << std::endl;
  //     for ( std::vector<Point>::const_iterator it_neighbors = neighbors.begin(),
  //             it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors )
  //       {
  //         Point q = *it_neighbors;
  //         double coef = ((double)r) - l2( p, q );
  //         Matrix33 vcm_q = pt2vcm[ q ];
  //         // Ne marche pas bien. Je ne sais pas pourquoi. Sans doute
  //         // que ce n'est pas la distance au germe, mais la distance
  //         // aux points de la cellule de Voronoi.
  //         vcm_q *= coef;
  //         vcm += vcm_q;
  //       }

  //     // On diagonalise le résultat.
  //     EigenVCM & evcm = pt2eigen_vcm[ p ];
  //     LinearAlgebraTool::getEigenDecomposition( vcm, evcm.vectors, evcm.values );

  //     neighbors.clear();
  //   }
  // trace.endBlock();

  trace.beginBlock ( "Intégration de VCM_r" );
  typedef std::map<Point,EigenVCM> Pt2EigenVCM;
  Pt2EigenVCM pt2eigen_vcm;
  int pts_size = vectPoints.size();
  int i = 0;
  functors::HatPointFunction< Point, double > chi_r( r, r );
  for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
        it != itE; ++it )
    {
      trace.progressBar(++i,pts_size);
      Point p = *it;
      Matrix33 m = vcm.measure( chi_r, p );
      // On diagonalise le résultat.
      EigenVCM & evcm = pt2eigen_vcm[ p ];
      LinearAlgebraTool::getEigenDecomposition( m, evcm.vectors, evcm.values );
    }
  trace.endBlock();


  trace.beginBlock ( "Calcul de l'application point -> surfels" );
  std::map< Point, std::set<Surfel > > pt2surfels;
  int surfels_size = surfels.size();
  for ( std::vector<Surfel>::const_iterator it = surfels.begin(), itE = surfels.end(); it != itE; ++it )
    {
      std::vector<Point> pts;
      Surfel surfel = *it;
      getPointsFromSurfel( std::back_inserter( pts ), ks, surfel );
      for ( int i = 0; i < 4; ++i )
        pt2surfels[ pts[ i ] ].insert( surfel );
    }
  trace.endBlock();

  // trace.beginBlock ( "Calcul de l'orientation moyenne des surfels" );
  // std::map< Surfel, Vector3 > trivial_normal;
  // i = 0; 
  // for ( std::vector<Surfel>::const_iterator it = surfels.begin(), itE = surfels.end(); it != itE; ++it )
  //   {
  //     trace.progressBar( ++i, surfels_size );
  //     std::vector<Point> pts;
  //     std::set<Surfel> neighbor_surfels;
  //     Surfel surfel = *it;
  //     getPointsFromSurfel( std::back_inserter( pts ), ks, surfel );
  //     Point p = pts[ 0 ];
  //     BallPredicate pred( p, r ); // ( l2, p, r )
  //     Point bucket = subdivision.bucket( p ); 
  //     subdivision.getPoints( neighbors, 
  //                            bucket - Point::diagonal(1),
  //                            bucket + Point::diagonal(1),
  //                            pred );
  //     Vector3 & n = trivial_normal[ surfel ];
      
  //     for ( std::vector<Point>::const_iterator it_neighbors = neighbors.begin(),
  //             it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors )
  //       {
  //         Point q = *it_neighbors;
  //         const std::set<Surfel> & pt_surfels = pt2surfels[ q ];
  //         neighbor_surfels.insert( pt_surfels.begin(), pt_surfels.end() );
  //       }
  //     neighbors.clear();
  //     for ( std::set<Surfel>::const_iterator it_surf_n = neighbor_surfels.begin(), 
  //             it_surf_end = neighbor_surfels.end(); it_surf_n != it_surf_end; ++it_surf_n )
  //       {
  //         Surfel neigh_surfel = *it_surf_n;
  //         Dimension k = ks.sOrthDir( neigh_surfel );
  //         bool dir =  ks.sDirect( neigh_surfel, k );
  //         Point out = ks.sCoords( ks.sIncident( neigh_surfel, k, ! dir ) );
  //         Point in = ks.sCoords( ks.sIncident( neigh_surfel, k, dir ) );
  //         Vector trivial = out - in;
  //         n += Vector3( trivial[ 0 ], trivial[ 1 ], trivial[ 2 ] );
  //       }
  //   }
  // trace.endBlock();

  // trace.beginBlock ( "Relaxation des vecteurs propres" );
  // i = 0;
  // for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
  //       it != itE; ++it )
  //   {
  //     trace.progressBar(++i,pts_size);
  //     Point p = *it;
  //     BallPredicate pred( p, 2 ); // look at neighborhood  // ( l2, p, r )
  //     Point bucket = subdivision.bucket( p ); 
  //     subdivision.getPoints( neighbors, 
  //                            bucket - Point::diagonal(1),
  //                            bucket + Point::diagonal(1),
  //                            pred );
  //     EigenVCM & evcm = pt2eigen_vcm[ p ];
  //     Matrix33 alignment;
  //     for ( std::vector<Point>::const_iterator it_neighbors = neighbors.begin(),
  //             it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors )
  //       {
  //         Point q = *it_neighbors;
  //         const EigenVCM & evcm_neighbor = pt2eigen_vcm[ q ];
  //         for ( int j = 0; j < 3; ++j )
  //           for ( int k = 0; k < 3; ++k )
  //             {
  //               double val = alignment( j, k );
  //               alignment.setComponent( j, k, 
  //                                       val + abs( evcm.vectors.column( j ).dot( evcm_neighbor.vectors.column( k ) ) ) );
  //             }
  //       }
  //     neighbors.clear();
  //     for ( int j = 0; j < 3; ++j )
  //       {
  //         int best = 0;
  //         for ( int k = 1; k < 3; ++k )
  //           if ( alignment( j, k ) > alignment( j, best ) )
  //             best = k;
  //         evcm.sigma[ j ] = best;
  //         if ( ( j == 2 ) && ( j != best ) )
  //           std::cout << "Point " << p << " is relaxed: " << j << " -> " << best << std::endl;
  //       }
  //   }
   
  // trace.endBlock();

  trace.beginBlock ( "Affichage des normales" );
  double size_n = vm[ "normal-dir" ].as<double>();
  double size_p1 = vm[ "first-principal-dir" ].as<double>();
  double size_p2 = vm[ "second-principal-dir" ].as<double>();
  for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
        it != itE; ++it )
    {
      Point p = *it;
      const EigenVCM & evcm = pt2eigen_vcm[ p ];
      RealPoint rp( p[ 0 ], p[ 1 ], p[ 2 ] );
      // last eigenvalue is the greatest.
      if ( size_n != 0.0 ) {
        Vector3 n = evcm.vectors.column( evcm.sigma[ 2 ] ); // troisième colonne
        n *= size_n;
        viewer.setLineColor( Color::Black );
        viewer.addLine( rp + n, rp - n, 0.1 );
      }
      if ( size_p1 != 0.0 ) {
        Vector3 n = evcm.vectors.column( evcm.sigma[ 1 ] ); // deuxième colonne
        n *= size_p1;
        viewer.setLineColor( Color::Blue );
        viewer.addLine( rp + n, rp - n, 0.1 );
      }
      if ( size_p2 != 0.0 ) {
        Vector3 n = evcm.vectors.column( evcm.sigma[ 0 ] ); // première colonne
        n *= size_p2;
        viewer.setLineColor( Color::Red );
        viewer.addLine( rp + n, rp - n, 0.1 );
      }

    }
  trace.endBlock();

  if ( vm.count( "export" ) )
    {
      std::string filename = vm[ "export" ].as<std::string>();
      std::ofstream output( filename.c_str() );
      std::set<Surfel> surfelSet;
      for ( std::vector<Surfel>::const_iterator it = surfels.begin(), itE = surfels.end(); it != itE; ++it )
        surfelSet.insert( *it );
      for ( std::set<Surfel>::const_iterator it = surfelSet.begin(), itE = surfelSet.end(); it != itE; ++it )
        {
          Dimension k = ks.sOrthDir( *it );
          Dimension i = (k+1)%3;
          Dimension j = (i+1)%3;
          SCell l1 = ks.sIncident( *it, i, true );
          SCell l2 = ks.sIncident( *it, i, false );
          Point p1 = ks.sCoords( ks.sIncident( l1, j, true ) );
          Point p2 = ks.sCoords( ks.sIncident( l1, j, false ) );
          Point p3 = ks.sCoords( ks.sIncident( l2, j, true ) );
          Point p4 = ks.sCoords( ks.sIncident( l2, j, false ) );
          const EigenVCM & evcm1 = pt2eigen_vcm[ p1 ];
          const EigenVCM & evcm2 = pt2eigen_vcm[ p2 ];
          const EigenVCM & evcm3 = pt2eigen_vcm[ p3 ];
          const EigenVCM & evcm4 = pt2eigen_vcm[ p4 ];
          bool dir =  ks.sDirect( *it, k );
          Point out = ks.sCoords( ks.sIncident( *it, k, ! dir ) );
          Point in = ks.sCoords( ks.sIncident( *it, k, dir ) );
          Vector trivial = out - in;
          Vector3 t( trivial[ 0 ], trivial[ 1 ], trivial[ 2 ] );
          // const Vector3 & t = trivial_normal[ *it ];
          Vector3 n1 = evcm1.vectors.column( evcm1.sigma[ 2 ] );
          if ( n1.dot( t ) > 0 ) n1 = -n1;
          Vector3 n2 = evcm2.vectors.column( evcm2.sigma[ 2 ] );
          if ( n2.dot( t ) > 0 ) n2 = -n2;
          Vector3 n3 = evcm3.vectors.column( evcm3.sigma[ 2 ] );
          if ( n3.dot( t ) > 0 ) n3 = -n3;
          Vector3 n4 = evcm4.vectors.column( evcm4.sigma[ 2 ] );
          if ( n4.dot( t ) > 0 ) n4 = -n4;
          n1 += n2; n1 += n3; n1 += n4;
          n1 /= 4.0;
          //if ( n1.dot( t ) > 0 ) n1 = -n1;
          output << "CellN " << ks.sKCoord( *it, 0 ) 
                 << " " << ks.sKCoord( *it, 1 )
                 << " " << ks.sKCoord( *it, 2 )
                 << " " << ks.sSign( *it )
                 << " 0.5 0.5 1.0" 
                 << " " << n1[ 0 ] << " " << n1[ 1 ] << " " << n1[ 2 ] << std::endl;
        }
      output.close();
    }
  // On affiche le vecteur vers le site le plus proche seulement si il est à distance <= 4.
  // for(Voronoi3D::Domain::ConstIterator it = voronoimap.domain().begin(),
  //     itend = voronoimap.domain().end(); it != itend; ++it)
  // {
  //   Point p = *it;
  //   Voronoi3D::Value q = voronoimap( p );   // site le plus proche de p
  //   if ( q != p )
  //     {
  //       double d = l2( q, p );
  //       if ( ( (double)(R-1) <= d ) && ( d <= (double)R ) ) // on affiche que la dernière couche.
  //       { 
  //         viewer.addLine( RealPoint( p[ 0 ], p[ 1 ], p[ 2 ] ),
  //                         RealPoint( q[ 0 ], q[ 1 ], q[ 2 ] ),
  //                         0.03 ); // width
  //       }
  //     }
  // }
  viewer << Viewer3D<>::updateDisplay;
  return application.exec();
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
