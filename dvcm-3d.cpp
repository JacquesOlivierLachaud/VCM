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
#include "DGtal/helpers/StdDefs.h"

//! [Voro2D-header]
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/SimpleMatrix.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/IntervalForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <DGtal/io/readers/GenericReader.h>
#include "DGtal/math/EigenValues3D.h"
#include "DGtal/topology/helpers/Surfaces.h"

//#include "DGtal/io/colormaps/HueShadeColorMap.h"
//#include "DGtal/io/boards/Board2D.h"
//! [Voro2D-header]

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
typedef DGtal::ExactPredicateLpSeparableMetric<Space, 3> Metric; // L2-metric
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
              pts.push_back( *its );
        }
      }
  }

  Domain myDomain;
  Coordinate mySize;
  Domain myReducedDomain;
  StorageArray myArray;
};

template <typename Metric, typename Point>
struct PointInBallPredicate {
  Metric metric;
  Point center;
  double radius;
  PointInBallPredicate( Clone<Metric> aMetric, Point aCenter, double aRadius )
    : metric( aMetric ), center( aCenter ), radius( aRadius )
  {}
  
  bool operator()( const Point & p ) const
  {
    double d = metric( center, p );
    return d <= radius;
  }
  
};

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
    ("R-radius,R", po::value<int>()->default_value( 5 ), "the parameter R in the VCM." )
    ("r-radius,r", po::value<int>()->default_value( 3 ), "the parameter r in the VCM." )
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
      typedef IntervalForegroundPredicate<Image> ThresholdedImage;
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
      std::back_insert_iterator< std::vector<Surfel> > outIt = std::back_inserter( surfels );
      Surfaces<KSpace>::sWriteBoundary( outIt,
                                        ks, thresholdedImage,
                                        domain.lowerBound(),
                                        domain.upperBound() );
      trace.info() << "Digital surface has " << surfels.size() << " surfels."
		 << std::endl;
      std::set<Point> pointSet;
      for ( std::vector<Surfel>::const_iterator it = surfels.begin(), itE = surfels.end(); it != itE; ++it )
        {
          Dimension k = ks.sOrthDir( *it );
          Dimension i = (k+1)%3;
          Dimension j = (i+1)%3;
          SCell l1 = ks.sIncident( *it, i, true );
          SCell l2 = ks.sIncident( *it, i, false );
          pointSet.insert( ks.sCoords( ks.sIncident( l1, j, true ) ) );
          pointSet.insert( ks.sCoords( ks.sIncident( l1, j, false ) ) );
          pointSet.insert( ks.sCoords( ks.sIncident( l2, j, true ) ) );
          pointSet.insert( ks.sCoords( ks.sIncident( l2, j, false ) ) );
        }
      vectPoints.resize( pointSet.size() );
      std::copy( pointSet.begin(), pointSet.end(), vectPoints.begin() );
      pointSet.clear();
      trace.endBlock();
    }
  else
    {
      trace.error() << " Input filename is required." << endl;      
      return 0;
    }
  
  int R = vm["R-radius"].as<int>();
  int r = vm["r-radius"].as<int>();

  // // Récupère les points.
  // vector<unsigned int> vPos;
  // vPos.push_back(0);
  // vPos.push_back(1);
  // vPos.push_back(2);
  // std::string inputSDP = vm["input"].as<std::string>();
  // trace.info() << "Reading input 3d discrete points file: " << inputSDP; 
  // std::vector<Point> vectPoints=  PointListReader<Point>::getPointsFromFile(inputSDP, vPos); 
  // trace.info() << " [done] " << std::endl ; 

  // Calcule le domaine englobant.
  Point lower = *vectPoints.begin();
  Point upper = *vectPoints.begin();
  for ( std::vector<Point>::const_iterator it = vectPoints.begin()+1, itE = vectPoints.end();
        it != itE; ++it )
    {
      lower = lower.inf( *it );
      upper = upper.sup( *it );
    }
  lower -= Point::diagonal( R );
  upper += Point::diagonal( R );
  Domain domain( lower, upper );

  // Définit les points de X.
  CharacteristicSet charSet( domain );
  for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
        it != itE; ++it )
    charSet.setValue( *it, true );

  // Display input set.
  Viewer3D<> viewer;
  viewer.setWindowTitle("Voronoi 3D viewer");
  viewer.show();

  if ( vm.count( "view" ) )
    for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
          it != itE; ++it )
      viewer << *it;

  // Le diagramme de Voronoi est calculé sur le complément de X.
  CharacteristicSetPredicate inCharSet( charSet );
  typedef NotPointPredicate<CharacteristicSetPredicate> NotPredicate;
  NotPredicate notSetPred( inCharSet);

  trace.beginBlock ( "Calcul du diagramme de Voronoi 3D" );
  typedef VoronoiMap<Z3i::Space, NotPredicate, Metric > Voronoi3D;
  Metric l2;
  Voronoi3D voronoimap(domain,notSetPred,l2);
  trace.endBlock();

  // Calcul du VCM
  trace.beginBlock ( "Calcul du VCM" );
  RegularSubdivision<Domain> subdivision( domain, r ); // on subdivise en cellules de taille R.
  typedef std::map<Point,Matrix33> Pt2VCM;
  Pt2VCM pt2vcm;   // mapping point -> VCM
  Matrix33 m; // zero matrix
  for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
        it != itE; ++it )
    {
      pt2vcm[ *it ] = m;      // on initialise le VCM.
      subdivision.store( *it ); // on stocke les points dans la structure de subdivision.
    }
  // On parcourt le domaine pour calculer le VCM.
  int domain_size = voronoimap.domain().size();
  int i = 0;
  for(Voronoi3D::Domain::ConstIterator it = voronoimap.domain().begin(),
        itend = voronoimap.domain().end(); it != itend; ++it)
    {
      trace.progressBar(++i,domain_size);
      Point p = *it;
      Voronoi3D::Value q = voronoimap( p );   // site le plus proche de p
      if ( q != p )
        {
          double d = l2( q, p );
          if ( d <= (double)R ) // on se limite à l'offset de taille R
            { 
              RealVector v( p[ 0 ] - q[ 0 ], p[ 1 ] - q[ 1 ],p[ 2 ] - q[ 2 ] );
              // calcul le produit tensoriel V^t x V
              for ( Dimension i = 0; i < 3; ++i ) 
                for ( Dimension j = 0; j < 3; ++j )
                  m.setComponent( i, j, v[ i ] * v[ j ] ); 
              pt2vcm[ q ] += m;
            }
        }
    }
  trace.endBlock();

  trace.beginBlock ( "Intégration de VCM_r" );
  typedef PointInBallPredicate<Metric,Point> BallPredicate;
  typedef std::map<Point,EigenVCM> Pt2EigenVCM;
  Pt2EigenVCM pt2eigen_vcm;
  std::vector<Point> neighbors;
  int pts_size = vectPoints.size();
  i = 0;
  for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
        it != itE; ++it )
    {
      trace.progressBar(++i,pts_size);
      Point p = *it;
      BallPredicate pred( l2, p, r );
      Point bucket = subdivision.bucket( p ); 
      subdivision.getPoints( neighbors, 
                             bucket - Point::diagonal(1),
                             bucket + Point::diagonal(1),
                             pred );
      Matrix33 vcm;
      // std::cout << *it << " has " << neighbors.size() << " neighbors." << std::endl;
      for ( std::vector<Point>::const_iterator it_neighbors = neighbors.begin(),
              it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors )
        {
          Point q = *it_neighbors;
          double coef = ((double)r) - l2( p, q );
          Matrix33 vcm_q = pt2vcm[ q ];
          // Ne marche pas bien. Je ne sais pas pourquoi. Sans doute
          // que ce n'est pas la distance au germe, mais la distance
          // aux points de la cellule de Voronoi.
          // vcm_q *= coef;
          vcm += vcm_q;
        }

      // On diagonalise le résultat.
      EigenVCM & evcm = pt2eigen_vcm[ p ];
      LinearAlgebraTool::getEigenDecomposition( vcm, evcm.vectors, evcm.values );

      neighbors.clear();
    }
  trace.endBlock();

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
        Vector3 n = evcm.vectors.column( 2 ); // troisième colonne
        n *= size_n;
        viewer.setLineColor( Color::Black );
        viewer.addLine( rp + n, rp - n, 0.1 );
      }
      if ( size_p1 != 0.0 ) {
        Vector3 n = evcm.vectors.column( 1 ); // deuxième colonne
        n *= size_p1;
        viewer.setLineColor( Color::Blue );
        viewer.addLine( rp + n, rp - n, 0.1 );
      }
      if ( size_p2 != 0.0 ) {
        Vector3 n = evcm.vectors.column( 0 ); // première colonne
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
      for ( std::vector<Surfel>::const_iterator it = surfels.begin(), itE = surfels.end(); it != itE; ++it )
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
          Vector3 n = pt2eigen_vcm[ p1 ].vectors.column( 2 )
            + pt2eigen_vcm[ p2 ].vectors.column( 2 )
            + pt2eigen_vcm[ p3 ].vectors.column( 2 )
            + pt2eigen_vcm[ p4 ].vectors.column( 2 );
          n /= 4.0;
          bool orth = ks.sDirect( *it, k );
          Point out = ks.sCoords( ks.sIndirectIncident( *it, k ) );
          Point in = ks.sCoords( ks.sDirectIncident( *it, k ) );
          Vector trivial = out - in;
          Vector3 t( trivial[ 0 ], trivial[ 1 ], trivial[ 2 ] );
          if ( n.dot( t ) > 0 ) n = -n;
          output << "CellN " << ks.sKCoord( *it, 0 ) 
                 << " " << ks.sKCoord( *it, 1 )
                 << " " << ks.sKCoord( *it, 2 )
                 << " " << ks.sSign( *it )
                 << " 0.5 0.5 1.0" 
                 << " " << n[ 0 ] << " " << n[ 1 ] << " " << n[ 2 ] << std::endl;
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
