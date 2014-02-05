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
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/io/viewers/Viewer3D.h"

//#include "DGtal/io/colormaps/HueShadeColorMap.h"
//#include "DGtal/io/boards/Board2D.h"
//! [Voro2D-header]

///////////////////////////////////////////////////////////////////////////////

using namespace std;

typedef DGtal::Z3i::Space Space;
typedef DGtal::Z3i::Vector Vector;
typedef DGtal::Z3i::Point Point;
typedef DGtal::Z3i::RealPoint RealPoint;
typedef DGtal::Z3i::RealVector RealVector;
typedef DGtal::HyperRectDomain<Space> Domain;
typedef DGtal::ImageContainerBySTLVector<Domain,bool> CharacteristicSet;
typedef DGtal::ExactPredicateLpSeparableMetric<Space, 3> Metric; // L2-metric
typedef DGtal::SimpleMatrix<double,3,3> Matrix3;

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
    ("input,i", po::value<std::string>(), "name of the file containing 3d discrete points (.sdp) " )
    ("big-radius,R", po::value<int>()->default_value( 4 ), "the parameter R in the VCM." )
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
      std::cout << "Usage: " << argv[0] << " -i [file.sdp]\n"
		<< "Reads a set of points and computes a Voronoi map."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "dvcm-3d -i ellipse.sdp \n";
      return 0;
    }
  if(! vm.count("input"))
    {
      trace.error() << " Input filename is required." << endl;      
      return 0;
    }
  
  int R = vm["big-radius"].as<int>();

  // Récupère les points.
  vector<unsigned int> vPos;
  vPos.push_back(0);
  vPos.push_back(1);
  vPos.push_back(2);
  std::string inputSDP = vm["input"].as<std::string>();
  trace.info() << "Reading input 3d discrete points file: " << inputSDP; 
  std::vector<Point> vectPoints=  PointListReader<Point>::getPointsFromFile(inputSDP, vPos); 
  trace.info() << " [done] " << std::endl ; 

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

  trace.beginBlock ( "Calcul du VCM" );
  // Calcul du VCM
  typedef std::map<Point,Matrix3> PT2VCM;
  PT2VCM pt2vcm;   // mapping point -> VCM
  Matrix3 m; // zero matrix
  for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
        it != itE; ++it )
    pt2vcm[ *it ] = m;
  // On parcourt le domaine pour calculer le VCM.
  for(Voronoi3D::Domain::ConstIterator it = voronoimap.domain().begin(),
        itend = voronoimap.domain().end(); it != itend; ++it)
    {
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
  // On diagonalise le VCM. ATTENTION il faudrait aussi utiliser le petit r !
  for ( std::vector<Point>::const_iterator it = vectPoints.begin(), itE = vectPoints.end();
        it != itE; ++it )
    {
      Matrix3 & vcm = pt2vcm[ *it ];
      // a diagonaliser ...
    }
  trace.endBlock();

  
  // On affiche le vecteur vers le site le plus proche seulement si il est à distance <= 4.
  for(Voronoi3D::Domain::ConstIterator it = voronoimap.domain().begin(),
      itend = voronoimap.domain().end(); it != itend; ++it)
  {
    Point p = *it;
    Voronoi3D::Value q = voronoimap( p );   // site le plus proche de p
    if ( q != p )
      {
        double d = l2( q, p );
        if ( ( (double)(R-1) <= d ) && ( d <= (double)R ) ) // on affiche que la dernière couche.
        { 
          viewer.addLine( RealPoint( p[ 0 ], p[ 1 ], p[ 2 ] ),
                          RealPoint( q[ 0 ], q[ 1 ], q[ 2 ] ),
                          0.03 ); // width
        }
      }
  }
  viewer << Viewer3D<>::updateDisplay;
  return application.exec();
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
