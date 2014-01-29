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
 * @file voronoimap2D.cpp
 * @ingroup Examples
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2013/01/17
 *
 * An example file named voronoimap2D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

//! [Voro2D-header]
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"

#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/boards/Board2D.h"
//! [Voro2D-header]

///////////////////////////////////////////////////////////////////////////////

using namespace std;

typedef DGtal::Z2i::Space Space;
typedef DGtal::Z2i::Vector Vector;
typedef DGtal::Z2i::Point Point;
typedef DGtal::HyperRectDomain<Space> Domain;
typedef DGtal::ImageContainerBySTLVector<Domain,bool> CharacteristicSet;
typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> Metric; // L2-metric

// Model of CPointPredicate
struct CharacteristicSetPredicate {
  typedef CharacteristicSetPredicate Self;
  typedef DGtal::Z2i::Point Point;
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

  // parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "name of the file containing 2d discrete points (.sdp) " )
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
		<< "dvcm-2d -i ellipse.sdp \n";
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
  std::string inputSDP = vm["input"].as<std::string>();
  trace.info() << "Reading input 2d discrete points file: " << inputSDP; 
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
  Board2D board;
  for ( Domain::ConstIterator it = domain.begin(), itE = domain.end();
        it != itE; ++it )
    if ( charSet( *it ) ) board << *it;

  // Le diagramme de Voronoi est calculé sur le complément de X.
  CharacteristicSetPredicate inCharSet( charSet );
  typedef NotPointPredicate<CharacteristicSetPredicate> NotPredicate;
  NotPredicate notSetPred( inCharSet);

  trace.beginBlock ( "Calcul du diagramme de Voronoi 2D" );
  typedef VoronoiMap<Z2i::Space, NotPredicate, Metric > Voronoi2D;
  Metric l2;
  Voronoi2D voronoimap(domain,notSetPred,l2);
  trace.endBlock();

  // On affiche le vecteur vers le site le plus proche seulement si il est à distance <= 4.
  // board.clear();
  // board << domain;
  for(Voronoi2D::Domain::ConstIterator it = voronoimap.domain().begin(),
      itend = voronoimap.domain().end(); it != itend; ++it)
  {
    Voronoi2D::Value site = voronoimap( *it );   //closest site to (*it)
    if (site != (*it))
      {
        double d = l2( site, *it );
        if( d <= (double)R )
        { 
            Display2DFactory::draw( board,   site - (*it), (*it)); //Draw an arrow
        }
      }
  }
  board.saveSVG("voronoimap-voro-R.svg");

  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
