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
 * @file delta-distance.cpp
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/11/16
 *
 * Computes the delta-distance to a gray-level image.
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
#include "DGtal/base/ConstAlias.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


template <typename Distance>
struct DistanceToPointFunctor {

  typedef typename Distance::Space   Space;
  typedef typename Distance::Value   Value;
  typedef typename Space::Point      Point;
  
  Point p;
  DistanceToPointFunctor( Clone<Distance> distance,
                          const Point& aP )
    : myDistance( distance ), p( aP ) {}

  Value operator()( const Point& q ) const
  {
    return myDistance( p, q );
  }
  Distance myDistance;
};

// A measure is a function
template <typename ImageFct>
class DistanceToMeasure {
public:
  typedef typename ImageFct::Value   Value;
  typedef typename ImageFct::Point   Point;
  typedef typename ImageFct::Domain  Domain;
  typedef typename Domain::Space     Space;
  typedef typename Space::RealVector RealVector;

public:
  DistanceToMeasure( Value m0, const ImageFct& measure, Value rmax = 10.0 )
    : myMass( m0 ), myMeasure( measure ), myDistance2( myMeasure.domain() ),
      myR2Max( rmax*rmax )
  {
    init( myMeasure );
  }
  
  void init( const ImageFct& measure )
  {
    double       nb = myDistance2.domain().size();
    unsigned int i  = 0;
    trace.progressBar( i, nb );

    for ( typename Domain::ConstIterator it = myDistance2.domain().begin(),
            itE = myDistance2.domain().end(); it != itE; ++it, ++i )
      {
        if ( ( i % 100 ) == 0 ) trace.progressBar( i, nb );
        myDistance2.setValue( *it, computeDistance2( *it ) );
      }
  }

  Value distance2( const Point& p ) const
  {
    return myDistance2( p );
  }

  Value safeDistance2( const Point& p ) const
  {
    if ( myDistance2.domain().isInside( p ) )
      return myDistance2( p );
    else return myDistance2( box( p ) );
  }
  
  Point box( const Point& p ) const
  {
    Point q = p.sup( myDistance2.domain().lowerBound() );
    return q.inf( myDistance2.domain().upperBound() );
  }

  RealVector projection( const Point& p ) const
  {
    Point p_left = box( p - Point( 1, 0 ) );
    Point p_right = box( p + Point( 1, 0 ) );
    Point p_down = box( p - Point( 0, 1 ) );
    Point p_up = box( p + Point( 0, 1 ) );
    Value d2_center = distance2( p );
    Value d2_left = distance2( p_left );
    Value d2_right = distance2( p_right );
    Value d2_down = distance2( p_down );
    Value d2_up = distance2( p_up );
    // Value gx = 
    //   // std::min( ( d2_right - d2_left ) / ( p_right[ 0 ] - p_left[ 0 ] ),
    //   std::min( ( d2_right - d2_center ) / ( p_right[ 0 ] - p[ 0 ] ),
    //             ( d2_center - d2_left ) / ( p[ 0 ] - p_left[ 0 ] ) ); //  );
    // Value gy = 
    //   // std::min( ( d2_up - d2_down ) / ( p_up[ 1 ] - p_down[ 1 ] ),
    //   std::min( ( d2_up - d2_center ) / ( p_up[ 1 ] - p[ 1 ] ),
    //             ( d2_center - d2_down ) / ( p[ 1 ] - p_down[ 1 ] ) ); // );
    bool right = abs( d2_right - d2_center ) >= abs( d2_center - d2_left );
    bool up    = abs( d2_up    - d2_center ) >= abs( d2_center - d2_down );
    Value gx = right ? ( d2_right - d2_center ) : ( d2_center - d2_left );
    Value gy = up    ? ( d2_up    - d2_center ) : ( d2_center - d2_down );
    return RealVector( -gx / 2.0, -gy / 2.0 );
    // Value gx = (distance2( px2 ) - distance2( px1 ))
    //   / ( 2.0 * ( px2[ 0 ] - px1[ 0 ] ) );
    // Value gy = (distance2( py2 ) - distance2( py1 ))
    //   / ( 2.0 * ( py2[ 1 ] - py1[ 1 ] ) );
    // return RealVector( -gx, -gy );
  }
  
  Value computeDistance2( const Point& p )
  {
    typedef ExactPredicateLpSeparableMetric<Space,2> Distance;
    typedef DistanceToPointFunctor<Distance>         DistanceToPoint;
    typedef MetricAdjacency<Space, 1>                Graph;
    typedef DistanceBreadthFirstVisitor< Graph, DistanceToPoint, std::set<Point> > 
      DistanceVisitor;
    typedef typename DistanceVisitor::Node MyNode;
    typedef typename DistanceVisitor::Scalar MySize;

    Value             m  = NumberTraits<Value>::ZERO;
    Value             d2 = NumberTraits<Value>::ZERO;
    Graph             graph;
    DistanceToPoint   d2pfct( Distance(), p );
    DistanceVisitor   visitor( graph, d2pfct, p );

    unsigned long nbSurfels = 0;
    Value last = d2pfct( p );
    MyNode node;
    while ( ! visitor.finished() )
    {
      node = visitor.current();
      if ( ( node.second != last ) // all the vertices of the same layer have been processed. 
           && ( m >= myMass ) ) break;
      if ( node.second > myR2Max ) { d2 = m * myR2Max; break; }
      if ( myMeasure.domain().isInside( node.first ) )
        {
          Value mpt  = myMeasure( node.first );
          d2        += mpt * node.second * node.second; 
          m         += mpt;
          last       = node.second;
          visitor.expand();
        }
      else
        visitor.ignore();
    }
    return d2 / m;
  }

public:
  Value myMass;
  const ImageFct& myMeasure;
  ImageFct myDistance2;
  Value myR2Max;
};

int main( int argc, char** argv )
{
  using namespace DGtal;
  using namespace DGtal::Z2i;
  
  typedef ImageContainerBySTLVector<Domain,unsigned char> GrayLevelImage2D;
  typedef ImageContainerBySTLVector<Domain,float>         FloatImage2D;
  typedef DistanceToMeasure<FloatImage2D>                 Distance;
  if ( argc <= 3 ) return 1;
  GrayLevelImage2D img  = GenericReader<GrayLevelImage2D>::import( argv[ 1 ] );
  double           mass = atof( argv[ 2 ] );
  double           rmax = atof( argv[ 3 ] );
  FloatImage2D     fimg( img.domain() );
  FloatImage2D::Iterator outIt = fimg.begin();
  for ( GrayLevelImage2D::ConstIterator it = img.begin(), itE = img.end();
        it != itE; ++it )
    {
      float v = ((float)*it) / 255.0;
      *outIt++ = v;
    }
  trace.beginBlock( "Computing delta-distance." );
  Distance     delta( mass, fimg, rmax );
  const FloatImage2D& d2 = delta.myDistance2;
  trace.endBlock();

  float m = 0.0f;
  for ( typename Domain::ConstIterator it = d2.domain().begin(),
          itE = d2.domain().end(); it != itE; ++it )
    {
      Point p = *it;
      float v = sqrt( d2( p ) );
      m = std::max( v, m );
    }

  GradientColorMap<float> cmap_grad( 0, m );
  cmap_grad.addColor( Color( 255, 255, 255 ) );
  cmap_grad.addColor( Color( 255, 255, 0 ) );
  cmap_grad.addColor( Color( 255, 0, 0 ) );
  cmap_grad.addColor( Color( 0, 255, 0 ) );
  cmap_grad.addColor( Color( 0,   0, 255 ) );
  cmap_grad.addColor( Color( 0,   0, 0 ) );
  Board2D board;
  board << SetMode( d2.domain().className(), "Paving" );
  

  for ( typename Domain::ConstIterator it = d2.domain().begin(),
          itE = d2.domain().end(); it != itE; ++it )
    {
      Point p = *it;
      float v = sqrt( d2( p ) );
      v = std::min( (float)m, std::max( v, 0.0f ) ); 
      board << CustomStyle( p.className(),
                            new CustomColors( Color::Black, cmap_grad( v ) ) )
            << p;

      RealVector grad = delta.projection( p );
      // / ( 1.1 - ( (double)img( *it ) ) / 255.0 ) ;
      board.drawLine( p[ 0 ], p[ 1 ], p[ 0 ] + grad[ 0 ], p[ 1 ] + grad[ 1 ], 0 );
    }
  std::cout << endl;
  board.saveEPS("delta2.eps");
  return 0;
}
