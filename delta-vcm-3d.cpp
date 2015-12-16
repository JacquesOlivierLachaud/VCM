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
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/kernel/Point2ScalarFunctors.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

// work-arounds for DGtal
namespace DGtal {
  typedef SimpleMatrix< double, 3, 3 > MatrixDouble;
  bool operator!=( const MatrixDouble& m1, const MatrixDouble& m2 )
  { return ! ( m1 == m2 ); }
  typedef SimpleMatrix< float, 3, 3 > MatrixFloat;
  bool operator!=( const MatrixFloat& m1, const MatrixFloat& m2 )
  { return ! ( m1 == m2 ); }
  namespace functors {
    bool operator==( Identity f1, Identity f2 )
    { return true; }
  }
}

/**
* Structure to store a traversal in a graph. Useful if you have a
* translation invariant graph and if you wish to repeat the traversal
* from another point.
*/
template <typename TVisitor>
struct TraversalReplay
{
  typedef TVisitor                  Visitor;
  typedef typename Visitor::Node    Node;
  typedef typename Visitor::Scalar  Scalar;
  typedef typename Visitor::Vertex  Vertex;
  typedef std::vector<Node>         Container;
  typedef typename Container::const_iterator ConstIterator;

  std::vector<Node> myNodes;
  
  TraversalReplay() {}

  void init( Visitor& visitor, Scalar dmax )
  {
    myNodes.clear();
    Node node;
    while ( ! visitor.finished() )
      {
        node = visitor.current();
        myNodes.push_back( node );
        if ( node.second > dmax ) break;
        visitor.expand();
      }
  }

  struct NodeLessComparator {
    bool operator()( const Node& n1, const Node& n2 ) const
    {
      return n1.second < n2.second;
    }
  };

  ConstIterator begin() const { return myNodes.begin(); }
  ConstIterator end() const   { return myNodes.end(); }
  // Returns an iterator pointing to the first element which does not compare less than val.
  ConstIterator find( Scalar val ) const 
  {
    Node dummy( Vertex(), val ); 
    return std::lower_bound( begin(), end(), dummy, NodeLessComparator() );
  }

};



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
template <typename TImageFct>
class DistanceToMeasure {
public:
  typedef TImageFct                  ImageFct;
  typedef typename ImageFct::Value   Value;
  typedef typename ImageFct::Point   Point;
  typedef typename ImageFct::Domain  Domain;
  typedef typename Domain::Space     Space;
  typedef typename Space::Vector     Vector;
  typedef typename Space::RealVector RealVector;

  typedef ExactPredicateLpSeparableMetric<Space,2> Distance;
  typedef DistanceToPointFunctor<Distance>         DistanceToPoint;
  typedef MetricAdjacency<Space, 1>                Graph;
  typedef DistanceBreadthFirstVisitor< Graph, DistanceToPoint, std::set<Point> > 
  DistanceVisitor;
  typedef TraversalReplay< DistanceVisitor >       DistanceTraversal;

public:
  DistanceToMeasure( Value m0, const ImageFct& measure, Value rmax = 10.0 )
    : myMass( m0 ), myMeasure( measure ), myDistance2( myMeasure.domain() ),
      myRMax( rmax )
  {
    init( myMeasure );
  }
  
  void init( const ImageFct& measure )
  {
    // Precompute traversal
    myP0 = *( myDistance2.domain().begin() );
    Value             m  = NumberTraits<Value>::ZERO;
    Value             d2 = NumberTraits<Value>::ZERO;
    Graph             graph;
    DistanceToPoint   d2pfct( Distance(), myP0 );
    DistanceVisitor   visitor( graph, d2pfct, myP0 );
    myTraversal.init( visitor, myRMax );

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

    inline const Domain& domain() const
  {
    return myMeasure.domain();
  }

  inline const ImageFct& measure() const
  {
    return myMeasure;
  }
  
  /// Distance function
  inline Value operator()( const Point& p ) const
  {
    return distance( p );
  }

  /// Gradient of distance^2 function
  inline RealVector normal( const Point& p ) const
  {
    return projection( p );
  }
  
  Value distance( const Point& p ) const
  {
    return sqrt( distance2( p ) );
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
    typedef DGtal::MetricAdjacency<Space, 1> Adjacency;
    std::vector<Point> neighborsP;
    std::back_insert_iterator<std::vector<Point> > outIterator(neighborsP);
    Adjacency::writeNeighbors(outIterator, p);

    typedef typename std::vector<Point>::iterator Iterator;
    Value distance_center = distance2( p );
    RealVector vectorToReturn;
    for (Iterator it = neighborsP.begin(), ite = neighborsP.end();
         it != ite; ++it) {
      Value distance = (myDistance2.domain().isInside(*it)) ? distance2( *it ) : distance_center;
      for (int d = 0; d < Point::dimension; d++) {
        if (p[d] < (*it)[d]) {
          Point otherPoint = *it;
          otherPoint[d] = p[d] + (p[d] - (*it)[d]);
          Value otherDistance = (myDistance2.domain().isInside(otherPoint)) ? distance2( otherPoint ) : distance_center;
          vectorToReturn[d] = ( abs( distance - distance_center) >= abs( distance_center - otherDistance) ) ? -(distance - distance_center) / 2.0 : -(distance_center - otherDistance) / 2.0;
        }		  		  		  
      }
    }
    return vectorToReturn;
    // Point p_left = box( p - Point( 1, 0 ) );
    // Point p_right = box( p + Point( 1, 0 ) );
    // Point p_down = box( p - Point( 0, 1 ) );
    // Point p_up = box( p + Point( 0, 1 ) );
    // Point p_front = box( p + Point
    // Value d2_center = distance2( p );
    // Value d2_left = distance2( p_left );
    // Value d2_right = distance2( p_right );
    // Value d2_down = distance2( p_down );
    // Value d2_up = distance2( p_up );
    // // Value gx = 
    // //   // std::min( ( d2_right - d2_left ) / ( p_right[ 0 ] - p_left[ 0 ] ),
    // //   std::min( ( d2_right - d2_center ) / ( p_right[ 0 ] - p[ 0 ] ),
    // //             ( d2_center - d2_left ) / ( p[ 0 ] - p_left[ 0 ] ) ); //  );
    // // Value gy = 
    // //   // std::min( ( d2_up - d2_down ) / ( p_up[ 1 ] - p_down[ 1 ] ),
    // //   std::min( ( d2_up - d2_center ) / ( p_up[ 1 ] - p[ 1 ] ),
    // //             ( d2_center - d2_down ) / ( p[ 1 ] - p_down[ 1 ] ) ); // );
    // bool right = abs( d2_right - d2_center ) >= abs( d2_center - d2_left );
    // bool up    = abs( d2_up    - d2_center ) >= abs( d2_center - d2_down );
    // Value gx = right ? ( d2_right - d2_center ) : ( d2_center - d2_left );
    // Value gy = up    ? ( d2_up    - d2_center ) : ( d2_center - d2_down );
    // return RealVector( -gx / 2.0, -gy / 2.0 );
    // // Value gx = (distance2( px2 ) - distance2( px1 ))
    // //   / ( 2.0 * ( px2[ 0 ] - px1[ 0 ] ) );
    // // Value gy = (distance2( py2 ) - distance2( py1 ))
    // //   / ( 2.0 * ( py2[ 1 ] - py1[ 1 ] ) );
    // // return RealVector( -gx, -gy );
  }

  Value computeDistance2( const Point& p )
  {
    typedef typename DistanceTraversal::Node          Node;
    typedef typename DistanceTraversal::ConstIterator ConstIterator;
    Value           last = NumberTraits<Value>::ZERO;
    Value             m  = NumberTraits<Value>::ZERO;
    Value             d2 = NumberTraits<Value>::ZERO;
    Vector shift = p - myP0;
    for ( ConstIterator it = myTraversal.begin(), itE = myTraversal.end(); it != itE; ++it )
      {
        const Node& node = *it;
        if ( ( node.second != last ) // all the vertices of the same layer have been processed. 
             && ( m >= myMass ) ) break;
        if ( node.second > myRMax ) { break; } // { d2 = m * myRMax; break; }
        Point q = node.first + shift;
        if ( myMeasure.domain().isInside( q ) )
        {
          Value mpt  = myMeasure( q );
          d2        += mpt * node.second * node.second; 
          m         += mpt;
          last       = node.second;
        }
      }
    return m > 0 ? d2 / m : myRMax * myRMax;
  }
  
  Value computeDistance2Naive( const Point& p )
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
    // trace.info() << p << endl;
    while ( ! visitor.finished() )
    {
      node = visitor.current();
      if ( ( node.second != last ) // all the vertices of the same layer have been processed. 
           && ( m >= myMass ) ) break;
      if ( node.second > myRMax ) { break; }
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
    return m > 0 ? d2 / m : myRMax * myRMax;
  }

public:
  Value myMass;
  const ImageFct& myMeasure;
  ImageFct myDistance2;
  Value myRMax;
  DistanceTraversal myTraversal;
  Point             myP0;
};


/**
 * A distance-like function d is a proper function (of infinite limit)
 * with the 1-semi concave property. The distance to a measure is a
 * typical example of such function.
 */
template <typename TDistanceLikeFunction>
struct DeltaVCM {
  typedef TDistanceLikeFunction                   DistanceLikeFunction;
  typedef typename DistanceLikeFunction::Value    Value;
  typedef typename DistanceLikeFunction::Point    Point;
  typedef typename DistanceLikeFunction::Domain   Domain;
  typedef typename DistanceLikeFunction::ImageFct ImageFct;
  typedef typename Domain::Space                  Space;
  typedef typename Space::Integer                 Integer;
  typedef typename Space::RealVector              RealVector;
  typedef typename DistanceLikeFunction::Value    Scalar;
  typedef DGtal::SimpleMatrix< Scalar,
                               Space::dimension,
                               Space::dimension > Matrix; ///< the type for nxn matrix of real numbers.
  // typedef typename Matrix::RowVector Vector;   ///< the type for N-vector of real numbers

  typedef ImageContainerBySTLVector<Domain,Matrix> MatrixField;

  DeltaVCM( const DistanceLikeFunction& delta, float R, float r )
    : myDelta( delta ), myR( R ), myr( r ),
      myVCM( delta.domain() ),
      myProjectedMeasure( delta.domain() )
  {
    init();
  }

  void init()
  {
    Matrix M;
    for ( typename Domain::ConstIterator it = myDelta.domain().begin(),
            itE = myDelta.domain().end(); it != itE; ++it )
      {
        Point      p = *it;
        // eliminates points too far away.
        if ( myDelta( p ) > myR ) continue;
        RealVector n = myDelta.projection( p );
        Point      q = Point( (Integer) round( p[ 0 ] + n[ 0 ] ),
                              (Integer) round( p[ 1 ] + n[ 1 ] ),
                              (Integer) round( p[ 2 ] + n[ 2 ] ) );
        // eliminates projections going outside the domain.
        if ( q != myDelta.box( q ) ) continue;
        for ( Dimension i = 0; i < Space::dimension; ++i ) 
          for ( Dimension j = 0; j < Space::dimension; ++j )
            M.setComponent( i, j, n[ i ] * n[ j ] ); 
        myVCM.setValue( q, myVCM( q ) + M ); // add tensor n x n
        myProjectedMeasure.setValue( q, myProjectedMeasure( q ) + myDelta.measure()( p ) );
      }
  }

  inline const Domain& domain() const
  {
    return myDelta.domain();
  }

  /**
     Computes the Voronoi Covariance Measure of the function \a chi_r.
     
     @tparam Point2ScalarFunction the type of a functor
     Point->Scalar. For instance functors::HatPointFunction and
     functors::BallConstantPointFunction are models of this type.
     
     @param chi_r the kernel function whose support is included in
     the cube centered on the origin with edge size 2r.
     
     @param p the point where the kernel function is moved. It must lie within domain.
  */
  template <typename Point2ScalarFunction>
  Matrix measure( Point2ScalarFunction chi_r, Point p ) const
  {
    Integer r = (Integer) ceil( myr );
    Point low = domain().lowerBound().sup( p - Point::diagonal( r ) );
    Point up = domain().upperBound().inf( p + Point::diagonal( r ) );
    //trace.info() << "r=" << r << " low=" << low << " up=" << up << std::endl;
    Domain local( low, up );
    Scalar mass = 0.0;
    Matrix M;
    for ( typename Domain::ConstIterator it = local.begin(), itE = local.end();
          it != itE; ++it )
      {
        Point q = *it;
        Scalar chi = chi_r( q - p );
        if ( chi <= 0.0 ) continue;
        // JOL: to check : I don't know if you should weight chi by the measure.
        // (0) no correction
        chi *= myDelta.measure()( q );     // (1) more stable than (2) and (0)
        // chi *= myProjectedMeasure( q ); // (2)
        //trace.info() << "chi=" << chi << " VCM=" << myVCM( q ) << endl;
        M += ::operator*(chi, myVCM( q ) ); // workaround simplematrix bug in DGtal.
      }
    return M;
  }

  // chi_r is normalized to have mass 1.
  template <typename Point2ScalarFunction>
  Matrix measure1( Point2ScalarFunction chi_r, Point p ) const
  {
    Integer r = (Integer) ceil( myr );
    Point low = domain().lowerBound().sup( p - Point::diagonal( r ) );
    Point up = domain().upperBound().inf( p + Point::diagonal( r ) );
    //trace.info() << "r=" << r << " low=" << low << " up=" << up << std::endl;
    Domain local( low, up );
    Scalar mass = 0.0;
    Matrix M;
    for ( typename Domain::ConstIterator it = local.begin(), itE = local.end();
          it != itE; ++it )
      {
        Point q = *it;
        Scalar chi = chi_r( q - p );
        if ( chi <= 0.0 ) continue;
        mass += chi;
        // JOL: to check : I don't know if you should weight chi by the measure.
        chi *= myDelta.measure()( q ); 
        //trace.info() << "chi=" << chi << " VCM=" << myVCM( q ) << endl;
        M += myVCM( q ) * chi; // else ::operator*(chi, myVCM( q )); workaround simplematrix bug in DGtal.
      }
    return mass > 0.0 ? M / mass : M;
  }

  
  DistanceLikeFunction myDelta;
  Scalar               myR;
  Scalar               myr;
  MatrixField          myVCM;
  ImageFct             myProjectedMeasure;
};




int main( int argc, char** argv )
{
  QApplication application(argc,argv);

  using namespace DGtal;
  using namespace DGtal::Z3i;
  
  typedef ImageContainerBySTLVector<Domain,unsigned char> GrayLevelImage3D;
  typedef ImageContainerBySTLVector<Domain,float>         FloatImage3D;
  typedef DistanceToMeasure<FloatImage3D>                 Distance;
  if ( argc <= 3 ) return 1;
  GrayLevelImage3D img  = GenericReader<GrayLevelImage3D>::import( argv[ 1 ] );
  float            mass = atof( argv[ 2 ] );
  float            rmax = atof( argv[ 3 ] );
  float            R    = atof( argv[ 4 ] );
  float            r    = atof( argv[ 5 ] );
  float            T1    = atof( argv[ 6 ] );
  float            T2    = atof( argv[ 7 ] );
  unsigned char    seuil = atof( argv[ 8 ] );

  {
    Viewer3D<> viewer;
    viewer.show();
    
    for ( Domain::ConstIterator it = img.domain().begin(), itE = img.domain().end();
          it != itE; ++it )
      {
        // std::cout << *it << " " << (int) img( *it ) << endl;
        if ( img( *it ) > seuil ) viewer << *it;
      }
    viewer << Viewer3D<>::updateDisplay;
    application.exec();
  }

  FloatImage3D     fimg( img.domain() );
  FloatImage3D::Iterator outIt = fimg.begin();
  for ( GrayLevelImage3D::ConstIterator it = img.begin(), itE = img.end();
        it != itE; ++it )
    {
      float v = 2.0 * ((float)*it) / seuil; // 255.0;
      v = std::min( 255.0f, v );
      *outIt++ = v;
    }
  trace.beginBlock( "Computing delta-distance." );
  Distance     delta( mass, fimg, rmax );
  const FloatImage3D& d2 = delta.myDistance2;
  trace.endBlock();

  float m = 0.0f;
  for ( typename Domain::ConstIterator it = d2.domain().begin(),
          itE = d2.domain().end(); it != itE; ++it )
    {
      Point p = *it;
      float v = sqrt( d2( p ) );
      m = std::max( v, m );
    }
  
  // GradientColorMap<float> cmap_grad( 0, m );
  // cmap_grad.addColor( Color( 255, 255, 255 ) );
  // cmap_grad.addColor( Color( 255, 255, 0 ) );
  // cmap_grad.addColor( Color( 255, 0, 0 ) );
  // cmap_grad.addColor( Color( 0, 255, 0 ) );
  // cmap_grad.addColor( Color( 0,   0, 255 ) );
  // cmap_grad.addColor( Color( 0,   0, 0 ) );
  // QApplication application(argc,argv);
  // Viewer3D<> viewer;
  // viewer.show();

  // for ( typename Domain::ConstIterator it = d2.domain().begin(),
  //         itE = d2.domain().end(); it != itE; ++it )
  //   {
  //     Point p = *it;
  //     float v = sqrt( d2( p ) );
  //     v = std::min( (float)m, std::max( v, 0.0f ) ); 
  //     viewer << CustomColors3D(Color(cmap_grad(v).red(), cmap_grad(v).green(), cmap_grad(v).blue(), 120), Color(cmap_grad(v).red(), cmap_grad(v).green(), cmap_grad(v).blue(),120) )
  //           << p;

  //     RealVector grad = delta.projection( p );
  //     viewer.addLine( p, p+grad );
  //   }
  // std::cout << endl;
  // viewer << Viewer3D<>::updateDisplay;
  // application.exec();

  trace.beginBlock( "Computing delta-VCM." );
  typedef DeltaVCM< Distance > DVCM;
  typedef DVCM::Matrix                     Matrix;
  DVCM dvcm( delta, R, r );
  trace.endBlock();


  typedef EigenDecomposition<3,float> LinearAlgebraTool;
  typedef functors::HatPointFunction<Point,float> KernelFunction;
  KernelFunction chi( 1.0, r );

  // Flat zones are metallic blue, slightly curved zones are white,
  // more curved zones are yellow till red.
  double size = 1.0;
  GradientColorMap<double> colormap( 0.0, T2 );
  colormap.addColor( Color( 128, 128, 255 ) );
  colormap.addColor( Color( 255, 255, 255 ) );
  colormap.addColor( Color( 255, 255, 0 ) );
  colormap.addColor( Color( 255, 0, 0 ) );
  Matrix vcm_r, evec, null;
  typedef PointVector<3,float> RealVector3f;
  RealVector3f eval;

  Viewer3D<> viewer;
  viewer.show();

  for ( Domain::ConstIterator it = dvcm.domain().begin(), itE = dvcm.domain().end();
        it != itE; ++it )
    {
      // Compute VCM and diagonalize it.
      Point p = *it;
      vcm_r = dvcm.measure( chi, p );
      if ( vcm_r == null ) continue;
      LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
      //double feature = eval[ 0 ] / ( eval[ 0 ] +  eval[ 1 ] );
      eval[ 0 ] = std::max( eval[ 0 ], 0.00001f );
      float tubular = ( eval[ 2 ] <= 0.00001f ) // (R*R/4.0) )
        ? 0
        : ( ( eval[ 1 ] + eval[ 2 ] ) / ( eval[ 0 ] + eval[ 1 ] + eval[ 2 ] ) );
      float bound = T1;
      float tubular2 = tubular * (eval[ 0 ] + eval[ 1 ] + eval[ 2 ] ) / (R*R*r*r*3.14f/12.0f);
      float display = tubular2 <= bound ? 0.0f : ( tubular2 - bound ) / (1.0f - bound);
      //: eval[ 1 ] / ( 1.0 + eval[ 0 ] ) / ( 1.0 + delta( p )*delta( p ) );
      //: eval[ 1 ] * eval[ 1 ] / ( 1.0 + eval[ 0 ] ) / ( 1.0 + delta( p ) );
      if (display > 0.01f)
        trace.info() << "l0=" << eval[ 0 ] << " l1=" << eval[ 1 ]
                     << " tub=" << tubular
                     << " tub2=" << tubular2
                     << " disp=" << display << std::endl;
      if (display > 0.5f*T2 )
        {
          viewer << CustomColors3D( Color::Black,
                                    colormap( display > T2 ? T2 : display ) )
                 << p;
          RealVector normal = evec.column( 0 );
          RealPoint rp( p[ 0 ], p[ 1 ] ); 
          viewer.addLine( rp - size*normal, rp + size*normal );
        }
    }      
  viewer << Viewer3D<>::updateDisplay;
  application.exec();
  return 0;
}
		 
