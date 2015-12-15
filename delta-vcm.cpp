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
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/kernel/Point2ScalarFunctors.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

// work-arounds for DGtal
namespace DGtal {
  typedef SimpleMatrix< double, 2, 2 > MatrixDouble;
  bool operator!=( const MatrixDouble& m1, const MatrixDouble& m2 )
  { return ! ( m1 == m2 ); }
  typedef SimpleMatrix< float, 2, 2 > MatrixFloat;
  bool operator!=( const MatrixFloat& m1, const MatrixFloat& m2 )
  { return ! ( m1 == m2 ); }
  namespace functors {
    bool operator==( Identity f1, Identity f2 )
    { return true; }
  }
}

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

  DeltaVCM( const DistanceLikeFunction& delta, double R, double r )
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
                              (Integer) round( p[ 1 ] + n[ 1 ] ) );
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
  using namespace DGtal;
  using namespace DGtal::Z2i;
  
  typedef ImageContainerBySTLVector<Domain,unsigned char> GrayLevelImage2D;
  typedef ImageContainerBySTLVector<Domain,double>         DoubleImage2D;
  typedef DistanceToMeasure<DoubleImage2D>                 Distance;
  if ( argc <= 3 ) return 1;
  GrayLevelImage2D img  = GenericReader<GrayLevelImage2D>::import( argv[ 1 ] );
  double           mass = atof( argv[ 2 ] );
  double           rmax = atof( argv[ 3 ] );
  double           R    = atof( argv[ 4 ] );
  double           r    = atof( argv[ 5 ] );
  double           T1    = atof( argv[ 6 ] );
  double           T2    = atof( argv[ 7 ] );
  DoubleImage2D     fimg( img.domain() );
  DoubleImage2D::Iterator outIt = fimg.begin();
  for ( GrayLevelImage2D::ConstIterator it = img.begin(), itE = img.end();
        it != itE; ++it )
    {
      double v = ((double)*it) / 255.0;
      *outIt++ = v;
    }
  trace.beginBlock( "Computing delta-distance." );
  Distance     delta( mass, fimg, rmax );
  const DoubleImage2D& d2 = delta.myDistance2;
  trace.endBlock();

  double m = 0.0f;
  for ( typename Domain::ConstIterator it = d2.domain().begin(),
          itE = d2.domain().end(); it != itE; ++it )
    {
      Point p = *it;
      double v = sqrt( d2( p ) );
      m = std::max( v, m );
    }

  GradientColorMap<double> cmap_grad( 0, m );
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
      double v = sqrt( d2( p ) );
      v = std::min( (double) m, std::max( v, 0.0 ) ); 
      board << CustomStyle( p.className(),
                            new CustomColors( Color::Black, cmap_grad( v ) ) )
            << p;

      RealVector grad = delta.projection( p );
      // / ( 1.1 - ( (double)img( *it ) ) / 255.0 ) ;
      board.drawLine( p[ 0 ], p[ 1 ], p[ 0 ] + grad[ 0 ], p[ 1 ] + grad[ 1 ], 0 );
    }
  std::cout << endl;
  board.saveEPS("dvcm-delta2.eps");
  board.clear();
  
  trace.beginBlock( "Computing delta-VCM." );
  typedef DeltaVCM< Distance > DVCM;
  typedef DVCM::Matrix                     Matrix;
  DVCM dvcm( delta, R, r );
  trace.endBlock();

  {
    GrayLevelImage2D pm_img( dvcm.myProjectedMeasure.domain() );
    DoubleImage2D::ConstIterator it    = dvcm.myProjectedMeasure.begin();
    DoubleImage2D::ConstIterator itE   = dvcm.myProjectedMeasure.end();
    GrayLevelImage2D::Iterator  outIt = pm_img.begin();
    for ( ; it != itE; ++it )
      {
        double v = std::max( 0.0, std::min( (*it) * 255.0, 255.0 ) );
        *outIt++ = v;
      }
    
    GenericWriter< GrayLevelImage2D >::exportFile( "dvcm-projmeasure.pgm", pm_img );
  }

  typedef EigenDecomposition<2,double> LinearAlgebraTool;
  typedef functors::HatPointFunction<Point,double> KernelFunction;
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
  RealVector eval;
  for ( Domain::ConstIterator it = dvcm.domain().begin(), itE = dvcm.domain().end();
        it != itE; ++it )
    {
      // Compute VCM and diagonalize it.
      Point p = *it;
      vcm_r = dvcm.measure( chi, p );
      if ( vcm_r == null ) continue;
      LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
      //double feature = eval[ 0 ] / ( eval[ 0 ] +  eval[ 1 ] );
      eval[ 0 ] = std::max( eval[ 0 ], 0.00001 );
      double tubular = ( eval[ 1 ] <= 0.00001 ) // (R*R/4.0) )
        ? 0
        : ( eval[ 1 ] / ( eval[ 0 ] + eval[ 1 ] ) );
      double bound = T1;
      double tubular2 = tubular * (eval[ 0 ] + eval[ 1 ]) / (R*R*r/12.0);
      double display = tubular2 <= bound ? 0.0 : ( tubular2 - bound ) / (1.0 - bound);
      //: eval[ 1 ] / ( 1.0 + eval[ 0 ] ) / ( 1.0 + delta( p )*delta( p ) );
      //: eval[ 1 ] * eval[ 1 ] / ( 1.0 + eval[ 0 ] ) / ( 1.0 + delta( p ) );
      trace.info() << "l0=" << eval[ 0 ] << " l1=" << eval[ 1 ]
                   << " tub=" << tubular
                   << " tub2=" << tubular2
                   << " disp=" << display << std::endl;
      board << CustomStyle( p.className(), 
                            new CustomColors( Color::Black,
                                              colormap( display > T2 ? T2 : display ) ) )
            << p;
      // Display normal
      RealVector normal = evec.column( 0 );
      RealPoint rp( p[ 0 ], p[ 1 ] ); 
      Display2DFactory::draw( board, size*normal, rp );
      Display2DFactory::draw( board, -size*normal, rp );
    }      
  board.saveEPS("dvcm-hat-r.eps");
  board.clear();
  
  return 0;
}
