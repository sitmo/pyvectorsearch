/*! \brief PkTree, a spatial index and search structure.
*  \author M.A. van den Berg, thijs@sitmo.com
*  \version 8.0
*  \date  2004-2023
*
*  Copyright 2004-2023 M.A. van den Berg, thijs@sitmo.com
*
*  Licensed to the Apache Software Foundation (ASF) under one or more contributor license agreements; and to You under
*  the Apache License, Version 2.0. "
*/

#ifndef _PKSPATIAL_
#define _PKSPATIAL_

#include <cassert>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>

namespace pkspatial {


/*!
* Integer data type to represent the level of a node inside the tree.
* The root node has level zero.
*/
typedef short int level_type;


/*!
* Integer data type to represent the number of elements (points or child nodes) attached to a node.
*/
typedef unsigned  int count_type;


/*!
* The pktree used for creating spatial indices and spatial searches
*/
template<typename TIndexPoint, typename TQueryPoint = TIndexPoint>
class pktree
{
private:
       /*!
       * A simple node of the pktree.
       * The node is connected to it's parent, it's first child, and it's next sibling.
       * The node contains a pointer to coordinates, pc->[0] is the first coordinate
       * The number of children is dynamic.
       */
       template<typename TPoint>
       class pknode
       {
           friend class pktree<TIndexPoint, TQueryPoint>;

       public:
           typedef pknode<TPoint>  TNode;
           typedef pknode<TPoint>* TNodePtr;

       private:
           // Tree related variables
           TNodePtr    _parent;                    ///< Pointer to the parent node.
           TNodePtr    _next_sibling;              ///< Pointer to the next sibling.
           TNodePtr    _first_child;               ///< Pointer to the first child.
           count_type  _children;                  ///< Number of children in this node.

           TPoint      _point;             ///< Pointer to the coordinates

           /*!
           * Level of the node (distance from the root).
           * Nodes represent spatial hyperbox. The level of the node is the only information that is needed to
           * derive the spatial extend of the node, i.e. the length of the hyperbox along various directions.
           * A level of std::numeric_limits<level_type>::max() is the maximum possible distance to the root and
           * represents the smallest possible hyperbox. This value is used to represent point.
           */
           level_type  _level;

       public:
           /*!
           * Constructor.
           */
           pknode() :  _parent(0),
                       _next_sibling(0),
                       _first_child(0),
                       _children(0),
                       _level( std::numeric_limits<level_type>::max() )
           {
               //std::cout << "node created\n";
           }

           /*!
           * Destructor.
           */
           ~pknode()
           {
               if (_next_sibling != 0) delete _next_sibling;
               if (_first_child != 0) delete _first_child;
               //std::cout << "node deleted\n";
           }

           /*!
           * Test is a node represent a region or a point in space.
           */
           inline bool is_point_node() const
           {
               return (_level == std::numeric_limits<level_type>::max());
           };

           /*!
           * Attaches a childnode to the current node.
           * \param child a pointer to the child node being attached.
           */
           void attach_child(const TNodePtr child)
           {
               assert(child != 0);
               child->_parent = this;
               child->_next_sibling  = _first_child;
               _first_child = child;
               ++_children;
           };

           /*!
           * Attaches a sibling to the current node.
           * \param sibling a pointer to the sibling node being attached.
           */
           void attach_sibling(const TNodePtr sibling)
           {
               sibling->_next_sibling = _next_sibling;
               sibling->_parent = _parent;
               _next_sibling = sibling;
               _parent->_children++;
           };


           /*!
           * Detach this node from it's parent.
           */
           void detach()
           {
               if (_parent!=0)
               {
                   if (this == _parent->_first_child)
                       _parent->_first_child = _parent->_first_child->_next_sibling;
                   else
                   {
                       // find the left sibling
                       TNodePtr left_child = _parent->_first_child;
                       while (this != left_child->_next_sibling)
                           left_child=left_child->_next_sibling;

                       left_child->_next_sibling = left_child->_next_sibling->_next_sibling;
                   }
                   // update parent childcount & ensure a valid center
                   --(_parent->_children);
                   // remove pointers to prevent delete cascade
                   _next_sibling = 0;
               } // parent!=0
           };


           /*!
           * Detach the first child of this node from the tree.
           * \return A pointer to the detached node or 0 when no child was present.
           */
           TNodePtr detach_first_child()
           {
               if (_first_child == 0) return 0;
               TNodePtr detached_child = _first_child;
               _first_child  = _first_child->_next_sibling;

               // update childcount & ensure a valid center
               --_children;

               return detached_child;
           };


           /*!
           * Detach the next sibling from the tree.
           * \return A pointer to the detached node or 0 when no sibling was present.
           */
           TNodePtr detach_next_sibling()
           {
               assert(_parent != 0);

               if (_next_sibling != 0)
               {
                   TNodePtr detached_node = _next_sibling;
                   _next_sibling = _next_sibling->_next_sibling;

                   // update parent childcount & ensure a valid center
                   --(_parent->_children);

                   return detached_node;
               }
               return 0;
           };


           /*!
           * Returns a pointer to the n-th child of a node.
           * \param n the index of the child
           * \param base the base used in counting n-th (1 based, 0 based)
           * \return A pointer to the n-th child node or 0 when no node was found.
           */
           TNodePtr nth_child( const count_type n,
                               const count_type base=1) const
           {
               assert(n >= base);
               assert(_children > 0);

               if (n == base) return _first_child;
               else
               {
                   TNodePtr child = _first_child;
                   count_type child_number = base;

                   while  ( (child!=0) && (child_number<n) )
                   {
                       child = child->_next_sibling;
                       ++child_number;
                   }
                   return child;
               }
           };


           /*!
           * Moves a set child nodes to a target node.
           *
           * \param target_node The target node where the child nodes are moved to as extra children.
           * \param nodes_to_skip The initial number of child nodes to skip before beginning the move, these child nodes will stay attached to this node.
           * \param nodes_to_move The (subsequent) number of child nodes moves. Remaing child nodes will stay attached to this node.
           */
           void move_child_nodes(  TNodePtr target_node,
                                   const count_type nodes_to_skip,
                                   const count_type nodes_to_move)
           {
               assert(target_node != 0);
               assert(nodes_to_skip < _children);
               assert(nodes_to_move > 0);
               assert(nodes_to_skip + nodes_to_move <= _children);

               if (nodes_to_skip>0)
               {
                   //find the last node to skip and move its siblings one by one
                   TNodePtr last_skipped_node = nth_child(nodes_to_skip,1);
                   for (count_type i=0; i<nodes_to_move; ++i)
                   {
                       TNodePtr node_to_move = last_skipped_node->detach_next_sibling();
                       if (node_to_move != 0) target_node->attach_child(node_to_move);
                       else return;
                   }
               }
               else
               {
                   //move the first n nodes
                   for (count_type i=0; i<nodes_to_move; ++i)
                   {
                       TNodePtr node_to_move = detach_first_child();
                       if (node_to_move != 0) target_node->attach_child(node_to_move);
                       else return;
                   }
               }
           }
       }; //end of: class pknode

private:
   /*!
    * Node data type.
    */
   typedef pknode<TIndexPoint> TNode;

   /*!
    * Pointer to a node data type.
    */
   typedef pknode<TIndexPoint >* TNodePtr;

public:


   /*!
    * List of pointers to spatial-points.
    * This data type is user to store the results of a range search.
    * The range search will return an unknown number of points.
    */
   typedef  std::list<TIndexPoint> TResultsetRange;

   /*!
    * An vector of multiple <point,distance> pairs.
    * This data type is user to store the results of a k-nearest neighbor (knn) search.
    * The knn search will return a vector of the k points closest to a given point.
    * The distance to that point is also given.
    */
   typedef  std::vector< std::pair<double, TIndexPoint> >  TResultsetKnn;  ///< a vector of point-distance pairs*/


private:
   TNodePtr        _root;          ///< The root node.
   count_type      _dim;           ///< Dimension of the spatial points.
   count_type      _rank;          ///< Minimum number of children (K) for a node.
   count_type      _pointcount;    ///< Number of points in the tree.
   count_type      _nodecount;     ///< Number of nodes in the tree.
   count_type      _rr;            ///< Round robin steps, number of dimensions that are sliced simultaneous per level.
   double          _ratio;         ///< Number of subdivisions per slice, (3 -> slice the cube into 3 sub-cubes).
   double          _width;         ///< Spatial bounds at level==0, the width of the tree.
   std::vector<double> _O;         ///< Origin of the root.

   int dim_div_rr;                 ///< Precalculated constant: floor(dim/rr).

private:
   /*!
    * integer modulus: 5 mod 3 = 2, -5 mod 3 = 1
    * \param a the dividend.
    * \param b the divisor.
    * \return modulus(a,b) a non-negative integer, the remainder of division of a by b.
    */
   inline int math_modulus(const int a, const int b) const
   {
       int rem = a%b;
       return (rem<0) ? rem+b : rem;
   }


   /*!
    * integer division floor: 4/2=2 5/2=2 -4/2=-2 -3/2=-2
    * \param a the dividend.
    * \param b the divisor.
    * \return floor(a,b) the result of a/b rounded down to an integer value.
    */
   inline int math_int_div_floor(const int a, const int b) const
   {
       return a/b - (a%b<0 ? 1:0);
   }


   /*!
    *  integer division ceil: 4/2=2 5/2=3 -4/2=-2 -3/2=-1
    * \param a the dividend.
    * \param b the divisor.
    * \return ciel(a,b) the result of a/b rounded up to an integer value.
    */
   inline int math_int_div_ciel(const int a, const int b) const
   {
       return a/b + (a%b>0 ? 1:0);
   }


   /*!
    * The number of slices along an axis(c) at a specific level(l).
    * \param l the level.
    * \param c the axis number index (zero based), a 3-d point has three axis with indices 0,1 and 2.
    * \return the number of slices.
    */
   inline int geom_slices_at_level(const int l, const int c) const
   {
       if (l==std::numeric_limits<level_type>::max())
           return std::numeric_limits<level_type>::max();

       return math_int_div_floor(l*_rr+_dim-c-1,_dim);
   }


   /*!
    * The width of the node along an axis(c) at a specific level(l).
    * \param l the level.
    * \param c the axis number index (zero based), a 3-d point has three axis with indices 0,1 and 2.
    * \return the length of the node along the axis. Zero for a point node.
    */
   // --------------------------------------------------
   inline double geom_width_at_level(const int l, const int c) const
   {
       if (l>=std::numeric_limits<level_type>::max()-1)
           return 0.0;

       return _width * std::pow(_ratio,  -math_int_div_floor(l*_rr+_dim-c-1,_dim));
   }


   /*!
    * Level closest to the root with an axis (c) sliced (s) a nr of times.
    * \param c the axis number index (zero based), a 3-d point has three axis with indices 0,1 and 2.
    * \param s number of slices.
    * \return the level number.
    */
   inline int geom_mindepth_with_slicecount(const int s, const int c) const
   {
       if (s==std::numeric_limits<level_type>::max())
           return std::numeric_limits<level_type>::max();

       return math_int_div_ciel((s-1)*_dim+1+c,_rr);
   }

   /*!
    * Level furthest from the root with an axis (c) sliced (s) a nr of times.
    * \param c the axis number index (zero based), a 3-d point has three axis with indices 0,1 and 2.
    * \param s number of slices.
    * \return the level number.
    */
   inline int geom_maxdepth_with_slicecount(const int s, const int c) const
   {
       if (s==std::numeric_limits<level_type>::max())
           return std::numeric_limits<level_type>::max();

       return math_int_div_ciel(s*_dim+1+c,_rr)-1;
   }

   /*!
    * Find the level-number of the smallest nodes that will contain both nodes.
    * \param p1 a node.
    * \param p2 another node.
    * \param min_level the minimal level size that can be assumed, needs to be a priori valid information.
    * \param max_level the maximum level size that can be assumed, needs to be a priori valid information.
    * \return value can be min_level, max_level or any of the levels in between.
    */
   level_type mutual_level(        const TIndexPoint& p1,
                                   const TIndexPoint& p2,
                                   const level_type min_level,
                                   const level_type max_level) const
   {
       if (max_level - min_level <= dim_div_rr)
       {

           // Levels are close together: only some coordinates are divided once...
           // Test if p1 and p2 fit in a node at level l+1 (assuming that they are
           // together a level l).
           for (level_type l=min_level; l<max_level; ++l)
               if ( !is_point_covered_by_node<TIndexPoint>( l, p1, p2, l+1) ) return l;

           return max_level;

       }
       else
       {

           // There are quite a few levels between min_level and max_level
           // some coordinates are divided multiple times during the descent...

           // First we find the coordinate with the biggest distance between the points
           double max_width = std::fabs( p1[0] - p2[0] );
           count_type max_coord = 0;

           for (count_type i=1; i<_dim; ++i)
           {
               double point_dist = std::fabs(p1[i] - p2[i]);

               if (point_dist > max_width)
               {
                   max_width = point_dist;
                   max_coord = i;
               }
           }

           // Calculate the deepest level that is wider that max_width
           level_type level = max_level; //we assume very close
           int divisions = static_cast<int>(floor(std::log(_width / max_width)/ std::log(_ratio) ));
           level = std::min<level_type>(level, geom_maxdepth_with_slicecount(divisions,max_coord));

           if (level<=min_level)
               return min_level; //point are far apart

           if (level>max_level)
               level = max_level; //points are very close

           //make grid at level,
           double width = geom_width_at_level(level, 0);
           count_type last_divided_coord = math_modulus(level*_rr-1, _dim);


           //examine all coordinates
           for (count_type i=0; i<_dim; ++i)
           {
               int gridcell_1 = static_cast<int>(floor( (p1[i]-_O[i]) / width));
               int gridcell_2 = static_cast<int>(floor( (p2[i]-_O[i]) / width));

               // increase the node size (and grid size) until both point are in the same grid-cell
               if (gridcell_1 != gridcell_2)
               {
                   int slices = geom_slices_at_level(level, i);
                   int slice_lb = geom_slices_at_level(min_level, i);

                   while  (gridcell_1!=gridcell_2)
                   {
                       gridcell_1 = static_cast<int>(gridcell_1 / _ratio);
                       gridcell_2 = static_cast<int>(gridcell_2 / _ratio);
                       slices--;

                       if (slices < slice_lb)
                           return min_level;
                   }

                   level = geom_maxdepth_with_slicecount(slices,i);
                   width = geom_width_at_level(level, i);

                   if (level<=min_level)
                       return min_level;

                   last_divided_coord = math_modulus(level*_rr-1, _dim);
               }

               // the next coordinate is divided one time less
               if (i==last_divided_coord)
                   width *= _ratio;

           } // end of examining all coordinates

           return level;
       } // end of "if levels are far apart"
   };



   /*!
    * Fix the (center) location of a node so that all children are contained in it.
    * After resizing, a node might have the good size but wrong location.  This function fixes this efficiently.
    * \param n a node to fix.
    */
   inline void fix_internal_point(TNodePtr n)
   {
       if (n->_first_child!=0)
           n->_point = n->_first_child->_point;
   }


   /*!
    * Test if a point is covered by a node.
    * We only need to examine the coordinates that are divided at level "mutual_level+1"
    * up to and including the coordinated that are divided at node_level.
    * \param mutual_level is the know level that contains both the point and the node. this is typically the level of
    * the node being tested.
    * \param p the point.
    * \param node_ip an internal point of the node, could be the center or a child.
    * \param node_level the level (size specifier) of the node.
    * \return true is the point is located inside the spatial node, false otherwise.
    */
   template <typename TPoint>
   bool is_point_covered_by_node(
               const level_type mutual_level,
               const TPoint& p,
               const TIndexPoint& node_ip,
               const level_type node_level) const
   {
       if ( node_level == std::numeric_limits<level_type>::max() )
           return false;
       // determine coordinate ranges
       // ----------------------------
       count_type first,last;
       last = math_modulus( node_level*_rr - 1, _dim);

       if(node_level - mutual_level <= dim_div_rr)
               first = math_modulus(mutual_level*_rr, _dim);
       else    first = math_modulus(node_level*_rr, _dim);

       // determine width of the coordinate "first" at level "
       // --------------------------------
       int divisions = 1 + math_int_div_floor( node_level*_rr-first-1, _dim);
       double width = _width * std::pow(_ratio, -divisions);

       // loop over coordinates, checking segment numbers
       // -----------------------------------------------
       for (count_type i=first; ; ++i)
       {
           if (i == _dim)
           {
               i = 0;
               width /= _ratio;
           }

           if (   floor(  (p[i]       - _O[i])/width )
               != floor(  (node_ip[i] - _O[i])/width )
               )
               return false;

           if (i==last) break;
       }

       return true;
   };

   // ------------------------------------------------------
   // PRE:  node->first_child is a newly inserted point
   // POST: if return true then the first _rank children of
   //       the node will fit in a smaller node (maybe more)
   // ------------------------------------------------------
   bool can_create_subnode_after_insert(
               TNodePtr n,
               const count_type inserted_cluster_size)
   {
       TIndexPoint& inserted_point = n->_first_child->_point;

       // get last node of cluster
       count_type cluster_size = inserted_cluster_size;
       TNodePtr last_node_of_cluster = n->nth_child(inserted_cluster_size, 1);
       TNodePtr left_child = last_node_of_cluster;

       // start analysing nodes
       for (count_type i=inserted_cluster_size; i<n->_children; ++i)
       {
           if (cluster_size + n->_children -i < _rank)
               return false;

           TNodePtr child = left_child->_next_sibling;
           assert(child!=0);

           // check if the node fits at the next level. if so reorder

           if ( is_point_covered_by_node(  n->_level, child->_point, inserted_point, n->_level+1) )
           {

               if (last_node_of_cluster == left_child)
               {
                   // It is already the node after the cluster
                   left_child = child;

               }
               else
               {
                   // detach node from chain of children, and and move it to the begin of the chain
                   TNodePtr c = left_child->detach_next_sibling();
                   fix_internal_point(left_child->_parent);

                   last_node_of_cluster->attach_sibling(c);
               }
               last_node_of_cluster = child;
               ++cluster_size;
           } else
               left_child = left_child->_next_sibling;

           if (cluster_size >= _rank)
               return true;
       }
       return false;
   };



   // ------------------------------------------------------
   // The first "cluster_size" children in "node" have been
   // added to this node. The cluster does not contain _rank
   // children, that's the reason it was move to a supernode
   // (its parent).
   // POST:    The node is shrunk to the smallest node containing
   //          the cluster & a total of at least _rank children.
   //          When shrunk, some nodes might have been moved to the parent
   // RET:     Returns the number of nodes moved to the parent
   // ------------------------------------------------------
   count_type fix_node(
               TNodePtr n,
               const count_type cluster_size,
               const level_type cluster_level)
   {
       assert(n != 0);
       assert(n->_children >= _rank);
       assert(cluster_size > 0);
       assert(cluster_level != std::numeric_limits<level_type>::max());

       count_type nodes_moved_to_parent = 0;

       // Reposition the node (without moving).
       // If we shrink, we want the cluster to stay inside
       //assert(n->_first_child->_point._point != 0);
       fix_internal_point(n);

       TNodePtr last_node_of_cluster = n->nth_child(cluster_size, 1);
       TNodePtr child = last_node_of_cluster->_next_sibling;

       // allocate memory to store distance information
       count_type nodes_to_analyze = n->_children - cluster_size;
       std::vector< std::pair<level_type, TNodePtr> > dist_info(nodes_to_analyze);
       // Compare the cluster with every other node
       for (count_type i=0; i<nodes_to_analyze; ++i)
       {

           dist_info[i].first =  mutual_level(
               last_node_of_cluster->_point,
               child->_point,
               n->_level,
               cluster_level
           );

           dist_info[i].second = child;
           child = child->_next_sibling;
       } // end comparing the cluster with every other node

       // sort on distance
       std::sort(dist_info.begin(), dist_info.end(), std::less<std::pair<level_type,TNodePtr> >() );

       // check if we can create a sub-node of al least _rank children...
       if ( dist_info[_rank-1 - cluster_size].first > n->_level)
       {
           // create a new sub-node
           TNodePtr new_sub_node = new TNode();
           _nodecount++;
           new_sub_node->_level =  dist_info[_rank-1 - cluster_size].first;
           new_sub_node->_point = n->_first_child->_point;

           // add all the cluster elements to the new sub-node
           n->move_child_nodes(new_sub_node, 0, cluster_size);

           //reset the current node
           n->_first_child = 0;
           n->_children = 0;

           // distribute the nodes from the sorted list among the current and new node
           for (count_type i=0; i<nodes_to_analyze; ++i)
           {
               if (dist_info[i].first >= new_sub_node->_level)
                   new_sub_node->attach_child(dist_info[i].second);
               else
                   n->attach_child(dist_info[i].second);
           } // end distributing nodes

           // don't forget to attach the new sub-node to the node....
           n->attach_child(new_sub_node);

           if (n->_children < _rank)
           {

               // the current node is depleted, move all child-nodes up, and remove the node
               if (n->_parent!=0)
               {

                   nodes_moved_to_parent = n->_children;
                   n->move_child_nodes(n->_parent, 0, n->_children);
                   assert(n->_children==0);
                   assert(n->_first_child==0);
                   n->detach();
                   delete n;
                   _nodecount--;
               } // end moving depleted elements to its parent

           } // end handling a depleted node

       } // end  creating  a sub-node
       return nodes_moved_to_parent;
   };


   // ------------------------------------------------------
   // ALG: Try to find prove a quickly as possible that the
   // node can't shrink.
   // ------------------------------------------------------
   void shrink_node(TNodePtr n)
   {

       // -------------------------------------------------------
       // a node should at least be bigger than its biggest child
       // -------------------------------------------------------
       level_type level = n->_first_child->_level - 1;
       TIndexPoint p = n->_first_child->_point;
       TNodePtr c = n->_first_child->_next_sibling;

       while (c != 0)
       {
           level = std::min<level_type>(level, c->_level-1);
           if (level <= n->_level) return;
           c = c->_next_sibling;
       }

       for (count_type i=0; i<_dim; ++i)
       { //all coordinated

           TNodePtr c = n->_first_child->_next_sibling;
           double grid_width = geom_width_at_level(level, i);

           while (c!=0)
           { // all child nodes

               // Move to a bigger grid if that's evident
               double dist = fabs( p[i] - c->_point[i] );
               if (dist > grid_width)
               {
                   int slices = static_cast<int>(floor(std::log(_width/dist) / std::log(_ratio) ));
                   level = math_int_div_ciel( slices * _dim + 1 + i, _rr )-1;
                   grid_width = geom_width_at_level(level, i);
               }

               // Discretize the nodes onto a grid
               unsigned int gridcell_1 = static_cast<unsigned int>(floor( (p[i]-_O[i]) / grid_width));
               unsigned int gridcell_2 = static_cast<unsigned int>(floor( ( c->_point[i] - _O[i]) / grid_width));

               // When we are not in the same level...
               if (gridcell_1!=gridcell_2)
               {
                   int slices = geom_slices_at_level(level, i);
                   int slice_lb = geom_slices_at_level(n->_level, i); // the lower bound for this coord

                   while  (gridcell_1!=gridcell_2)
                   {
                       gridcell_1 = static_cast<unsigned int>(gridcell_1 / _ratio);
                       gridcell_2 = static_cast<unsigned int>(gridcell_2 / _ratio);
                       slices--;
                       if (slices < slice_lb) return;
                   }

                   level_type l = geom_maxdepth_with_slicecount(slices, i);

                   if (l<level)
                   {
                       level=l;
                       grid_width = geom_width_at_level(level, i);
                   }
               }

               c = c->_next_sibling;
           } // end node loop
       } // end coordinate loop
       n->_level = level;
       n->_point = n->_first_child->_point;
   };


   // ------------------------------------------------------
   // PRE: the point is inside the tree, covered by the root
   // POST: returns the smallest node that covers the point
   // ------------------------------------------------------
   // TODO: order children of a node: nodes first, than points, exit when we find a point
   template <typename TPoint>
   TNodePtr find_smallest_covering_node(const TPoint& p) const {
       TNodePtr n = _root;
       TNodePtr guess = n->_first_child;
       while (guess!=0) {
           if (
               (guess->_level != std::numeric_limits<level_type>::max()) && //must be a node
               is_point_covered_by_node(n->_level, p, guess->_point, guess->_level) //most cover the point
               ) {
                   n = guess;
                   guess = n->_first_child;
           } // if
           else guess = guess->_next_sibling;
       }
       return n;
   }


   // ======================================================
   // Increase the size of the root until if covers a point
   // ======================================================
   void grow_root_to_cover_point(const TIndexPoint &point) {
       TIndexPoint p = point;  //make a copy

       // If the root is empty: just add the point
       // and initialize the lower-bound of the grid
       // -------------------------------------
       if (_root->_children == 0)
       {
           _root->_level = std::numeric_limits<level_type>::max();
           _root->_point = point;
           for (count_type i=0; i<_dim; ++i)
               _O[i] = point[i];
           return;
       }

       // If this is second point: initialize the grid to the min
       // of the two points, ensuring [..)
       if (_root->_children == 1)
       {
           for (count_type i=0; i<_dim; ++i)
           {
               double dist = std::abs(p[i] - _O[i]);
               if (dist > 0)
               {
                   int slices = static_cast<int>(floor(std::log(_width / dist)/ std::log(_ratio) ));
                   level_type l = geom_maxdepth_with_slicecount(slices, i);
                   if (l < _root->_level) _root->_level = l;
               }
               _O[i] = std::min<double>(point[i], _root->_point[i]);
           }
           return;
       }


       // Calculate the the lowest level (biggest size)
       level_type new_level = _root->_level;
       for (count_type i=0; i<_dim; ++i)
       {
           double root_width = geom_width_at_level(_root->_level,i);
           double dist = 0;

           if  (p[i] < _O[i])
               dist = _O[i] - p[i] + root_width;
           else if  (p[i] >= _O[i] + root_width)
                   dist = p[i] - _O[i];

           if (dist>0)
           {
               int slices = static_cast<int>(floor(std::log(_width / dist)/ std::log(_ratio) ));
               level_type l = geom_maxdepth_with_slicecount(slices, i);
               if (l <  new_level)
                   new_level = l;
           }
       }

       // If the point is outside a node  we need to calculate
       // the new level (increase the node width) , and possibly
       // shift the lower-bound of the grid to the left.
       // --------------------------------------------------
       //std::cout << "_root->_level = " << _root->_level << std::endl;
       //std::cout << "new_level = " << new_level << std::endl;
       assert(new_level < _root->_level); //if not: why are we doing this?

       if (new_level < _root->_level)
       {

           // move the origin of the grid
           for (count_type i=0; i<_dim; ++i)
               if (point[i] < _O[i])
                   _O[i] += geom_width_at_level(_root->_level, i)
                          - geom_width_at_level(new_level, i);

           // if the previous root was valid, make it a sub
           if (_root->_children >= _rank)
           {
               TNodePtr new_root = new pknode<TIndexPoint>();
               ++_nodecount;
               new_root->attach_child(_root);
               new_root->_level = new_level;
               new_root->_point = _root->_point;
               _root = new_root;
           }
           else
           {
               // If the previous root was not valid, just grow it
               _root->_level = new_level;
           }
       }
   }


public:
   void search_range(
               const TQueryPoint& p,
               const double range,
               TResultsetRange& result) const
   {
       recursive_find_range(p, range*range, result, _root);
   }

   //-------------------------------------------------------------------------------------------------
   //Title:     Optimal Multi-Step k-Nearest Neighbor Search
   //Author:    Thomas Seidl, Hans-Peter Kriegel
   //Published: Proc. ACM SIGMOD Int. Conf. on Management of Data, Seattle, Washington, June 1998
   //Source:    http://knight.cis.temple.edu/~vasilis/Courses/CIS750/Papers/Knearestneighbor-Seidl.pdf
   // result is a std::vector<std::pair<double, T>>. The fist value of the pair is the square distance.
   //-------------------------------------------------------------------------------------------------
   void search_knn(
               const TQueryPoint& p,
               const count_type k,
               TResultsetKnn& result,
               const bool sort=false) const
   {
        double bound;
        count_type bound_index;
        TNodePtr start_node  = find_smallest_covering_node<TQueryPoint>(p);

        result.clear();
        recursive_add_points(p, start_node, result, k, bound, bound_index, 0, true);

        // sort results on distance
        if (sort==true)
            std::sort(
                result.begin(),
                result.end(),
                [](
                const std::pair<double, TIndexPoint> &left,
                const std::pair<double, TIndexPoint> &right
                ) {
                    return left.first < right.first;
                }
            );
   }



private:
   void recursive_add_points(
               const           TQueryPoint& p,
               const           TNodePtr n,
               TResultsetKnn&  result,
               const count_type k,
               double          &bound,
               count_type&     bound_index,
               const TNodePtr  node_to_skip,
               const bool      traverse_up) const
   {
       for (TNodePtr c = n->_first_child; c!=0; c=c->_next_sibling)
       {
           if (c != node_to_skip)
           {
               if (c->is_point_node())
               {
                   // always calculate distance
                   double d = 0;
                   for (count_type i=0; i<_dim; ++i)
                       d += std::pow(p[i] - c->_point[i], 2);

                   // two options: we have less that K, or exact K neighbors
                   if (result.size() < k)
                   {

                       // Less than K nn so far: always add points till we reach K nn's
                       result.push_back( std::pair<double,TIndexPoint>(d,c->_point) );

                       // initialize or update bound, bound keeps track of the maximum distance
                       if ( (result.size() == 1) || (d > bound) )
                       {
                           bound_index = result.size()-1;
                           bound = d;
                       }

                   }
                   else
                   {
                       //we have K NN,  replace the furthest neighbor when this is a better neighbor
                       if (d < bound)
                       {
                           //replace furthest neighbor
                           result[bound_index].first = d;
                           result[bound_index].second = c->_point;

                           //find new furthest neighbor
                           bound = d; //better that "bound = 0;"
                           for (count_type i = 0; i<result.size(); ++i)
                               if (result[i].first > bound)
                               {
                                   bound_index = i;
                                   bound = result[i].first;
                               }
                       }
                   } // end of the case of K matches

               }

               else
               {
                   //handle non-point nodes, i.e. nodes with sub-nodes

                   // we always need to go in the node & get points when we don't have enough matches
                   if (result.size() < k)
                       recursive_add_points(p, c, result, k, bound, bound_index, 0, false);
                   else
                   {
                       // we need to check the minimal distance between the point and this node, and so check
                       // if a potential better match exists in the node or not

                       double min_dist = 0.0;
                       TIndexPoint p2 = c->_point;

                       for (unsigned int i=0; i<_dim; ++i)
                       {
                           double w = geom_width_at_level(c->_level,i);    //grid width
                           double lb = _O[i] + std::floor((p2[i] -_O[i])/w)*w; //lower bound of the node grid-cell
                           if (p[i]<lb)
                               min_dist += std::pow(lb - p[i],2); // left of node
                           else if (p[i]>lb+w)
                               min_dist += std::pow(p[i] - lb - w,2); // right of node

                           if (min_dist >= bound) break; // no need for further refinement of distance calculation
                       }

                       if (min_dist < bound)
                           recursive_add_points(p, c, result, k, bound, bound_index, 0, false);
                   }
               } // if (c->_point.is_point_node() else ...
           } // if (c != node_to_skip)
       }  // end: for
       // after processing all children, go to the parent and process the remainder of the tree
       if ((traverse_up==true) && (n->_parent != 0))
           recursive_add_points(p, n->_parent, result, k, bound, bound_index, n, traverse_up);
   }


   /*!
    * Add all points in the node (and it's children) to the result.
    * \param n The node.
    * \param result The result.
    */
   void recursive_add_points(const TNodePtr n,  TResultsetRange& result) const
   {
       for (TNodePtr c = n->_first_child; c!=0; c=c->_next_sibling)
       {
           if (c->is_point_node())
               result.push_back(c->_point);
           else
               recursive_add_points(c, result);
       }
   }


   /*!
    * Find all points that are within range distance to a point.
    *
    * There are three different cases when looking at the distance between a node and a point:
    * \li The whole region occupied by the node is inside the range. In this case we select all point that
    *     reside in the node
    * \li The whole region occupied by the node is outside the range. We ignore the node and all it's children
    * \li Part of the node is withing the range. In this case we examine all children of the node.
    *
    * \param p Center point of the range.
    * \param range The range (Euclidian distance).
    * \param result The result set.
    * \param n The node (and its children) to examine for points that are within the given range.
    */
   void recursive_find_range(
               const TQueryPoint& p,
               const double range,
               TResultsetRange& result,
               const TNodePtr n) const
   {

       // the simple case: we have to examine a point-point distance.
       if (n->is_point_node())
       {
               double dist = 0.0;
               TIndexPoint p2 = n->_point;
               for(count_type i=0; i<_dim; ++i)
               {
                   dist += std::pow( p[i] - p2[i], 2);
                   if (dist>range) return; // we know enough, exit
               }
               result.push_back(n->_point);
       }

       // the more complicated case: we have to examine the point-node distance
       else
       {
           double min_dist = 0.0;      // the minimal point-node distance
           double max_dist = 0.0;      // the maximal point-node distance

           TIndexPoint p2 = n->_point;
           for (unsigned int i=0; i<_dim; ++i)
           {
               double w = geom_width_at_level(n->_level,i);
               double lb = _O[i] + std::floor((p2[i] -_O[i])/w)*w;

               if (p[i] < lb)
               {
                   // point is left of the node
                   min_dist += std::pow(lb - p[i],2);
                   max_dist += std::pow(lb - p[i] + w,2);
               }
               else if (p[i] > lb+w)
               {

                   // point is right of the node
                   min_dist += std::pow(p[i] - lb - w,2);
                   max_dist += std::pow(p[i] - lb,2);
               }
               else
               {

                   // point is inside node
                   double dx = std::max(p[i] - lb, lb + w - p[i]);
                   max_dist += std::pow(dx,2);
               }

               // we can stop calculating the min /max distance
               // per coordinate if the node is outside the range.
               if (min_dist > range) return;

           } //end for

           // there are two possible cases: then whole node is
           // within the range, or part of the node is inside the range.
           if (max_dist <= range)
               recursive_add_points(n, result);
           else
               for (TNodePtr c = n->_first_child; c!=0; c=c->_next_sibling)
                   recursive_find_range(p, range, result, c);

       } //end of: if
   } //end of: recursive_find_range

public:
   /*!
    * Default Constructor.
    */
   pktree() : _dim(0),_rank(3),_rr(2),_ratio(2),_width(1.0), _pointcount(0), _nodecount(0), _root(0)
   {
   };

   /*!
    * Constructor.
    * \param dim Dimension of the points.
    * \param rank Minimal number of child-nodes.
    * \param round_robin_step Number of dimensions to segment in one go.
    * \param dividing_ratio Number of segments when dividing.
    */
   pktree(
           int dim,
           int rank=3,
           int round_robin_step=2,
           int dividing_ratio=2
           ) :
           _dim(dim),
           _rank(rank),
           _rr(round_robin_step),
           _ratio(dividing_ratio),
           _width(1.0), _pointcount(0),
           _nodecount(1)
   {
       _root =  new TNode();
       _O.resize(_dim);
       for (unsigned int i=0; i<_dim; ++i)
           _O[i] = 0.0;
       dim_div_rr = math_int_div_floor(_dim , _rr);
   }

   ~pktree()
   {
       if (_root != 0) delete _root;
   }

   /*!
    * Initialize or resets the tree.
    * \param dim Dimension of the points.
    * \param rank Minimal number of childnodes.
    * \param round_robin_step Number of dimensions to segment in one go.
    * \param dividing_ratio Number of segments when dividing.
    */
   void init(
       const int dim,
       const int rank=3,
       const int round_robin_step=2,
       const int dividing_ratio=2
       )
   {
       //std::cout << "pktree::init()" << std::endl;

       _dim    = dim;
       _rank   = rank;
       _rr     = round_robin_step;
       _ratio  = dividing_ratio;
       _width  = 1.0;
       _pointcount = 0;
       _nodecount = 1;

       if (_root != 0) delete _root;
       _root =  new TNode();

       _O.clear();
       _O.resize(_dim);

       for (unsigned int i=0; i<_dim; ++i)
           _O[i] = 0.0;
       dim_div_rr = math_int_div_floor(_dim , _rr);

   }

   void clear()
   {
       init(_dim, _rank, _rr, _ratio);
   }


   /*!
    * Insert a point into the tree.
    * \param p The pointer to a point that will be added to the index.
    */
   void insert(const TIndexPoint& p)
   {
       ++_pointcount;

       // find the leaf-node that contains the point
       TNodePtr enclosing_node;

       if ( !is_point_covered_by_node<TIndexPoint>(
                   std::numeric_limits<level_type>::min(),
                   p,
                   _root->_point,
                   _root->_level)
           )
       {
           //std::cout << "p = "; for (int i=0; i<_dim; ++i) std::cout << " " << p[i]; std::cout << std::endl;
           grow_root_to_cover_point(p);
           enclosing_node = _root;
       }
       else
           enclosing_node = find_smallest_covering_node(p);

       // make a new node containing the point, and attach it to the leaf.
       TNodePtr node_new = new pknode<TIndexPoint>();
       node_new->_point = p;
       enclosing_node->attach_child(node_new);

       // Check if the leaf-node is crowded.
       // Create a sub-node with_rank children if possible
       if (can_create_subnode_after_insert(enclosing_node, 1))
       {

           // We want to move _rank children to a new sub-node,
           // effectively reducing the number of children with _rank-1.
           // If the node has enough children or if it's the root we move the
           // children to a sub-node.
           // If the node has only a few children, we shrink it and move
           // outliers to the parent. This saves us a new+delete action
           if (
               (enclosing_node->_children < 2*_rank-1) &&
               (enclosing_node->_parent != 0)
           )
           {
               // Move all but the first _rank children to the parent
               // We know that the remaining _rank children will fit in a node at
               // least one level deeper.
               count_type nodes_to_parent = enclosing_node->_children - _rank;
               level_type nodes_from_level = enclosing_node->_level;
               enclosing_node->move_child_nodes(enclosing_node->_parent, _rank, nodes_to_parent);
               enclosing_node->_level = enclosing_node->_level+1;
               shrink_node(enclosing_node);

               TNodePtr node_to_fix = enclosing_node->_parent;

               // All but the first _rank nodes have been moved to the parent
               if ((nodes_from_level - node_to_fix->_level > 1))
               {
                   if (can_create_subnode_after_insert(node_to_fix, nodes_to_parent))
                   {
                       // fix parent
                       while (
                               (nodes_to_parent>0) &&
                               (node_to_fix!=0) )
                       {
                           level_type level_before_fix = node_to_fix->_level;
                           TNodePtr next_node_to_fix = node_to_fix->_parent;

                           // fix node: nodes_to_parent>0 implies that is it removed
                           nodes_to_parent = fix_node(node_to_fix, nodes_to_parent, nodes_from_level);
                           nodes_from_level = level_before_fix;
                           node_to_fix = next_node_to_fix;
                       } // end: while there are nodes moved to the parent

                   } // end: if can_create_subnode_after_insert
               }// end: if more that one level between node and its parent

           }
           else
           {
               // Move _rank children from the current node to a new sub-node.
               // The current node stays valid
               TNodePtr sub_node = new TNode();
               ++_nodecount;

               sub_node->_point = enclosing_node->_first_child->_point;
               sub_node->_level = enclosing_node->_level + 1;

               enclosing_node->move_child_nodes(sub_node, 0, _rank);
               fix_internal_point(enclosing_node);

               shrink_node(sub_node);
               enclosing_node->attach_child(sub_node);
           }
       }
   } // end of: insert


}; // en of: pktree



template<typename TKey, typename TVal>
class pkmap_traits
{
public:

   typedef TKey key_type;
   typedef TVal value_type;

   // Key Value storage
   TKey    key;
   TVal    value;

   //Constructors
   pkmap_traits() : key(TKey()), value(TVal()) {};
   pkmap_traits(const TKey& Key, const TVal& Value) : key(Key),value(Value) {};

   // Trait operator
   double operator[] (const count_type& i) const { return key[i]; }

};


template<typename TKey,typename TVal>
class pkmap : public pktree<pkmap_traits<TKey,TVal>, TKey >
{

public:
   typedef TKey key_type;
   typedef TVal value_type;
   typedef typename pktree<pkmap_traits<TKey,TVal>, TKey >::TResultsetKnn TResultsetKnn;
private:
   typedef pkmap_traits<TKey,TVal> map_type;

public:

   pkmap() : pktree<pkmap_traits<TKey,TVal>, TKey >()
   {
   }

   pkmap(int dim, int rank=3, int round_robin_step=2, int dividing_ratio=2) :
     pktree<pkmap_traits<TKey,TVal>, TKey >(dim, rank, round_robin_step, dividing_ratio)
   {
   }


   /*!
    * Insert a key-value pair in the tree.
    */
   void insert(const TKey& k,const TVal& v)
   {
       pktree<pkmap_traits<TKey,TVal>, TKey >::insert( map_type(k,v) );
   }


}; // end of pkmap

} // end of namespace
#endif