
/** This is the implementation of several types of extensional fuzzy numbers with application examples 
 *
 * @author Krzysztof Siminski
 * @date   2024-09-06
 * @date   2026-07-30
 *
 * CITATION REQUEST:
 * Please cite this if you use this code:
 * 
 *
 *
 *
 * */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <random>
#include <sstream>
#include <unordered_map>
#include <vector>

#define debug(x) std::cerr << __FILE__ << " (" << __LINE__ << ") " << #x << " : " << (x) << std::endl

///////////////////////////////////////


namespace  ksi 
{
   class tnorm 
   {
      public: 
         virtual ~tnorm() = default;
         virtual double value (const double a, const double b) const = 0;
   };

   class tnorm_product : public tnorm
   {
      public: 
         virtual double value (const double a, const double b) const override 
         {
            return a * b;
         }
   };

   class fuzzy_set
   {
      public :
         virtual ~fuzzy_set() = default;
         virtual double membership (const double x) const = 0;
   };

   class fuzzy_set_triangle : public fuzzy_set
   {
      double _supp_left, _core, _supp_right;
      public :
      fuzzy_set_triangle (const double support_left, const double core, const double support_right) : _supp_left(support_left), _core(core), _supp_right(support_right) {}

      virtual double membership (const double x) const override
      {
         if (x < _supp_left or x > _supp_right)
            return 0.0;
         if (x < _core)
         {
            return (x - _supp_left) / (_core - _supp_left);
         }
         else 
         {
            return (_supp_right - x) / (_supp_right - _core);
         }
         return 0.0;
      }
   };
   //////////////////////////

   std::map<double, double> fuzzy_operation (const fuzzy_set & A, const fuzzy_set & B, const tnorm & t, const auto op, const double mini_x, const double maxi_x, const double step)
   {
      std::map<double, double> result;

      for (double x = mini_x; x <= maxi_x; x += step)
      {
         for (double y = mini_x; y <= maxi_x; y += step)
         {
            double memb = t.value(A.membership(x), B.membership(y));
            double z = op(x, y);
            double prev = result[z];
            // debug(x); debug(y); debug(z); debug(value); debug(prev);
            result[z] = std::max(prev, memb);
         }
      }
      return result;
   }

   /////////////////////////
   template <class T> 
   void print (const T & a)
   {
      std::cout << a << " ";
   }
   template <class T> 
   void print (const std::vector<T> & A)
   {
      std::cout << "[ ";
      for (const auto & a : A)
         print(a);
      std::cout << "]" << std::endl;
   }
   /////////////////////////

   class extensional_fuzzy_number
   { 
      public:
         virtual ~extensional_fuzzy_number() = default;
   };

   class triangular_efn : public extensional_fuzzy_number
   {
      public:
         triangular_efn () : triangular_efn(0.0, 0.0) {} 
         triangular_efn (const double number, const double p) : _core(number), _p(p) {}
         static double relation (const double x, const double core, const double p)
         {
            return std::max (0.0, 1 - std::abs(x - core) / p);
         }
         static double equal (const triangular_efn & l, const triangular_efn & r)
         {
            double d = std::abs(l._core - r._core);
            double p = std::max(l._p, r._p); 
            return relation(0.0, d, p);
         }
         static double less (const triangular_efn & l, const triangular_efn & r)
         {
            if (l._core > r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }
         }
         static double greater (const triangular_efn & l, const triangular_efn & r)
         {
            if (l._core < r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }
         }
         static double less_equal (const triangular_efn & l, const triangular_efn & r)
         {
            if (l._core > r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }
         static double greater_equal (const triangular_efn & l, const triangular_efn & r)
         {
            if (l._core < r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }

         //////////////////
         static triangular_efn zero () 
         {
            return {0.0, 0.0};
         }
         static triangular_efn infinity () 
         {
            return {std::numeric_limits<double>::infinity(), 0.0};
         }
         //////////////////
         triangular_efn operator * (const double & r) const
         {
            return triangular_efn (this->_core * r, this->_p);
         }
         triangular_efn operator * (const triangular_efn & r) const
         {
            return triangular_efn (this->_core * r._core, std::max(this->_p, r._p));
         }
         triangular_efn operator + (const triangular_efn & r) const
         {
            return triangular_efn (this->_core + r._core, std::max(this->_p, r._p));
         }
         triangular_efn operator - (const triangular_efn & r) const
         {
            return triangular_efn (this->_core - r._core, std::max(this->_p, r._p));
         }
         triangular_efn operator / (const triangular_efn & r) const
         {
            return triangular_efn (this->_core / r._core, std::max(this->_p, r._p));
         }
         friend std::ostream & operator << (std::ostream & sos, const triangular_efn & t)
         {
            return sos << "(" << t._core << ", " << t._p << ")";
         }
         //////////////////
      protected:
         double _core; // core of the number
         double _p; // half of the support

   };

   class trapezoidal_efn : public extensional_fuzzy_number 
   {
      public:
         trapezoidal_efn () : trapezoidal_efn(0.0, 0.0) {} 
         trapezoidal_efn (const double number, const double p = 0) : _core(number), _p(p) {}
         static double relation (const double x, const double core, const double p, const double k)
         {
            return std::min (1.0, std::max (0.0, 1 - std::abs(x - core) / p) / (1 - k));
         }
         static double equal (const trapezoidal_efn & l, const trapezoidal_efn & r)
         {
            double d = std::abs(l._core - r._core);
            double p = std::max(l._p, r._p); 
            return relation(0.0, d, p, _k);
         }
         static double less (const trapezoidal_efn & l, const trapezoidal_efn & r)
         {
            if (l._core > r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }

         }
         static double greater (const trapezoidal_efn & l, const trapezoidal_efn & r)
         {
            if (l._core < r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }
         }
         static double less_equal (const trapezoidal_efn & l, const trapezoidal_efn & r)
         {
            if (l._core > r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }
         static double greater_equal (const trapezoidal_efn & l, const trapezoidal_efn & r)
         {
            if (l._core < r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }

         //////////////////
         static trapezoidal_efn zero () 
         {
            return {0.0, 0.0};
         }
         static trapezoidal_efn infinity () 
         {
            return {std::numeric_limits<double>::infinity(), 0.0};
         }
         //////////////////
         trapezoidal_efn operator * (const double & r) const
         {
            return trapezoidal_efn (this->_core * r, this->_p);
         }
         trapezoidal_efn operator * (const trapezoidal_efn & r) const
         {
            return trapezoidal_efn (this->_core * r._core, std::max(this->_p, r._p));
         }
         trapezoidal_efn operator + (const trapezoidal_efn & r) const
         {
            return trapezoidal_efn (this->_core + r._core, std::max(this->_p, r._p));
         }
         trapezoidal_efn operator - (const trapezoidal_efn & r) const
         {
            return trapezoidal_efn (this->_core - r._core, std::max(this->_p, r._p));
         }
         trapezoidal_efn operator / (const trapezoidal_efn & r) const
         {
            return trapezoidal_efn (this->_core / r._core, std::max(this->_p, r._p));
         }
         friend std::ostream & operator << (std::ostream & sos, const trapezoidal_efn & t)
         {
            return sos << "(" << t._core << ", " << t._p << ")";
         }
         //////////////////
      protected:
         double _core; // core of the number
         double _p; // half of the support
         static constexpr double _k = 0.5; // length of core in comparison to the length of the support, 0 <= k < 1 

   };

   class triangular_asymmetric_efn  : public extensional_fuzzy_number
   {
      public:
         triangular_asymmetric_efn () : triangular_asymmetric_efn(0.0, 0.0, 0.0) {} 
         triangular_asymmetric_efn (const double number, const double p_left = 0, const double p_right = 0) : _core(number), _p_left(p_left), _p_right(p_right) {}
         static double relation (const double x, const double core, const double p_left, const double p_right)
         {
            if (x < core) 
               return std::max (0.0, 1 - std::abs(x - core) / p_left);
            else 
               return std::max (0.0, 1 - std::abs(x - core) / p_right);
         }
         static double equal (const triangular_asymmetric_efn & l, const triangular_asymmetric_efn & r)
         {
            double d = std::abs(l._core - r._core);
            double p_left = std::max(l._p_left, r._p_left); 
            double p_right = std::max(l._p_right, r._p_right); 
            return relation(0.0, d, p_left, p_right);
         }
         static double less (const triangular_asymmetric_efn & l, const triangular_asymmetric_efn & r)
         {
            if (l._core > r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }
         }
         static double greater (const triangular_asymmetric_efn & l, const triangular_asymmetric_efn & r)
         {
            if (l._core < r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }
         }
         static double less_equal (const triangular_asymmetric_efn & l, const triangular_asymmetric_efn & r)
         {
            if (l._core > r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }
         static double greater_equal (const triangular_asymmetric_efn & l, const triangular_asymmetric_efn & r)
         {
            if (l._core < r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }
         //////////////////
         static triangular_asymmetric_efn zero () 
         {
            return {0.0, 0.0, 0.0};
         }
         static triangular_asymmetric_efn infinity () 
         {
            return {std::numeric_limits<double>::infinity(), 0.0, 0.0};
         }
         //////////////////
         triangular_asymmetric_efn operator * (const double & r) const
         {
            return triangular_asymmetric_efn (this->_core * r, this->_p_left, this->_p_right);
         }
         triangular_asymmetric_efn operator * (const triangular_asymmetric_efn & r) const
         {
            return triangular_asymmetric_efn (this->_core * r._core, std::max(this->_p_left, r._p_left), std::max(this->_p_right, r._p_right));
         }
         triangular_asymmetric_efn operator + (const triangular_asymmetric_efn & r) const
         {
            return triangular_asymmetric_efn (this->_core + r._core, std::max(this->_p_left, r._p_left), std::max(this->_p_right, r._p_right));
         }
         triangular_asymmetric_efn operator - (const triangular_asymmetric_efn & r) const
         {
            return triangular_asymmetric_efn (this->_core - r._core, std::max(this->_p_left, r._p_left), std::max(this->_p_right, r._p_right));
         }
         triangular_asymmetric_efn operator / (const triangular_asymmetric_efn & r) const
         {
            return triangular_asymmetric_efn (this->_core / r._core, std::max(this->_p_left, r._p_left), std::max(this->_p_right, r._p_right));
         }
         friend std::ostream & operator << (std::ostream & sos, const triangular_asymmetric_efn & t)
         {
            return sos << "(" << t._core << ", " << t._p_left << ", " << t._p_right << ")";
         }
         //////////////////
      protected:
         double _core; // core of the number
         double _p_left, _p_right; // left and right parts of the support (with regard to the core)

   };

   class gaussian_efn  : public extensional_fuzzy_number
   {
      public:
         gaussian_efn () : gaussian_efn(0.0, 0.0) {} 
         gaussian_efn (const double number, const double p = 0) : _core(number), _p(p) {}; 
         static double relation (const double x, const double core, const double p)
         {
            double x_core_2 = (x - core) * (x - core);
            double _2_s_2   = 2 * p * p;
            return std::exp (- x_core_2/ _2_s_2);
         }
         static double equal (const gaussian_efn & l, const gaussian_efn & r)
         {
            double d = std::abs(l._core - r._core);
            double p = std::max(l._p, r._p); 
            return relation(0.0, d, p);
         }
         static double less (const gaussian_efn & l, const gaussian_efn & r)
         {
            if (l._core > r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }
         }
         static double greater (const gaussian_efn & l, const gaussian_efn & r)
         {
            if (l._core < r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }
         }
         static double less_equal (const gaussian_efn & l, const gaussian_efn & r)
         {
            if (l._core > r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }
         static double greater_equal (const gaussian_efn & l, const gaussian_efn & r)
         {
            if (l._core < r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }
         //////////////////
         static gaussian_efn zero () 
         {
            return {0.0, 0.0};
         }
         static gaussian_efn infinity () 
         {
            return {std::numeric_limits<double>::infinity(), 0.0};
         }
         //////////////////
         gaussian_efn operator * (const double & r) const 
         {
            return gaussian_efn (this->_core * r, this->_p);
         }
         gaussian_efn operator * (const gaussian_efn & r) const 
         {
            return gaussian_efn (this->_core * r._core, std::max(this->_p, r._p));
         }
         gaussian_efn operator + (const gaussian_efn & r) const 
         {
            return gaussian_efn (this->_core + r._core, std::max(this->_p, r._p));
         }
         gaussian_efn operator - (const gaussian_efn & r) const 
         {
            return gaussian_efn (this->_core - r._core, std::max(this->_p, r._p));
         }
         gaussian_efn operator / (const gaussian_efn & r) const 
         {
            return gaussian_efn (this->_core / r._core, std::max(this->_p, r._p));
         }
         //////////////////
         friend std::ostream & operator << (std::ostream & sos, const gaussian_efn & t)
         {
            return sos << "(" << t._core << ", " << t._p << ")";
         }
      protected:
         double _core; // core of the number
         double _p; // half of the support

   };

   class expabs_efn  : public extensional_fuzzy_number
   {
      public:
         expabs_efn () : expabs_efn(0.0, 0.0) {} 
         expabs_efn (const double number, const double p = 0) : _core(number), _p(p) {}; 
         static double relation (const double x, const double core, const double p)
         {
            double abs_x_core = std::abs(x - core);
            return std::exp (- abs_x_core/ p);
         }
         static double equal (const expabs_efn & l, const expabs_efn & r)
         {
            double d = std::abs(l._core - r._core);
            double p = std::max(l._p, r._p); 
            return relation(0.0, d, p);
         }
         static double less (const expabs_efn & l, const expabs_efn & r)
         {
            if (l._core > r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }
         }
         static double greater (const expabs_efn & l, const expabs_efn & r)
         {
            if (l._core < r._core)
               return 0.0;
            else 
            {
               return 1 - equal (l, r);
            }
         }
         static double less_equal (const expabs_efn & l, const expabs_efn & r)
         {
            if (l._core > r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }
         static double greater_equal (const expabs_efn & l, const expabs_efn & r)
         {
            if (l._core < r._core)
               return equal(l, r);
            else 
            {
               return 1;
            }
         }
         //////////////////
         static expabs_efn zero () 
         {
            return {0.0, 0.0};
         }
         static expabs_efn infinity () 
         {
            return {std::numeric_limits<double>::infinity(), 0.0};
         }
         //////////////////
         expabs_efn operator * (const double & r) const 
         {
            return expabs_efn (this->_core * r, this->_p);
         }
         expabs_efn operator * (const expabs_efn & r) const 
         {
            return expabs_efn (this->_core * r._core, std::max(this->_p, r._p));
         }
         expabs_efn operator + (const expabs_efn & r) const 
         {
            return expabs_efn (this->_core + r._core, std::max(this->_p, r._p));
         }
         expabs_efn operator - (const expabs_efn & r) const 
         {
            return expabs_efn (this->_core - r._core, std::max(this->_p, r._p));
         }
         expabs_efn operator / (const expabs_efn & r) const 
         {
            return expabs_efn (this->_core / r._core, std::max(this->_p, r._p));
         }
         //////////////////
         friend std::ostream & operator << (std::ostream & sos, const expabs_efn & t)
         {
            return sos << "(" << t._core << ", " << t._p << ")";
         }
      protected:
         double _core; // core of the number
         double _p; // parameter

   };
}

namespace ksi 
{
   template <class T>
      void sort(typename std::vector<T> & numbers)
      {
         std::size_t size = numbers.size();
         for (std::size_t i = 1; i < size; i++)
         {
            T minimum = numbers[i];
            std::size_t j = i;
            while (j > 0 and minimum < numbers[j - 1]) 
            {
               numbers[j] = numbers[j - 1];
               j--;
            }
            numbers[j] = minimum;
         }
      }

   template <class T>
      void sort(typename std::vector<T> & numbers, const double xi)
      {
         std::size_t size = numbers.size();
         for (std::size_t i = 1; i < size; i++)
         {
            T minimum = numbers[i];
            std::size_t j = i;

            while (j > 0 and T::less(minimum, numbers[j - 1]) > xi) 
            {
               numbers[j] = numbers[j - 1];
               j--;
            }
            numbers[j] = minimum;
         }
      }

   template <typename T>
      T owa (const typename std::vector<T> & values, const std::vector<double> & weights)
      {
         if (values.size() != weights.size())
         {
            std::stringstream sos;
            sos << "Number of values (" << values.size() << ") does not match number of weights (" << weights.size() << ").";
            throw sos.str(); 
         }

         typename std::vector<T> X (values);
         ksi::sort(X);
         ksi::print (X); // tymczasowo
         T results;
         for (std::size_t i = 0; i < values.size(); ++i)
         {
            results = results + X[i] * weights[i];
         }
         return results;
      }

   template <typename T>
      T owa (const typename std::vector<T> & values, const std::vector<double> & weights, const double xi)
      {
         if (values.size() != weights.size())
         {
            std::stringstream sos;
            sos << "Number of values (" << values.size() << ") does not match number of weights (" << weights.size() << ").";
            throw sos.str(); 
         }

         typename std::vector<T> X (values);
         ksi::sort(X, xi);
         ksi::print (X); // tymczasowo
         T result;
         for (std::size_t i = 0; i < values.size(); ++i)
         {
            result = result + X[i] * weights[i];
         }
         return result;
      }

   template <typename T>
      class graph
      {
         public: 
            std::unordered_map<int, std::vector<std::pair<int, T>>> nodes;

         public:
            void add_edge (const int start_node, const int end_node, const T &  weight)
            {
               nodes[start_node].push_back({end_node,  weight});

               if (nodes.count(end_node) == 0)
                  nodes.insert({end_node, {}});
            }
         public:
            template <typename Y>
               friend std::ostream & operator << (std::ostream &, const graph<Y> &);
      };

   template <typename T>
      std::ostream & operator << (std::ostream & sos, const graph<T> & g)
      {
         for (const auto & n : g.nodes)
         {
            for (const auto & e : n.second)
            {
               sos << n.first << " -> " << e.first << " : " << e.second << std::endl;
            }
         }
         return sos;
      }

   template <typename T>
      void print (std::ostream & sos, const std::unordered_map<int, std::unordered_map<int, T>> & d) 
      {
         for (const auto [i, edges] : d)
         {
            sos << i << " | ";
            for (const auto [j, w] : edges)
            {
               sos << w << " ";             
            }
            sos << std::endl;
         }
      }

   template <typename T>
      void print_paths (std::ostream & sos, const std::unordered_map<int, std::unordered_map<int, T>> & d, const std::unordered_map<int, std::unordered_map<int, int>> & p)
      {
         // najpierw sczytuje indeksy wezlow:
         std::vector<int> node_indices;
         for (const auto [index, _] : d)
            node_indices.push_back(index);

         // wypisuje sciezki:
         for (const int i : node_indices)
         {
            for (const int j : node_indices)
            {
               if (i != j) 
               {
                  sos << "d[" << i <<  "][" << j << "]: " << d.at(i).at(j) << ", p: ";
                  sos << i;
                  int predecessor = p.at(i).at(j);
                  auto prev = i;
                  while (predecessor != prev and predecessor != j and predecessor != -1)
                  {
                     auto distance = d.at(prev).at(predecessor);
                     // sos << " ---[" << distance << "]--> " << predecessor;
                     sos << " ---" << distance << "--> " << predecessor;
                     prev = predecessor;
                     predecessor = p.at(predecessor).at(j);
                  }
                  auto distance =  d.at(prev).at(j);
                  sos << " ---[" << distance << "]--> " << j << std::endl;
               }
            }
         }
      }

   template <typename T>
      std::pair< std::unordered_map<int, std::unordered_map<int, T>>, std::unordered_map<int, std::unordered_map<int, int>>>
      floyd_warshall (const graph<T> & G, const double xi)
      {
         std::vector<int> node_indices;
         for (const auto n : G.nodes)
         {
            node_indices.push_back(n.first); 
         }
         // Mam indeksy wszystkich wierzcholkow.
         std::unordered_map<int, std::unordered_map<int, T>> d; // distances 
         std::unordered_map<int, std::unordered_map<int, int>> p; // predecessors 

         // inicjalizacja macierzy p i d:

         for (const int i : node_indices)
         {
            for (const int j : node_indices)
            {
               if  (i == j)
               {
                  d[i][j] = T::zero();  
                  p[i][j] = j;
               }
               else 
               {
                  d[i][j] = T::infinity();
                  p[i][j] = std::numeric_limits<int>::min(); // nothing 
               }
            }
         }
         for (const auto [i, edges] : G.nodes)
         {
            for (const auto [j, w] : edges)
            {
               d[i][j] = w;
               p[i][j] = j;
            }
         }

         ///////
         {
            std::cout << "distance matrix" << std::endl; 
            for (const auto [i, edges] : d)
            {
               std::cout << i << " | ";
               for (const auto [j, w] : edges)
               {
                  std::cout << w << " ";             
               }
               std::cout << std::endl;
            }
         }
         {
            std::cout << "predecessor matrix" << std::endl; 
            for (const auto [i, edges] : p)
            {
               std::cout << i << " | ";
               for (const auto [j, w] : edges)
               {
                  std::cout << w << " ";             
               }
               std::cout << std::endl;
            }
         }
         ///////////////////////////// 

         // no to jestemy gotowi!
         {
            for (const int k : node_indices)
            {
               for (const int i : node_indices)
               {
                  for (const int j : node_indices)
                  {
                     T suma = d[i][k] + d[k][j];
                     if (T::less(suma, d[i][j]) > xi)
                     {
                        d[i][j] = suma;
                        p[i][j] = p[i][k];
                     }
                  }
               }
            }
         }
         return {d, p};
      }
}

namespace ksi 
{
   template <typename T>
   class linear_regression
   {
      std::vector<std::vector<T>> transpose (const std::vector<std::vector<T>> & matrix)
      {
         auto num_rows = matrix.size();
         auto num_cols = matrix[0].size(); 

         std::vector<std::vector<T>> transposed (num_cols, std::vector<T>(num_rows));
         for (std::size_t i = 0; i < num_rows; ++i)
         {
            for (std::size_t j = 0; j < num_cols; ++j)
            {
               transposed[j][i] = matrix[i][j];
            }
         }
         return transposed;
      }

      std::vector<std::vector<T>> multiply (const std::vector<std::vector<T>> & A, const std::vector<std::vector<T>> & B)
      {
         auto num_rows_A = A.size();
         auto num_cols_A = A[0].size();
         auto num_rows_B = B.size();
         auto num_cols_B = B[0].size();

         if (num_cols_A != num_rows_B)
         {
            throw "Number of columns in A must match number of rows in B.";
         }

         std::vector<std::vector<T>> result (num_rows_A, std::vector<T>(num_cols_B, T()));
         for (std::size_t i = 0; i < num_rows_A; ++i)
         {
            for (std::size_t j = 0; j < num_cols_B; ++j)
            {
               for (std::size_t k = 0; k < num_cols_A; ++k)
               {
                  result[i][j] = result[i][j] + A[i][k] * B[k][j];
               }
            }
         }
         return result;
      }

      std::vector<std::vector<T>> invert (const std::vector<std::vector<T>> & matrix)
      {
         // inversja macierzy, np. metodą Gaussa-Jordana 
         // metoda Gaussa-Jordana: 
         
         auto num_rows = matrix.size();
         auto extented_matrix = matrix; // rozszerzona macierz [A | I]
         for (std::size_t i = 0; i < num_rows; ++i)
         {
            extented_matrix[i].resize(2 * num_rows, T());
            extented_matrix[i][num_rows + i] = T() + T(1); // dodajemy macierz jednostkową
         }
         // etap drugi: redukcja do postaci schodkowej 
         for (std::size_t i = 0; i < num_rows; ++i)
         {
            // eliminacja Gaussa
            for (std::size_t j = i + 1; j < num_rows; ++j)
            {
               T factor = extented_matrix[j][i] / extented_matrix[i][i];
               for (std::size_t k = i; k < 2 * num_rows; ++k)
               {
                  extented_matrix[j][k] = extented_matrix[j][k] - factor * extented_matrix[i][k];
               }
            }
         }
         // etap trzeci: redukcja do postaci diagonalnej  
         for (std::size_t i = num_rows; i > 0; --i)
         {
            // eliminacja Gaussa-Jordana
            for (std::size_t j = i - 1; j > 0; --j)
            {
               T factor = extented_matrix[j - 1][i - 1] / extented_matrix[i - 1][i - 1];
               for (std::size_t k = i - 1; k < 2 * num_rows; ++k)
               {
                  extented_matrix[j - 1][k] = extented_matrix[j - 1][k] - factor * extented_matrix[i - 1][k];
               }
            }
         }

         // etap czwarty: normalizacja do postaci jednostkowej
         for (std::size_t i = 0; i < num_rows; ++i)
         {
            T factor = extented_matrix[i][i];
            for (std::size_t k = i; k < 2 * num_rows; ++k)
            {
               extented_matrix[i][k] = extented_matrix[i][k] / factor;
            }
         }
         
         // etap piąty: wyodrębnienie macierzy odwrotnej 

         for (std::size_t i = 0; i < num_rows; ++i)
         {
            for (std::size_t j = 0; j < num_rows; ++j)
            {
               extented_matrix[i][j] = extented_matrix[i][num_rows + j];
            }
            extented_matrix[i].resize(num_rows);
         }

         return extented_matrix; 
      }


      public: 
         std::vector<std::vector<T>> elaborate_coefficients (const std::vector<std::vector<T>> & points, const std::vector<std::vector<T>> & output)
         {
            // points: matrix of input points, 
            // output: vector of output values
            // maths formula : points * coefficients = output 
            //
            std::vector<T> coefficients;
            // maths: coefficients = (points^T * points)^(-1) * points^T * output  
            auto pointsT = transpose(points);
            auto pointsT_points = multiply(pointsT, points);
            auto pointsT_points_inv = invert(pointsT_points);
            auto pointsT_output = multiply(pointsT, output);

            return multiply(pointsT_points_inv, pointsT_output);
         }

   };
}

//---------------------
void example_1()
{
  {
         std::cout << "# example 1" << std::endl; 

         // Floyd-Warshall algorithm 
         ksi::graph<ksi::triangular_asymmetric_efn> G; // graph 

         /////////////////
         G.add_edge(	 1, 2, {820,20,20 });
         G.add_edge(	 1, 3, {361,11,9 });
         G.add_edge(	 1, 6, {677,27,6 });
         G.add_edge(	 1, 9, {300,10,50 });
         G.add_edge(	 1, 10, {450,30,20 });
         G.add_edge(	 2, 3, {186,6,7 });
         G.add_edge(	 2, 5, {510,15,15 });
         G.add_edge(	 2, 9, {930,30,30 });
         G.add_edge(	 3, 4, {667,17,216 });
         G.add_edge(	 3, 5, {748,18,22 });
         G.add_edge(	 3, 8, {443,18,22 });
         G.add_edge(	 4, 5, {199,9,11 });
         G.add_edge(	 4, 6, {340,30,20 });
         G.add_edge(	 4, 11, {740,30,30 });
         G.add_edge(	 5, 6, {660,50,30 });
         G.add_edge(	 6, 11, {242,12,18 });
         G.add_edge(	 7, 6, {410,20,30 });
         G.add_edge(	 7, 11, {472,22,18 });
         G.add_edge(	 8, 4, {730,20,5 });
         G.add_edge(	 8, 7, {242,12,13 });
         G.add_edge(	 9, 8, {137,7,8 });
         G.add_edge(	 9, 7, {130,10,20 });
         G.add_edge(	 9, 10, {242,12,18 });
         G.add_edge(	 10, 7, {342,12,8 });
         G.add_edge(	 10, 11, {1310, 60, 130});
         /////////////////

         std::cout << "graph:" << std::endl;
         std::cout << G << std::endl;

         const double xi {0.3};
         auto [d, p] = ksi::floyd_warshall(G, xi);
         std::cout << "xi: " << xi << std::endl;
         std::cout << "distance matrix" << std::endl;
         ksi::print (std::cout, d);
         std::cout << "predecessor matrix" << std::endl;
         ksi::print (std::cout, p);
         std::cout  << "paths" << std::endl;
         print_paths(std::cout, d, p);
         std::cout << std::endl;
         std::cout << "----------------------" << std::endl;
         std::cout << std::endl;

      } 
}
void example_2()
{
  {
         std::cout << "# example 2" << std::endl; 

         // Floyd-Warshall algorithm 
         ksi::graph<ksi::triangular_asymmetric_efn> G; // graph 
         G.add_edge(1, 2, { 2, 1, 1});
         G.add_edge(1, 3, { 7, 2, 2});
         G.add_edge(2, 3, { 4, 3, 5});
         G.add_edge(2, 4, { 11, 1, 1});
         G.add_edge(2, 5, { 6, 1, 1});
         G.add_edge(3, 4, { 9, 1, 1});
         G.add_edge(4, 5, { -8, 1, 1});
         G.add_edge(4, 6, { 13, 2, 1});
         G.add_edge(5, 6, { 9, 1, 1});

         std::cout << "graph:" << std::endl;
         std::cout << G << std::endl;

         const double xi {0.3};
         auto [d, p] = ksi::floyd_warshall(G, xi);
         std::cout << "xi: " << xi << std::endl;
         std::cout << "distance matrix" << std::endl;
         ksi::print (std::cout, d);
         std::cout << "predecessor matrix" << std::endl;
         ksi::print (std::cout, p);
         std::cout  << "paths" << std::endl;
         print_paths(std::cout, d, p);
         std::cout << std::endl;
         std::cout << "----------------------" << std::endl;
         std::cout << std::endl;

      } 
}
//---------------------
void fuzzy_operations()
{
   ksi::fuzzy_set_triangle A (2, 3, 5);
   ksi::fuzzy_set_triangle B (6, 7, 9);
   ksi::tnorm_product t;
   

   {
      std::cout << "------------------------------" << std::endl; 
      auto product = [](auto a, auto b){return a * b;};

      auto C = fuzzy_operation(A, B, t, product, 1.0, 10.0, 0.011);

      for (const auto [x, m] : C)
      {
         std::cout << std::fixed << std::setprecision(6) << x << " " << m << std::endl;
      }
      std::cout << "------------------------------" << std::endl; 
   }
   {
      std::cout << "------------------------------" << std::endl; 
      auto suma = [](auto a, auto b){return a + b;};

      auto C = fuzzy_operation(A, B, t, suma, 1.0, 10.0, 0.011);

      for (const auto [x, m] : C)
      {
         std::cout << std::fixed << std::setprecision(6) << x << " " << m << std::endl;
      }
      std::cout << "------------------------------" << std::endl; 
   }
}

void example_comparison ()
{
   {
      {
         std::cout << "gaussian EFN" << std::endl;

         std::vector<double> fuzzifications;
         for (double f = 0; f < 10; f += 1)
            fuzzifications.push_back(f);
         for (double f = 10; f < 50; f += 5)
            fuzzifications.push_back(f);
         for (double f = 50; f < 110; f += 10)
            fuzzifications.push_back(f);

         for (double f : fuzzifications)
         {
            ksi::gaussian_efn _70 { 70, f};
            ksi::gaussian_efn _80 { 80, f};
            std::cout << "" << _70 << " = " << _80 << " : " << ksi::gaussian_efn::equal(_70, _80) << "\t";
            std::cout << "" << _70 << " < " << _80 << " : " << ksi::gaussian_efn::less(_70, _80) << std::endl;
         }
      }
      {
         std::cout << "expabs_efn EFN" << std::endl;

         std::vector<double> fuzzifications;
         for (double f = 0; f < 10; f += 1)
            fuzzifications.push_back(f);
         for (double f = 10; f < 50; f += 5)
            fuzzifications.push_back(f);
         for (double f = 50; f < 110; f += 10)
            fuzzifications.push_back(f);

         for (double f : fuzzifications)
         {
            ksi::expabs_efn _70 { 70, f};
            ksi::expabs_efn _80 { 80, f};
            std::cout << "" << _70 << " = " << _80 << " : " << ksi::expabs_efn::equal(_70, _80) << "\t";
            std::cout << "" << _70 << " < " << _80 << " : " << ksi::expabs_efn::less(_70, _80) << std::endl;
         }
      }
      {
         std::cout << "triangular_efn EFN" << std::endl;

         std::vector<double> fuzzifications;
         for (double f = 0; f < 10; f += 1)
            fuzzifications.push_back(f);
         for (double f = 10; f < 50; f += 5)
            fuzzifications.push_back(f);
         for (double f = 50; f < 110; f += 10)
            fuzzifications.push_back(f);

         for (double f : fuzzifications)
         {
            ksi::triangular_efn _70 { 70, f};
            ksi::triangular_efn _80 { 80, f};
            std::cout << "" << _70 << " = " << _80 << " : " << ksi::triangular_efn::equal(_70, _80) << "\t";
            std::cout << "" << _70 << " < " << _80 << " : " << ksi::triangular_efn::less(_70, _80) << std::endl;
         }
      }
      {
         std::cout << "trapezoidal_efn EFN" << std::endl;

         std::vector<double> fuzzifications;
         for (double f = 0; f < 10; f += 1)
            fuzzifications.push_back(f);
         for (double f = 10; f < 50; f += 5)
            fuzzifications.push_back(f);
         for (double f = 50; f < 110; f += 10)
            fuzzifications.push_back(f);

         for (double f : fuzzifications)
         {
            ksi::trapezoidal_efn _70 { 70, f};
            ksi::trapezoidal_efn _80 { 80, f};
            std::cout << "" << _70 << " = " << _80 << " : " << ksi::trapezoidal_efn::equal(_70, _80) << "\t";
            std::cout << "" << _70 << " < " << _80 << " : " << ksi::trapezoidal_efn::less(_70, _80) << std::endl;
         }
      }
      {
         const int PREC {4};
         const int WID  {8};

         std::vector<double> fuzzifications;
         for (double f = 0; f < 20; f += 1)
            fuzzifications.push_back(f);
         for (double f = 20; f < 50; f += 5)
            fuzzifications.push_back(f);
         for (double f = 50; f < 110; f += 10)
            fuzzifications.push_back(f);

         std::cout << "fuzzification triangular trapezoidal gaussian expabs" << std::endl;
         for (double f : fuzzifications)
         {
            ksi::triangular_efn  t_70 { 70, f};
            ksi::triangular_efn  t_80 { 80, f};
            ksi::trapezoidal_efn r_70 { 70, f};
            ksi::trapezoidal_efn r_80 { 80, f};
            ksi::gaussian_efn    g_70 { 70, f};
            ksi::gaussian_efn    g_80 { 80, f};
            ksi::expabs_efn      e_70 { 70, f};
            ksi::expabs_efn      e_80 { 80, f};


            std::cout << std::setw(WID) << std::setprecision(PREC) << std::fixed << f 
               << "   " << ksi::triangular_efn::equal(t_70, t_80) 
               << "   " << ksi::trapezoidal_efn::equal(r_70, r_80) 
               << "   " << ksi::gaussian_efn::equal(g_70, g_80) 
               << "   " << ksi::expabs_efn::equal(e_70, e_80) 
               << std::endl;
         }
      }

   }

}

void example_sorting()
{
   std::cout << std::endl;
   std::cout << "======================" << std::endl;
   std::cout << "EXAMPLE: SORTING" << std::endl;
   std::cout << "======================" << std::endl;
   std::cout << std::endl;

   {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "type: int" << std::endl;
      std::cout << "-------------------------------------------" << std::endl;
      std::vector<int> A {4, -4, 5, 0, 7,4}; 
      std::cout << "input:" << std::endl;
      ksi::print(A);
      std::cout << "output:" << std::endl;
      ksi::sort(A);
      ksi::print(A);
   }

   {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "type: double" << std::endl;
      std::cout << "-------------------------------------------" << std::endl;
      std::vector<double> A {4.5, -4.3, 5.6, 0.4, 7.7,4.2}; 
      std::cout << "input:" << std::endl;
      ksi::print(A);
      std::cout << "output:" << std::endl;
      ksi::sort(A);
      ksi::print(A);
   }

   {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "type: ksi::triangular_efn" << std::endl;
      std::cout << "-------------------------------------------" << std::endl;
      std::vector<ksi::triangular_efn> A {{4.5, 1.0}, {-4.3, 1.0}, {5.6, 1.0}, {0.4, 1.0}, {7.7, 1.0}, {4.2, 1.0}}; 
      std::cout << "input:" << std::endl;
      ksi::print(A);
      std::cout << "output for parameter xi:" << std::endl;
      for (double xi = 0; xi <= 1.0; xi += 0.1)
      {
         auto K = A;
         debug(xi);
         ksi::sort(K, xi);
         ksi::print(K);
      }
   }
   {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "type: ksi::gaussian_efn" << std::endl;
      std::cout << "-------------------------------------------" << std::endl;
      std::vector<ksi::gaussian_efn> A {{4.5, 1.0}, {-4.3, 1.0}, {5.6, 1.0}, {0.4, 1.0}, {7.7, 1.0}, {4.2, 1.0}}; 
      std::cout << "input:" << std::endl;
      ksi::print(A);
      std::cout << "output for parameter xi:" << std::endl;
      for (double xi = 0; xi <= 1.0; xi += 0.1)
      {
         auto K = A;
         debug(xi);
         ksi::sort(K, xi);
         ksi::print(K);
      }
   }
   {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "type: ksi::gaussian_efn" << std::endl;
      std::cout << "-------------------------------------------" << std::endl;
      std::vector<ksi::gaussian_efn> A {{5.5, 0.5}, {6.0, 3.0}, {4.5, 1.0}, {-4.3, 1.0}, {5.6, 1.0}, {0.4, 2.0}, {7.7, 1.0}, {4.2, 1.5}}; 
      std::cout << "input:" << std::endl;
      ksi::print(A);
      std::cout << "output for parameter xi:" << std::endl;
      for (double xi = 0; xi <= 1.0; xi += 0.1)
      {
         auto K = A;
         debug(xi);
         ksi::sort(K, xi);
         ksi::print(K);
      }
   }

}

void example_owa()
{
   {
      std::cout << "-------------------------------------------" << std::endl;
      std::vector<double> X {4.5, -3.2, 5.7, 2.6, 0.5};
      std::vector<double> W {0.1,  0.2, 0.2, 0.1, 0.4};
      try 
      {
         double result = ksi::owa<double>(X, W);
         ksi::print(X);
         ksi::print(W);
         debug(result);
      }
      catch (const std::string & str)
      {
         std::cout << str << std::endl;
      }
   }

   {
      std::cout << "-------------------------------------------" << std::endl;
      std::vector<ksi::gaussian_efn> X {{4.5, 0.4}, {-3.2, 0.1}, {5.7, 0.3}, {2.6, 0.1}, {0.5, 0.2}};
      std::vector<double> W {0.1,  0.2, 0.2, 0.1, 0.4};
      try 
      {
         const double xi {0.3}; 
         auto result = ksi::owa(X, W, xi);
         ksi::print(X);
         ksi::print(W);
         debug(xi);
         debug(result);
      }
      catch (const std::string & str)
      {
         std::cout << str << std::endl;
      }
   }

   {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "minimum" << std::endl;
      std::cout << "-------------------------------------------" << std::endl;

      std::vector<ksi::gaussian_efn> X {{4.5, 0.4}, {-3.1, 1.1}, {-3.2, 0.1}, {5.7, 0.3}, {2.5, 0.2}, {2.8, 0.1}, {0.5, 0.2}, {2.7, 0.3}, {2.6, 0.1}};
      std::vector<double> W            {1.0,  0, 0, 0, 0, 0, 0, 0, 0};
      try 
      {
         ksi::print(X);
         ksi::print(W);
         for (double xi = 0.0; xi <= 1.0; xi += 0.1)
         {
            auto result = ksi::owa(X, W, xi);
            debug(xi);
            debug(result);
         }
      }
      catch (const std::string & str)
      {
         std::cout << str << std::endl;
      }

   }
   {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "average" << std::endl;
      std::cout << "-------------------------------------------" << std::endl;

      std::vector<ksi::gaussian_efn> X {{4.5, 0.4}, {-3.1, 1.1}, {-3.2, 0.1}, {5.7, 0.3}, {2.5, 0.2}, {2.8, 0.1}, {0.5, 0.2}, {2.7, 0.3}, {2.6, 0.1}};
      std::vector<double> W            {1.0/9.0,  1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0};
      try 
      {
         ksi::print(X);
         ksi::print(W);
         for (double xi = 0.0; xi <= 1.0; xi += 0.1)
         {
            auto result = ksi::owa(X, W, xi);
            debug(xi);
            debug(result);
         }
      }
      catch (const std::string & str)
      {
         std::cout << str << std::endl;
      }
   }

   {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "median" << std::endl;
      std::cout << "-------------------------------------------" << std::endl;

      std::vector<ksi::gaussian_efn> X {{4.5, 0.4}, {-3.1, 1.1}, {-3.2, 0.1}, {5.7, 0.3}, {2.5, 0.2}, {2.8, 0.1}, {0.5, 0.2}, {2.7, 0.3}, {2.6, 0.1}};
      std::vector<double> W            {0,  0, 0, 0, 1, 0, 0, 0, 0};
      try 
      {
         ksi::print(X);
         ksi::print(W);
         for (double xi = 0.0; xi <= 1.0; xi += 0.1)
         {
            auto result = ksi::owa(X, W, xi);
            debug(xi);
            debug(result);
         }
      }
      catch (const std::string & str)
      {
         std::cout << str << std::endl;
      }
   }

   {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "random numbers" << std::endl;
      std::cout << "-------------------------------------------" << std::endl;

      const int n {11};
      std::vector<ksi::gaussian_efn> X;
      std::vector<double> W; 

      std::default_random_engine silnik (std::chrono::system_clock::now().time_since_epoch().count());
      std::uniform_real_distribution<double> core_distro (-10, 10);
      std::uniform_real_distribution<double> fuzzyfication_distro (0.01, 3);

      double sum {0.0};
      for (int i = 0; i < n; ++i)
      {
         X.push_back({core_distro(silnik), fuzzyfication_distro(silnik)});
         std::uniform_real_distribution<double> weight_distro (0.0, 1.0 - sum);
         auto liczba = weight_distro(silnik);
         sum += liczba;
         W.push_back(liczba);
      }
      std::shuffle(W.begin(), W.end(), silnik);

      debug(sum);
      for (auto & w : W)
      {
         w /= sum; 
      }
      try 
      {
         ksi::print(X);
         ksi::print(W);
         for (double xi = 0.0; xi <= 1.0; xi += 0.1)
         {
            auto result = ksi::owa(X, W, xi);
            debug(xi);
            debug(result);
         }
      }
      catch (const std::string & str)
      {
         std::cout << str << std::endl;
      }
   }

}

void example_floydwarshall()
{
   {
      std::cout << "### example 1" << std::endl; 

      // Floyd-Warshall algorithm 
      ksi::graph<ksi::triangular_efn> G; // graph 
      G.add_edge(1, 2, {25, 1});
      G.add_edge(1, 3, {15, 1});
      G.add_edge(2, 1, {12, 1});
      G.add_edge(2, 3, {19, 2});
      G.add_edge(3, 1, {10, 3});
      G.add_edge(3, 4, {12, 2});
      G.add_edge(4, 2, {8, 1});
      G.add_edge(4, 1, {9, 2});

      std::cout << "graph:" << std::endl;
      std::cout << G << std::endl;

      const double xi {0.3};
      auto [d, p] = ksi::floyd_warshall(G, xi);
      std::cout << "distance matrix" << std::endl;
      ksi::print (std::cout, d);
      std::cout << "predecessor matrix" << std::endl;
      ksi::print (std::cout, p);
      std::cout << "paths" << std::endl;
      print_paths(std::cout, d, p);
      std::cout << std::endl;
      std::cout << "----------------------" << std::endl;
      std::cout << std::endl;

   }
   {
      std::cout << "### example 2" << std::endl; 

      // Floyd-Warshall algorithm 
      ksi::graph<ksi::triangular_efn> G; // graph 
      G.add_edge(1, 3, { 6, 8});
      G.add_edge(2, 1, { 9, 5});
      G.add_edge(2, 4, { 2, 9});
      G.add_edge(3, 4, {13, 9});
      G.add_edge(3, 6, { 5, 7});
      G.add_edge(4, 1, { 5, 9});
      G.add_edge(4, 6, { 17, 8});
      G.add_edge(4, 7, { 3, 7});
      G.add_edge(5, 2, { 5, 8});
      G.add_edge(5, 4, { 8, 8});
      G.add_edge(6, 7, {12, 7});
      G.add_edge(7, 5, { 5, 8});

      std::cout << "graph:" << std::endl;
      std::cout << G << std::endl;

      const double xi {0.1};
      auto [d, p] = ksi::floyd_warshall(G, xi);
      std::cout << "xi: " << xi << std::endl;
      std::cout << "distance matrix" << std::endl;
      ksi::print (std::cout, d);
      std::cout << "predecessor matrix" << std::endl;
      ksi::print (std::cout, p);
      std::cout  << "paths" << std::endl;
      print_paths(std::cout, d, p);
      std::cout << std::endl;
      std::cout << "----------------------" << std::endl;
      std::cout << std::endl;

   }

   {
      std::cout << "### example 3" << std::endl; 

      // Floyd-Warshall algorithm 
      ksi::graph<ksi::gaussian_efn> G; // graph 
      G.add_edge(1, 3, { 6, 8});
      G.add_edge(2, 1, { 9, 5});
      G.add_edge(2, 4, { 2, 9});
      G.add_edge(3, 4, {13, 9});
      G.add_edge(3, 6, { 5, 7});
      G.add_edge(4, 1, { 5, 9});
      G.add_edge(4, 6, { 17, 8});
      G.add_edge(4, 7, { 3, 7});
      G.add_edge(5, 2, { 5, 8});
      G.add_edge(5, 4, { 8, 8});
      G.add_edge(6, 7, {12, 7});
      G.add_edge(7, 5, { 5, 8});

      std::cout << "graph:" << std::endl;
      std::cout << G << std::endl;

      const double xi {0.5};
      auto [d, p] = ksi::floyd_warshall(G, xi);
      std::cout << "xi: " << xi << std::endl;
      std::cout << "distance matrix" << std::endl;
      ksi::print (std::cout, d);
      std::cout << "predecessor matrix" << std::endl;
      ksi::print (std::cout, p);
      std::cout  << "paths" << std::endl;
      print_paths(std::cout, d, p);
      std::cout << std::endl;
      std::cout << "----------------------" << std::endl;
      std::cout << std::endl;

   }

}

void example_linear_regression()
{
   std::cout << std::endl;
   std::cout << "==========================" << std::endl;
   std::cout << "EXAMPLE: LINEAR REGRESSION" << std::endl;
   std::cout << "==========================" << std::endl;
   std::cout << std::endl;

   {
      std::cout << "double" << std::endl;
      using T = double; 
      ksi::linear_regression<T> lr; 

      std::vector<std::vector<T>> points {{2, 1}, {3, 1}, {5, 1}, {7, 1}, {11, 1}}; 
      std::vector<std::vector<T>> output {{1}, {2}, {3}, {4}, {5}}; 

      auto coefficients = lr.elaborate_coefficients(points, output); 
      std::cout << "points:" << std::endl; 
      ksi::print(points); 
      std::cout << "output:" << std::endl; 
      ksi::print(output); 
      std::cout << "coefficients:" << std::endl; 
      ksi::print(coefficients); 

      std::cout << "----------------------" << std::endl; 
   }

   {
      std::cout << "ksi::gaussian_efn" << std::endl;
      using T = ksi::gaussian_efn;
      ksi::linear_regression<T> lr;

      std::vector<std::vector<T>> points {{ {2, 0.1}, {1, 0.1} }, { {3, 0.2}, {1, 0.2} }, { {5, 0.3}, {1, 0.1} }, { {7, 0.4}, {1, 0.1} }, { {11, 0.4}, {1, 0.2} }};
      std::vector<std::vector<T>> output {{ {1, 0.2} }, { {2, 0.1} }, { {3, 0.2} }, { {4, 0.1} }, { {5, 0.3} }};

      auto coefficients = lr.elaborate_coefficients(points, output);
      std::cout << "points:" << std::endl;
      ksi::print(points);
      std::cout << "output:" << std::endl;
      ksi::print(output);
      std::cout << "coefficients:" << std::endl;
      ksi::print(coefficients);

      std::cout << "----------------------" << std::endl;
   }
   {
      std::cout << "ksi::gaussian_efn" << std::endl;
      using T = ksi::gaussian_efn;
      ksi::linear_regression<T> lr;

      std::vector<std::vector<T>> points 
      {
         { {1.0, 0.1}, {1, 0.1} },
            { {2.0, 0.2}, {1, 0.2} }, 
            { {3.0, 0.3}, {1, 0.1} }, 
            { {4.0, 0.4}, {1, 0.1} }, 
            { {5.0, 0.4}, {1, 0.2} }
      };
      std::vector<std::vector<T>> output 
      {
         { {1.5, 0.2} }, 
            { {4.5, 0.1} }, 
            { {5.5, 0.2} }, 
            { {8.2, 0.1} }, 
            { {9.5, 0.3} }
      };

      auto coefficients = lr.elaborate_coefficients(points, output);
      std::cout << "points:" << std::endl;
      ksi::print(points);
      std::cout << "output:" << std::endl;
      ksi::print(output);
      std::cout << "coefficients:" << std::endl;
      ksi::print(coefficients);

      std::cout << "----------------------" << std::endl;
   }
}

//---------------------
int main ()
{
   // fuzzy_operations();
   // example_1();
   // example_2();
   // return 0;

   const bool comparison {true};
   const bool sortowanie {true};
   const bool floydwarshall {true};
   const bool owa {true};
   const bool linear_regression {true};

   if (comparison)
      example_comparison();
   if (sortowanie)
      example_sorting();
   if (floydwarshall)
      example_floydwarshall();
   if (owa)
      example_owa();
   if (linear_regression)
      example_linear_regression();
   return 0;
}

