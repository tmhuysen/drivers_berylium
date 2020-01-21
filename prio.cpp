/**
 *  Berelyium0 Driver
 */

/**
 *  Auxiliary Structure
 */ 


#include <cmath>
#include <iomanip>      // std::setprecision
#include <queue>
#include <iostream>

struct Pair {
    double coeff;
    size_t index;
};

class myComparator 
{ 
public: 
    int operator() (const Pair& p1, const Pair& p2) 
    { 
        return std::abs(p1.coeff) > std::abs(p2.coeff); 
    } 
}; 
  

int main() {
   // MIN HEAP
   std::priority_queue <Pair, std::vector<Pair>, myComparator > pq;
   pq.push({0.1, 1}); 
   pq.push({0.2, 1}); 
   pq.push({0.5, 1}); 
   pq.push({0.3, 1}); 
   pq.push({0.4, 1}); 


   pq.push({0.6, 1});
   pq.pop();
   pq.push({0.6, 1});
   pq.pop();
   pq.push({0.6, 1});
   pq.pop();
   pq.push({0.6, 1});
   pq.pop();
   pq.push({0.1, 1});
   pq.pop();

   while(!pq.empty()) {
      const Pair& x = pq.top();
      std::cout << x.index << " : "<< x.coeff << std::endl;
      pq.pop();
   }

   return 0;
   }
