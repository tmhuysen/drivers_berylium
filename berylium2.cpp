/**
 *  Berelyium0 Driver
 */

/**
 *  Auxiliary Structure
 */ 

#include <iomanip>      // std::setprecision
#include <gqcp/gqcp.hpp>
#include <queue>
struct Pair {
    double coeff;
    size_t index;
};

class myComparator 
{ 
public: 
    int operator() (const Pair& p1, const Pair& p2) 
    { 
        return std::abs(p1.coeff) >std::abs(p2.coeff); 
    } 
}; 
  

int main() {	}
