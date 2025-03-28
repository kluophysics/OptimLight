#include "types.hpp"

#include <iostream>
#include <map>
#include <string>

namespace OptimLight
{
    std::string  algorithm_to_string(const Algorithm & algorithm)
    {
        switch(algorithm)
        {
            case Algorithm::ALGORITHM_SD: return "Steepest Descent";
            case  Algorithm::ALGORITHM_CG: return "Conjugate Gradient";
            case  Algorithm::ALGORITHM_BFGS: return "Broyden-Fletcher-Goldfarb-Shanno";
            case  Algorithm::ALGORITHM_LBFGS: return "Limited-memory BFGS";
            case  Algorithm::ALGORITHM_TNEWTON: return "Truncated Newton";
            case  Algorithm::ALGORITHM_NEWTON: return "Newton";
            case  Algorithm::ALGORITHM_TR: return "Trust Region";
        }
        return "";
    }

    Algorithm algorithm_from_string(std::string & name)
    {
        for (int i = 0; i <  sizeof(Algorithm); i++)
        {
            if( name == algorithm_to_string( static_cast<Algorithm>(i) ) ) 
                return static_cast<Algorithm>(i);
        }
        std::cout << "Algorithm name not found" << std::endl;
        return Algorithm::ALGORITHM_SD;
    }


   

} // namespace OptimLight