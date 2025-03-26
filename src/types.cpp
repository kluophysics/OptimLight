#include "types.hpp"

#include <map>
#include <string>

namespace OptimLight
{
    std::string  algorithm_to_string(algorithm algorithm)
    {
        switch(algorithm)
        {
            case ALGORITHM_SD: return "Steepest Descent";
            case ALGORITHM_CG: return "Conjugate Gradient";
            case ALGORITHM_BFGS: return "Broyden-Fletcher-Goldfarb-Shanno";
            case ALGORITHM_LBFGS: return "Limited-memory BFGS";
            case ALGORITHM_TNEWTON: return "Truncated Newton";
            case ALGORITHM_NEWTON: return "Newton";
            case ALGORITHM_TR: return "Trust Region";
            case ALGORITHM_LENGTH: return "Unknown";
        }
        return NULL;
    }
    algorithm algorithm_from_string(std::string & name)
    {
        int i;
        if (name == "Unknown" )
            return static_cast<algorithm>(-1);
        for (i = 0; i < ALGORITHM_LENGTH; i++)
        {
            if( name == algorithm_to_string( static_cast<algorithm>(i) ) ) 
                return static_cast<algorithm>(i);
        }
        return static_cast<algorithm>(-1);
    }


   

} // namespace OptimLight