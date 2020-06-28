#include <fstream>
#include "is_empty.h"

bool is_empty(std::fstream& pFile)
{
    return pFile.peek() == std::ofstream::traits_type::eof();
}
