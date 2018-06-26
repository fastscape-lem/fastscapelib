#pragma once

#include <fstream>

template <class T, class S>
void dbg_out(std::string fname, int num, T& array, const S& shape)
{
    std::stringstream ss;
    ss << fname << num << ".bin";
    std::ofstream file(ss.str(), std::ios_base::binary);

	if (!file)
	{
		std::cerr << "Impossible to open file " << ss.str() << std::endl;
		return;
	}

    size_t nrow = shape[0];
    size_t ncols = shape[1];
    file.write((char*)&nrow, sizeof(size_t));
    file.write((char*)&ncols, sizeof(size_t));

    size_t count = 1;
    for(size_t i = 0; i<shape.size(); ++i)	
        count*=shape[i];

    for(size_t i = 0; i< count; ++i)
    {
        double value = array(i);
        file.write((char*)&value, sizeof(double));
    }
}
