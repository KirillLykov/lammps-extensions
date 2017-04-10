#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <tuple>
using namespace std;

#include <boost/algorithm/string.hpp>
using namespace boost;

#include <box.h>
#include <serializerHDF5.h>
#include <basic_access_strategy.h>
using namespace ls;
using namespace geometry_utils;

typedef ls::BasicReadAccessStrategy< double > VRAStrategy;
typedef ls::Grid3D< double > VGrid;
typedef io::GridSerializerHDF5< VGrid, VRAStrategy > VSerializerHDF5;

void computeIndex(const VGrid& grid, const double relativePosition[3], size_t* index) 
{
    assert(relativePosition[0] >= 0.0 && relativePosition[1] >= 0.0 && relativePosition[2] >= 0.0);
    // index_x = floor( p_x / h_x )
    for (size_t i = 0; i < 3; ++i) {
      // coef := 1 / h
      double coef = (grid.size(i) - 1.0) / static_cast<double>(grid.getBoundingBox().getIthSize(i));
      index[i] = static_cast<int>(relativePosition[i] * coef);
      assert(index[i] < grid.size(i));
    }
}

void parseLammpsAveSpatial(size_t targetFrame, const std::string& fileName, double h[3], VGrid& grid, const string& datasets)
{
    ifstream file(fileName.c_str());
    if (!file.is_open()) throw runtime_error("Could not open file");

    double leftBottom[3];

    string line;
    size_t iline = 0;
    size_t iframe = 0;
    while (getline(file, line)) {
        if (line[0] == '#') {
            if (line.substr(2, 6) == "Chunk ") {
            //# Chunk Coord1 Coord2 Coord3 Ncount vx vy vz
                //trim(line);
                //vector<string> tokens;
                //split(tokens, line, is_any_of(" "), token_compress_on);
                //datasets.resize(tokens.size() - 6);
                //std::copy(tokens.begin() + 6, tokens.end(), datasets.begin());
            }
            continue;
        }
        trim(line);
        vector<string> tokens;
        split(tokens, line, is_any_of(" "), token_compress_on);
        if (tokens.size() == 3) { // new frame
            ++iframe;
            iline = 0;
            continue;
        }
        if (iframe-1 != targetFrame)
            continue;

        assert(tokens.size() > 5);
        std::vector<double> values(tokens.size());    
        std::transform(tokens.begin(), tokens.end(), values.begin(), [](const string& val) {return stod(val);});

        if (iline == 0) {
            leftBottom[0] = values[1] - h[0]/2.0;
            leftBottom[1] = values[2] - h[1]/2.0;
            leftBottom[2] = values[3] - h[2]/2.0;
            
            double l[] = {-2.0*(leftBottom[0]), -2.0*(leftBottom[1]), -2.0*(leftBottom[2])};            
            size_t n[3];
            for (int i = 0; i < 3; ++i) {
                n[i] = static_cast<size_t>(l[i]/h[i]);
            }
            std::cout << "Grid size: " << n[0] << " " <<  n[1] << " " << n[2] << endl;
            std::cout << "Box size: " << l[0] << " " <<  l[1] << " " << l[2] << endl;
            grid.resize(n[0], n[1], n[2]);
            grid.setBoundingBox(Box3D(l[0], l[1], l[2]));
        }
        double point[] = {values[1] - leftBottom[0], values[2] - leftBottom[1], values[3] - leftBottom[2]};
        size_t index[3];
        computeIndex(grid, point, index);
       
        // x 
        grid(index[0], index[1], index[2]) = values[5];

        ++iline;
    }
    file.close();
    std::cout << "Only 0th frame converter. N Frames: " << iframe << std::endl;
}

void printHeaderStat(const VGrid& grid)
{
    auto bb = grid.getBoundingBox();
    for (int i = 0; i < 3; ++i)
        std::cout << "Box size in " << i << "-direction: " <<bb.getIthSize(i) << ", grid sz:" << grid.size(i) << std::endl;
}

int main(int argc, char ** argv)
{
    if (argc != 7)
    {
        printf("usage: ./dat2hdf5 <input-file> <output-file> <frame> <hx> <hy> <hz>");
        return 1;
    }
    double h[] = {stod(argv[4]), stod(argv[5]), stod(argv[6])};
    
    std::cout << "Parsing input file..." << std::endl;
    VGrid grid;
    parseLammpsAveSpatial(stoi(argv[3]), argv[1], h, grid, "data");

    printHeaderStat(grid);    
    
    std::cout << "Saving hdf5..." << std::endl;
    VSerializerHDF5 writer(grid, argv[2], "data");
    writer.run();
    std::cout << "Convertation finished!" << std::endl;
    return 0;
}

