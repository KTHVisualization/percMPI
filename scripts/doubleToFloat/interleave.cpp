#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    int numComps = argc-4;
    int numBytesHeader = atoi(argv[1]);
    int numBytesFooter = atoi(argv[2]);
    argv+=3;

std::cout << "Omitting header of " << numBytesHeader << " bytes\n";
    std::ofstream outfile(argv[numComps], std::ios::binary);
    std::ifstream* inputs = new std::ifstream[numComps];
    for (int c = 0; c < numComps; ++c){
        inputs[c].open(argv[c], std::ios::binary);
        inputs[c].seekg(numBytesHeader, std::ios::beg);
std::cout << "Opening file " << argv[c] << "\n";
    }

    inputs[0].seekg(0, std::ios::end);
    int numBytes = (int)inputs[0].tellg() - numBytesFooter;
    inputs[0].seekg(numBytesHeader, std::ios::beg);
    
std::cout << "Reading " << numBytes/sizeof(double) << " doubles\n";

    for (int i = 0; i < numBytes / sizeof(double); ++i) {
        for (int c = 0; c < numComps; ++c) {
            double val;
            inputs[c].read((char*)(&val), sizeof(double));
            outfile << (float)val;
        }
    }
}