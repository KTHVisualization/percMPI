#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>

using ind = std::int64_t;

void makeVTI()
{
	std::ifstream inStream("../Data/iso512x512x512.raw", std::ios::in | std::ios::binary);
	std::ifstream wrapStream("../Data/header.combined", std::ios::in | std::ios::binary);
	std::ofstream outStream("../Data/iso512x512x512.vti", std::ios::out | std::ios::binary);


	// Write 
	std::string line;
	while(getline(wrapStream, line))
	{
		outStream << line << '\n';
	}
	outStream << "   _";

	// Get data size.
	inStream.seekg(0, std::ios::end);
	ind size = inStream.tellg();
	char* data = new char[size];
	inStream.seekg(0, std::ios::beg);

	// Write file size.
	std::int64_t size64 = size;
	std::cout << size64 << std::endl;
	outStream.write(reinterpret_cast<char*>(&size64), sizeof(std::int64_t));

	// Write data.
	inStream.read(data, size);
	outStream.write(data, size);

	// Write end.
	outStream << "\n  </AppendedData>\n</VTKFile>\n";

	delete[] data;
	inStream.close();
	wrapStream.close();
	outStream.close();
}

void resize512()
{
	const ind MULTIPLIER = 8;
	std::ifstream inStream("../Data/Iso4096/iso512x512x512.raw", std::ios::in | std::ios::binary);
	std::ofstream outStream("../Data/Iso4096/iso4096.raw", std::ios::out | std::ios::binary);
	
	inStream.seekg(0, std::ios::end);
	ind size = inStream.tellg();
	std::cout << "Filesize: " << size << std::endl;
	std::cout << "Expected: " << 512*512*512*sizeof(float) << std::endl; 
	assert(size == 512 * 512 * 512 * sizeof(float) && "Incorrect input size!");
	char* dataRaw = new char[size];
	inStream.seekg(0, std::ios::beg);

	// Read data.
	inStream.read(dataRaw, size);
	float* data = reinterpret_cast<float*>(dataRaw);

	// Out data.
	ind sideLength = (512 - 1) * MULTIPLIER + 1;

	ind totalSize = 0;
	float* result = new float[sideLength * sideLength * MULTIPLIER];
	//std::vector<float> result;
	//result.reserve(sideLength * sideLength * sideLength);

	for (ind z = 0; z < 512; ++z)
	{
		//for (ind test = 0; test < sideLength * sideLength * MULTIPLIER; ++test)
		//	result[test] = (test % 2) ? 0 : 100;

		for (ind y = 0; y < 512; ++y)
			for (ind x = 0; x < 512; ++x)
			{
				ind index = ind(x) + 512 * y + ind(512 * 512 * z);

				// Get all the neighbors we need for inderpolation.
				float corners[8];
				for (ind zB = 0; zB < 2; ++zB)
					for (ind yB = 0; yB < 2; ++yB)
						for (ind xB = 0; xB < 2; ++xB)
						{
							if (x + xB == 512 || y + yB == 512 || z + zB == 512) continue;

							ind cornerID = ind(x) + xB + ind(y + yB) * 512 + ind(z + zB) * 512 * 512;
							corners[xB + yB * 2 + zB * 4] = data[cornerID];
						}

				for (ind xOff = 0; xOff < MULTIPLIER; ++xOff)
					for (ind yOff = 0; yOff < MULTIPLIER; ++yOff)
						for (ind zOff = 0; zOff < MULTIPLIER; ++zOff)
						{
							//if (xOff == 0 && yOff == 0 && zOff == 0) continue;
							if ((x * MULTIPLIER + xOff >= sideLength) ||
								(y * MULTIPLIER + yOff >= sideLength) ||
								(z * MULTIPLIER + zOff >= sideLength)) continue;

							float xF = float(xOff) / MULTIPLIER;
							float yF = float(yOff) / MULTIPLIER;
							float zF = float(zOff) / MULTIPLIER;
							float xI = 1.0f - xF;
							float yI = 1.0f - yF;
							float zI = 1.0f - zF;

							float value =
								zI * yI * zI * corners[0] +
								zI * yI * zF * corners[1] +
								zI * yF * zI * corners[2] +
								zI * yF * zF * corners[3] +

								zF * yI * zI * corners[4] +
								zF * yI * zF * corners[5] +
								zF * yF * zI * corners[6] +
								zF * yF * zF * corners[7];

							ind indexBig =
								ind(x * MULTIPLIER + xOff) +
								ind(y * MULTIPLIER + yOff) * sideLength +
								ind(zOff) * sideLength * sideLength;
							result[indexBig] = value;


						}
			}
		ind numElements = ind(sideLength) * sideLength  * sizeof(float);
		if (z < 511) numElements *= MULTIPLIER;

		outStream.write(reinterpret_cast<char*>(result), numElements);
		totalSize += numElements;
	}


	/*
	for (ind z = 0; z < sideLength; ++z)
		for (ind y = 0; y < sideLength; ++y)
			for (ind x = 0; x < sideLength; ++x)
			{
				ind index = x + sideLength * y + sideLength * sideLength * z;
				ind b
			}
			*/
	delete[] result;
	inStream.close();
	outStream.close();
	std::cout << "Done\n";
	getchar();
}

int main()
{
	//makeVTI();
	resize512();
}
