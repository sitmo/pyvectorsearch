#include <iostream>
#include <fstream>
#include <sstream>
#include <sstream>
#include <string>
#include <array>
#include <vector>

#include "../include/pktree.hpp"

struct CityInfo {
    float latlon[2];
   std::string name;
};

void read_cities(const char* filename, std::vector<CityInfo>& cities)
{
    std::string line, city_name, latitude_str, longitude_str;
    float latitude, longitude;

    std::fstream file(filename, std::ios::in);
	if(file.is_open()) {
	    getline(file, line); // skip the first header row

	    // process data rows
		while(getline(file, line)) {

			std::stringstream str(line);

			// get the first 3 columns
			getline(str, city_name, ',');
			getline(str, latitude_str, ',');
			getline(str, longitude_str, ',');

			latitude = std::stof(latitude_str);
			longitude = std::stof(longitude_str);

            cities.push_back(CityInfo({latitude, longitude, city_name}));
		}
	}
}


int main(int argc, const char *argv[])
{

    // Spatial search variables
    typedef pkspatial::pkmap<float*, std::string> MapType;
	MapType::TResultsetKnn knn_search_results;
	MapType::TResultsetRange range_search_results;

    MapType map(2); // a map with 2D coordinates

    // Read the CSV file with raw city into into a vector of CityInfo structs.
    std::vector<CityInfo> cities;

    if (argc > 1)
        read_cities(argv[1], cities);
    else
        read_cities("../example_cpp/dutch_cities.csv", cities);


    // Insert all city data to the map
    for (auto i=cities.begin(); i!=cities.end(); ++i)
            map.insert(i->latlon, i->name);

    // Search key
	float search_key[2] = {52.0117, 4.3592};


	// Do a KNN search
	map.search_knn(search_key, 5, knn_search_results, true);

    std::cout << std::endl << "knn=5 search results:" << std::endl;
	for (auto i = knn_search_results.begin(); i != knn_search_results.end(); ++i) {

	    double distance = std::sqrt(i->first);
	    float* latlon = i->second.key;
	    std::string city_name = i->second.value;

	    std::cout << distance << ": (" << latlon[0] << ", " << latlon[1] << ") "<< city_name << std::endl;
	}

    //  Do a RANGE search
	map.search_range(search_key, 0.1, range_search_results);

    std::cout << std::endl << "range < 0.1 search results:" << std::endl;

	for (auto i = range_search_results.begin(); i != range_search_results.end(); ++i) {

	    float* latlon = i->key;
	    std::string city_name = i->value;

	    std::cout <<  "(" << latlon[0] << ", " << latlon[1] << ") " << city_name << std::endl;
	}
    return 0;
}