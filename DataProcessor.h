#ifndef DATAPROCESSOR_H
#define DATAPROCESSOR_H

#include <string>
#include <vector>

struct Point {
    double x;
    double y;
};

class DataProcessor {
public:
    DataProcessor(const std::string& inputFilePath, const std::string& outputFilePath);
    void process();
    void process2(std::vector<double>& distances, std::vector<double>& angles, int& anglesnum);

private:
    std::string inputFilePath_;
    std::string outputFilePath_;

    std::vector<Point> read_input_file();  //¶ÁµãµÄÎ»ÖÃ
    std::vector<double> calculate_distances(const std::vector<Point>& points);
    std::vector<double> calculate_angles(const std::vector<Point>& points);
    void write_to_file(const std::vector<double>& distances, const std::vector<double>& angles);
    void display_file_content(const std::string& filePath);

};

#endif // DATAPROCESSOR_H

