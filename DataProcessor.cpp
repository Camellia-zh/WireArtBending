#include "DataProcessor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

DataProcessor::DataProcessor(const std::string& inputFilePath, const std::string& outputFilePath)
    : inputFilePath_(inputFilePath), outputFilePath_(outputFilePath) {}

std::vector<Point> DataProcessor::read_input_file() {
    std::vector<Point> inputVertices;
    std::ifstream fin(inputFilePath_);
    int numOfPoints = -1;
    fin >> numOfPoints;
    inputVertices.resize(numOfPoints);
    for (int i = 0; i < numOfPoints; i++) {
        fin >> inputVertices[i].x >> inputVertices[i].y;
    }
    fin.close();
    return inputVertices;
}

std::vector<double> DataProcessor::calculate_distances(const std::vector<Point>& points) {
    std::vector<double> distances;
    for (size_t i = 0; i < points.size() - 1; ++i) {
        double dx = points[i + 1].x - points[i].x;
        double dy = points[i + 1].y - points[i].y;
        double distance = round(std::sqrt(dx * dx + dy * dy));
        distances.push_back(distance);
    }
    return distances;
}

std::vector<double> DataProcessor::calculate_angles(const std::vector<Point>& points) {
    std::vector<double> angles;
    for (size_t i = 0; i < points.size() - 2; ++i) {
        double v1_x = points[i].x - points[i + 1].x;
        double v1_y = points[i].y - points[i + 1].y;
        double v2_x = points[i + 1].x - points[i + 2].x;
        double v2_y = points[i + 1].y - points[i + 2].y;

        double cross = v1_x * v2_y - v1_y * v2_x;
        double dot = v1_x * v2_x + v1_y * v2_y;
        double norm_v1 = std::sqrt(v1_x * v1_x + v1_y * v1_y);
        double norm_v2 = std::sqrt(v2_x * v2_x + v2_y * v2_y);
        double angle = round(std::acos(dot / (norm_v1 * norm_v2)) * 180.0 / 3.1415926535);

        if (cross < 0) {
            angle = -angle;
        }

        angles.push_back(angle);
    }

    return angles;
}

void DataProcessor::write_to_file(const std::vector<double>& distances, const std::vector<double>& angles) {
    std::ofstream outfile(outputFilePath_);
    outfile << "Distances: ";
    for (const auto& d : distances) {
        outfile << d << " ";
    }
    outfile << std::endl;
    outfile << "Angles: ";
    for (const auto& a : angles) {
        outfile << a << " ";
    }
    outfile << std::endl;
    outfile.close();
}

void DataProcessor::display_file_content(const std::string& filePath) {
    std::ifstream infile(filePath);
    std::string line;
    while (std::getline(infile, line)) {
        std::cout << line << std::endl;
    }
    infile.close();
}

void DataProcessor::process() {
    std::vector<Point> inputVertices = read_input_file();
    std::vector<double> distances = calculate_distances(inputVertices);
    std::vector<double> angles = calculate_angles(inputVertices);
    write_to_file(distances, angles);
    display_file_content(outputFilePath_);
}

void DataProcessor::process2(std::vector<double>& distances, std::vector<double>& angles, int& anglesnum) {
    std::ifstream infile(outputFilePath_);
    if (!infile) {
        std::cout << "无法打开文件" << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string keyword;
        ss >> keyword;

        if (keyword == "Distances:") {
            double num;
            while (ss >> num) {
                distances.push_back(num);
            }
        }
        else if (keyword == "Angles:") {
            double num;
            while (ss >> num) {
                angles.push_back(num);
                anglesnum++;
            }
        }
    }
    infile.close();

}