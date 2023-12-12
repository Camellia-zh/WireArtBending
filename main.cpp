#include "DataProcessor.h"
#include "WireBending.h"
#include <iostream>
#include <igl/readPLY.h>
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <utility>
#include "GCoptimization.h" 
#include "GeneAlg.h"
#include "BeamSearch.h"
#include <chrono>
#include "build/fitting.h"
#include "lisan.cpp"
#include <CGAL/intersections.h>


std::chrono::duration<double> timeTry = std::chrono::duration<double>::zero();
std::chrono::duration<double> timeCC1 = std::chrono::duration<double>::zero();
std::chrono::duration<double> timeCC2 = std::chrono::duration<double>::zero();////cro碰撞检测时间

int main(int argc, char* argv[])
{
	//auto start_input = std::chrono::high_resolution_clock::now();
	////读距离和角度
	//std::vector<double> distances;
	//std::vector<double> angles, angles2, angles3;
	//int anglesnum = 0;
	//if (argc < 2) {
	//	std::cout << "请提供sample参数！" << std::endl;
	//	return 1;
	//}
	//std::string sample = argv[1];
	//std::string inputFilePath = "..\\file\\bending_input\\" + sample + "\\" + sample + ".txt";
	//std::string outputFilePath = "..\\file\\bending_input\\" + sample + "\\distances_Angles.txt";
	////std::string inputFilePath = "E:\\WireArtPjt\\input\\bending_input\\" + sample + "\\" + sample + ".txt";
	////std::string outputFilePath = "E:\\WireArtPjt\\input\\bending_input\\" + sample + "\\distances_Angles.txt";
	//DataProcessor processor(inputFilePath, outputFilePath);
	////processor.process();    //读入点 计算距离和角度 输出到文件并显示
	//processor.process2(distances, angles, anglesnum);
	//auto end_input = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> duration_input = end_input - start_input;
	//std::cout << "读数据经过的时间: " << duration_input.count() << " 秒" << std::endl;

	auto start = std::chrono::high_resolution_clock::now();
	std::vector<ShapeElement> ShapeElements;
	std::vector<double> shapeAngles;
	
	// Create instances of FittingElement for straight lines
	ShapeElement a1("segment", 27);
	ShapeElement a2("segment", 33);
	ShapeElement a4("segment", 31);
	ShapeElement a5("segment", 22);
	ShapeElement a7("segment", 30);
	ShapeElement a8("segment", 27);
	ShapeElement a9("segment", 33);
	ShapeElement a10("segment", 31);
	// Create instances of CurvePoint and FittingElement for arcs
	ShapeElement a3("arc", 63, 0.02);
	ShapeElement a6("arc", 1, 0.01);
	ShapeElement a11("arc", 38, 0.01);


	// Add elements to the vector
	ShapeElements.push_back(a1);
	ShapeElements.push_back(a2);
	ShapeElements.push_back(a3);
	ShapeElements.push_back(a4);
	//ShapeElements.push_back(a5);
	//ShapeElements.push_back(a6);
	//ShapeElements.push_back(a7);
	//ShapeElements.push_back(a8);
	//ShapeElements.push_back(a9);
	//ShapeElements.push_back(a10);
	//ShapeElements.push_back(a11);


	shapeAngles.push_back(100.0);
	shapeAngles.push_back(65.0);
	shapeAngles.push_back(40);
/*	
	shapeAngles.push_back(60);
	shapeAngles.push_back(100);
	shapeAngles.push_back(95);
	shapeAngles.push_back(-60);
	shapeAngles.push_back(23);
	shapeAngles.push_back(23);
	shapeAngles.push_back(45)*/;




	//找到最长片段 去重 
	/*{
		auto start = std::chrono::high_resolution_clock::now();
		//找到最长片段 去重
		std::vector<std::pair<int, int>> allSegments;
		for (int i = 0; i <= anglesnum; ++i) {
			if (i != anglesnum)
			{
				std::pair<int, int> result = findLongestNonCollidingSegment(distances, angles, i);
				allSegments.push_back(result);
				std::cout << i;
			}
			else {
				std::pair<int, int> result = findLongestNonCollidingSegment2(distances, angles, i);
				allSegments.push_back(result);
				std::cout << i;
			}
		}
		std::set<std::pair<int, int>> uniqueSegments(allSegments.begin(), allSegments.end());
		std::vector<double> collisionAngles;
		std::vector<double> collisionAnglesSeq(angles.size(), 0.0);
		std::cout << std::endl;
		std::cout << "Unique longest non-colliding segments:" << std::endl;
		for (const auto& segment : uniqueSegments) {
			std::cout << "L" << segment.first + 1 << "---- L" << segment.second + 1 << std::endl;
			std::cout << "A" << segment.first << "---- A" << segment.second + 1 << std::endl;
			collisionAngles.push_back(segment.first - 1);//0开始
			collisionAngles.push_back(segment.second);
		}

		//权重不同
		// 排序 collisionAngles
		std::sort(collisionAngles.begin(), collisionAngles.end());
		double firstValue = *collisionAngles.begin();
		double lastValue = *std::prev(collisionAngles.end());
		// 去除第一个值
		if (!collisionAngles.empty()) {
			collisionAngles.erase(collisionAngles.begin());
		}
		for (auto it = collisionAngles.begin(); it != collisionAngles.end() && *it == firstValue;) {
			it = collisionAngles.erase(it);
		}
		// 去除最后一个值
		if (!collisionAngles.empty()) {
			collisionAngles.erase(std::prev(collisionAngles.end()));
		}
		while (!collisionAngles.empty() && collisionAngles.back() == lastValue) {
			collisionAngles.erase(std::prev(collisionAngles.end()));
		}

		////权重都为1
		//// 排序 collisionAngles
		//std::sort(collisionAngles.begin(), collisionAngles.end());
		//// 去除重复元素
		//collisionAngles.erase(std::unique(collisionAngles.begin(), collisionAngles.end()), collisionAngles.end());
		//// 去除第一个值
		//if (!collisionAngles.empty()) {
		//	collisionAngles.erase(collisionAngles.begin());
		//}
		//// 去除最后一个值
		//if (!collisionAngles.empty()) {
		//	collisionAngles.pop_back();
		//}

		//输出 collisionAngles 的值
		std::cout << "collisionAngles: ";
		for (const double& point : collisionAngles) {
			std::cout << point << " ";
		}
		for (const double point : collisionAngles) {
			if (point >= 0 && point < collisionAnglesSeq.size()) {
				collisionAnglesSeq[static_cast<size_t>(point)] += 1.0;
			}
		}
		//输出 collisionAnglesSeq 的值
		std::cout << std::endl << "collisionAnglesSeq: ";
		for (const double& point : collisionAnglesSeq) {
			std::cout << point << " ";
		}
		std::cout << std::endl;


		double max_angle = *std::max_element(collisionAnglesSeq.begin(), collisionAnglesSeq.end());
		if (max_angle != 0.0) {
			for (double& angle : collisionAnglesSeq) {
				angle /= max_angle;
			}
		}
		std::cout << std::endl << "collisionAnglesSeq: ";
		for (const double& point : collisionAnglesSeq) {
			std::cout << point << " ";
		}
		std::cout << std::endl;


		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> duration = end - start;
		std::cout << "找到最长片段经过的时间: " << duration.count() << " 秒" << std::endl;
	}*/

	//beam search without weight
	std::cout << "anglewithoutweight" << std::endl;
	BeamSearch beamSearch(ShapeElements, shapeAngles);

	std::vector<std::vector<double>> results = beamSearch.runSearch(10);

	std::cout << "Beam Search Results:" << std::endl;
	for (const std::vector<double>& result : results) {
		for (double value : result) {
			std::cout << value << " ";
		}
		std::cout << std::endl;
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "beam search经过的时间: " << duration.count() << " 秒" << std::endl;
	return 0;
}


