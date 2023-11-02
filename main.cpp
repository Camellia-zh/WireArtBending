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

int main(int argc, char* argv[])
{
	auto start = std::chrono::high_resolution_clock::now();
	//读距离和角度
	std::vector<double> distances;
	std::vector<double> angles, angles2, angles3;
	int anglesnum = 0;
	if (argc < 2) {
		std::cout << "请提供sample参数！" << std::endl;
		return 1;
	}
	std::string sample = argv[1];
	std::string inputFilePath = "..\\file\\bending_input\\" + sample + "\\" + sample + ".txt";
	std::string outputFilePath = "..\\file\\bending_input\\" + sample + "\\distances_Angles.txt";
	//std::string inputFilePath = "E:\\WireArtPjt\\input\\bending_input\\" + sample + "\\" + sample + ".txt";
	//std::string outputFilePath = "E:\\WireArtPjt\\input\\bending_input\\" + sample + "\\distances_Angles.txt";
	DataProcessor processor(inputFilePath, outputFilePath);
	//processor.process();    //读入点 计算距离和角度 输出到文件并显示
	processor.process2(distances, angles, anglesnum);

	////初始化
	//Wire wire;
	//wire.initialize();
	////wire.modifyPoint(9, 100.0, -50.0, 75.0);
	//Eigen::MatrixXd V1, V2, V3;
	//Eigen::MatrixXi F1, F2, F3;
	//igl::readPLY("..\\file\\PLY\\cube.PLY", V1, F1);
	//igl::readPLY("..\\file\\PLY\\rack_pin.PLY", V2, F2);
	//igl::readPLY("..\\file\\PLY\\bender_gear.PLY", V3, F3);
	////动画
	//angles2 = angles; angles3 = angles;
	//Bending bending(V1, V2, V3, F1, F2, F3, wire.points, distances, angles, anglesnum, angles2);
	//bending.drawAnimation();
	//// 输出 angles 的值
	//std::cout << "angles: ";
	//for (const double& angle : angles) {
	//	std::cout << angle << " ";
	//}
	//std::cout << std::endl;
	//// 输出 angles2 的值
	//std::cout << "angles2: ";
	//for (const double& angle : angles2) {
	//	std::cout << angle << " ";
	//}
	//std::cout << std::endl;
	//Bending bending2(V1, V2, V3, F1, F2, F3, wire.points, distances, angles2, anglesnum, angles3);
	//bending2.drawAnimation();


	//找到最长片段 去重 graph cut 遗传算法
	/*
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
		}
	}
	std::set<std::pair<int, int>> uniqueSegments(allSegments.begin(), allSegments.end());
	std::cout << "Unique longest non-colliding segments:" << std::endl;
	for (const auto& segment : uniqueSegments) {
		std::cout << "L" << segment.first + 1 << "---- L" << segment.second + 1 << std::endl;
	}

	//Graph Cut 
	std::size_t num_labels = uniqueSegments.size();
	int num_pixels = anglesnum + 1;
	int* smooth = new int[num_labels * num_labels]; //smooth term
	int* data = new int[num_pixels * num_labels];	  //data term
	//smooth term，这个设置暂时不进行修改
	for (int l1 = 0; l1 < num_labels; l1++)
		for (int l2 = 0; l2 < num_labels; l2++)
			if (l1 == l2)//该矩阵，对角线上的均设置为0，其余设置为1
				smooth[l1 + l2 * num_labels] = 0;
			else
				smooth[l1 + l2 * num_labels] = 1;
	//data term，这里需要设置
	std::vector<std::vector<int>> pixels_data_terms(num_pixels, std::vector<int>(num_labels, 10000.0));
	std::set<std::pair<int, int>>::iterator it = uniqueSegments.begin();
	std::vector<int> graphcut_index(num_pixels);
	for (int i = 0; i < num_labels; i++)
	{
		std::pair<int, int> firstElement = *it;
		int firstNumber = firstElement.first;
		int secondNumber = firstElement.second;
		for (int j = firstNumber; j <= secondNumber; j++)
		{
			pixels_data_terms[j][i] = 0;
		}
		std::advance(it, 1);
	}
	for (int i = 0; i < num_pixels; i++)
		for (int j = 0; j < num_labels; j++)
		{
			data[i * num_labels + j] = pixels_data_terms[i][j];//更新矩阵
			//std::cout << data[i * num_labels + j] << std::endl;
		}
	try
	{
		//input necessary term to graph-cut 将必要的term输入到graph-cut中
		GCoptimizationGeneralGraph* gc = new GCoptimizationGeneralGraph(num_pixels, num_labels);
		gc->setDataCost(data);
		gc->setSmoothCost(smooth);
		//Set the weight between two neighboring pixel 在两个相邻点之间设置权重
		//400 is another parameter which can be tuned according to other considerations,
		// 400是另一个可以根据情况来调整的参数
		//which indicate the weight comparison between the data term and the weight between two neighboring pixel
		//它表示data term间的weight comparison以及两个相邻顶点之间的权重
		//先设置一个轮廓上相邻点的边
		for (int i = 0; i < num_pixels-1; i++) {
			gc->setNeighbors(i, i+1, 1);
		}
		//optimization process 最大化过程
		std::cout << "Before optimization energy is " << gc->compute_energy() << std::endl;
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		std::cout << "After optimization energy is " << gc->compute_energy() << std::endl;
		
		//get results 获得结果
		for (int i = 0; i < num_pixels; i++)
		{
			graphcut_index[i] = gc->whatLabel(i);//返回一个int型的label值
			std::cout << graphcut_index[i] << " ";
		}
		std::cout << std::endl;
		//cout << "这里可以进来" << endl;
	}
	catch (GCException e) {
		e.Report();
	}

	std::vector<int> angle_indices;//改变的角的索引（0开始
	int prev_value = graphcut_index[0];
	for (int i = 1; i < graphcut_index.size(); i++) {
		if (graphcut_index[i] != prev_value) {
			angle_indices.push_back(i);
			prev_value = graphcut_index[i];
		}
	}
	for (int& value : angle_indices) {
		value--;
	}
	std::vector<int> genesArray; // 传递给遗传算法的基因数组  改变的数值
	for (int index : angle_indices) {
		genesArray.push_back(angles[index]);
	}
	for (int i : angle_indices) {
		std::cout << i << " ";
	}
	std::cout << std::endl;
	for (int angle : genesArray) {
		std::cout << angle << " ";
	}
	std::cout << std::endl;

	GeneAlg geneAlgorithm(distances,angles, angle_indices, genesArray);
	geneAlgorithm.Run(5); // 运行遗传算法，传递基因数组和迭代次数
	*/



	//碰撞检测测试
    //Collision_check collisionChecker;
	//bool a=collisionChecker.checkCollision(distances, angles, 0, 4);//3个角
	//std::cout << a;

    //beam search
    BeamSearch beamSearch(distances, angles);

	//beam search
	std::vector<std::vector<double>> results = beamSearch.runSearch(5);
	std::cout << "Beam Search Results:" << std::endl;
	for (const std::vector<double>& result : results) {
		for (double value : result) {
			std::cout << value << " ";
		}
		std::cout << std::endl;
	}

	////-150~150
	//std::vector<double> bendingsequence = beamSearch.getBendingSequence();
	//for (const double& value : bendingsequence) {
	//	std::cout << value << " ";
	//}
	//std::cout << std::endl;
	
	////space
	//std::vector<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<int, int>>> bendingrange= beamSearch.getBendingSpace(0, 1);
	//std::ofstream outputFile("output.txt");
	//if (outputFile.is_open()) {
	//	for (const auto& element : bendingrange) {
	//		/*outputFile << "元素: ((" << std::get<0>(element).first << ", " << std::get<0>(element).second << "), ("
	//			<< std::get<1>(element).first << ", " << std::get<1>(element).second << "), ("
	//			<< std::get<2>(element).first << ", " << std::get<2>(element).second << "))" << std::endl;*/
	//		outputFile << std::get<1>(element).second << " "
	//			<< std::get<2>(element).first << "  " << std::get<2>(element).second << std::endl;
	//	}
	//	// 关闭文件
	//	outputFile.close();
	//	std::cout << "数据已写入文件 'output.txt'" << std::endl;
	//}
	//else {
	//	std::cerr << "无法打开文件 'output.txt' 以写入数据" << std::endl;
	//}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "经过的时间: " << duration.count() << " 秒" << std::endl;

	return 0;
}


