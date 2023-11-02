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
	//������ͽǶ�
	std::vector<double> distances;
	std::vector<double> angles, angles2, angles3;
	int anglesnum = 0;
	if (argc < 2) {
		std::cout << "���ṩsample������" << std::endl;
		return 1;
	}
	std::string sample = argv[1];
	std::string inputFilePath = "..\\file\\bending_input\\" + sample + "\\" + sample + ".txt";
	std::string outputFilePath = "..\\file\\bending_input\\" + sample + "\\distances_Angles.txt";
	//std::string inputFilePath = "E:\\WireArtPjt\\input\\bending_input\\" + sample + "\\" + sample + ".txt";
	//std::string outputFilePath = "E:\\WireArtPjt\\input\\bending_input\\" + sample + "\\distances_Angles.txt";
	DataProcessor processor(inputFilePath, outputFilePath);
	//processor.process();    //����� �������ͽǶ� ������ļ�����ʾ
	processor.process2(distances, angles, anglesnum);

	////��ʼ��
	//Wire wire;
	//wire.initialize();
	////wire.modifyPoint(9, 100.0, -50.0, 75.0);
	//Eigen::MatrixXd V1, V2, V3;
	//Eigen::MatrixXi F1, F2, F3;
	//igl::readPLY("..\\file\\PLY\\cube.PLY", V1, F1);
	//igl::readPLY("..\\file\\PLY\\rack_pin.PLY", V2, F2);
	//igl::readPLY("..\\file\\PLY\\bender_gear.PLY", V3, F3);
	////����
	//angles2 = angles; angles3 = angles;
	//Bending bending(V1, V2, V3, F1, F2, F3, wire.points, distances, angles, anglesnum, angles2);
	//bending.drawAnimation();
	//// ��� angles ��ֵ
	//std::cout << "angles: ";
	//for (const double& angle : angles) {
	//	std::cout << angle << " ";
	//}
	//std::cout << std::endl;
	//// ��� angles2 ��ֵ
	//std::cout << "angles2: ";
	//for (const double& angle : angles2) {
	//	std::cout << angle << " ";
	//}
	//std::cout << std::endl;
	//Bending bending2(V1, V2, V3, F1, F2, F3, wire.points, distances, angles2, anglesnum, angles3);
	//bending2.drawAnimation();


	//�ҵ��Ƭ�� ȥ�� graph cut �Ŵ��㷨
	/*
	//�ҵ��Ƭ�� ȥ��
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
	//smooth term�����������ʱ�������޸�
	for (int l1 = 0; l1 < num_labels; l1++)
		for (int l2 = 0; l2 < num_labels; l2++)
			if (l1 == l2)//�þ��󣬶Խ����ϵľ�����Ϊ0����������Ϊ1
				smooth[l1 + l2 * num_labels] = 0;
			else
				smooth[l1 + l2 * num_labels] = 1;
	//data term��������Ҫ����
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
			data[i * num_labels + j] = pixels_data_terms[i][j];//���¾���
			//std::cout << data[i * num_labels + j] << std::endl;
		}
	try
	{
		//input necessary term to graph-cut ����Ҫ��term���뵽graph-cut��
		GCoptimizationGeneralGraph* gc = new GCoptimizationGeneralGraph(num_pixels, num_labels);
		gc->setDataCost(data);
		gc->setSmoothCost(smooth);
		//Set the weight between two neighboring pixel ���������ڵ�֮������Ȩ��
		//400 is another parameter which can be tuned according to other considerations,
		// 400����һ�����Ը�������������Ĳ���
		//which indicate the weight comparison between the data term and the weight between two neighboring pixel
		//����ʾdata term���weight comparison�Լ��������ڶ���֮���Ȩ��
		//������һ�����������ڵ�ı�
		for (int i = 0; i < num_pixels-1; i++) {
			gc->setNeighbors(i, i+1, 1);
		}
		//optimization process ��󻯹���
		std::cout << "Before optimization energy is " << gc->compute_energy() << std::endl;
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		std::cout << "After optimization energy is " << gc->compute_energy() << std::endl;
		
		//get results ��ý��
		for (int i = 0; i < num_pixels; i++)
		{
			graphcut_index[i] = gc->whatLabel(i);//����һ��int�͵�labelֵ
			std::cout << graphcut_index[i] << " ";
		}
		std::cout << std::endl;
		//cout << "������Խ���" << endl;
	}
	catch (GCException e) {
		e.Report();
	}

	std::vector<int> angle_indices;//�ı�Ľǵ�������0��ʼ
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
	std::vector<int> genesArray; // ���ݸ��Ŵ��㷨�Ļ�������  �ı����ֵ
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
	geneAlgorithm.Run(5); // �����Ŵ��㷨�����ݻ�������͵�������
	*/



	//��ײ������
    //Collision_check collisionChecker;
	//bool a=collisionChecker.checkCollision(distances, angles, 0, 4);//3����
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
	//		/*outputFile << "Ԫ��: ((" << std::get<0>(element).first << ", " << std::get<0>(element).second << "), ("
	//			<< std::get<1>(element).first << ", " << std::get<1>(element).second << "), ("
	//			<< std::get<2>(element).first << ", " << std::get<2>(element).second << "))" << std::endl;*/
	//		outputFile << std::get<1>(element).second << " "
	//			<< std::get<2>(element).first << "  " << std::get<2>(element).second << std::endl;
	//	}
	//	// �ر��ļ�
	//	outputFile.close();
	//	std::cout << "������д���ļ� 'output.txt'" << std::endl;
	//}
	//else {
	//	std::cerr << "�޷����ļ� 'output.txt' ��д������" << std::endl;
	//}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "������ʱ��: " << duration.count() << " ��" << std::endl;

	return 0;
}


