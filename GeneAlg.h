#ifndef GENEALG_H
#define GENEALG_H

#include <iostream>
#include <ctime>
#include <vector>
#include <list>
#include <algorithm> 

class GeneAlg
{
public:
	GeneAlg(std::vector<double> distances, std::vector<double> angles, std::vector<int> angle_indices, std::vector<int>& genesArray);
	void Run(int numGenerations);

	std::vector<double> distances;
	std::vector<double> angles;
	std::vector<int> angle_indices;
	std::vector<int> genesArray;
	std::vector<int> genesArrayFir;
	int population;
	int genes;
	int top;
	int geneval;

	struct Individual
	{
		std::vector<int> chromosome;//int chromosome[genes];
		int fitness;

		Individual(int genes) {
			chromosome.resize(genes);
		}
	};

	std::vector<Individual> individuals;
	std::vector<Individual> topIndividuals;

	int Random(int start, int end);
	//static bool Comp(const Individual& a, const Individual& b);
	static bool Comp(const Individual& a, const Individual& b) {
		return a.fitness > b.fitness;
	}
	int GetFitness(Individual* in);
	Individual Single_Point_Crossover(Individual* a, Individual* b);
	Individual Two_Point_Crossover(Individual* a, Individual* b);
	void Initialize();
	void Selection();
	void Crossover();
	//void Mutation();
	void Reset();
	void Print();

private:

};


#endif