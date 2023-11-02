#include "GeneAlg.h"
#include "WireBending.h"

GeneAlg::GeneAlg(std::vector<double> distances, std::vector<double> angles, std::vector<int> angle_indices, std::vector<int>& genesArray)
    :distances(distances), angles(angles), angle_indices(angle_indices), genesArray(genesArray), genesArrayFir(genesArray),
    population(30),
    genes(genesArray.size()),
    top(10),
    geneval(300)
{
    srand(static_cast<unsigned int>(time(0)));
    //srand(time(0));
    Initialize();
    Print();
}


int GeneAlg::Random(int start, int end) {
    int range = (end - start) + 1;
    int random_int = start + (rand() % range);
    return random_int;
}

//static bool GeneAlg::Comp(const Individual& a, const Individual& b) {
//    return a.fitness > b.fitness;
//}

int GeneAlg::GetFitness(Individual* in) {
    int change=1;
    int collision = 0;
    std::vector<double> newAngle=angles;
    for (int i = 0; i < angle_indices.size(); i++) {
        newAngle[angle_indices[i]] = in->chromosome[i];
    }
    for (int i = 0; i < genes; i++) {
        if (in->chromosome[i] != genesArrayFir[i]) {
            change++;
        }
    }
    Collision_check collisionChecker;
    if (collisionChecker.checkCollision(distances, angles, 0, angles.size()))
    {
        collision = -50;
    }

    int a = genes/change  + collision;
    in->fitness = a;
    return a;
}

GeneAlg::Individual GeneAlg::Single_Point_Crossover(Individual* a, Individual* b) {
    int r = Random(genes / 2, genes - 1);
    Individual child(genes);
    for (int i = 0; i < genes; i++) {
        child.chromosome[i] = a->chromosome[i];
    }
    for (int j = r; j < genes; j++) {
        child.chromosome[j] = b->chromosome[j];
    }
    GetFitness(&child);
    return child;
}

GeneAlg::Individual GeneAlg::Two_Point_Crossover(Individual* a, Individual* b) {
    int r = Random(0, genes / 2);
    int endpoint = Random(genes / 2, genes - 1);
    Individual child(genes);
    for (int i = 0; i < genes; i++) {
        child.chromosome[i] = a->chromosome[i];
    }
    for (int j = r; j < endpoint; j++) {
        child.chromosome[j] = b->chromosome[j];
    }
    GetFitness(&child);
    return child;
}

void GeneAlg::Initialize() {
    for (int start = 0; start < population; start++) {
        Individual ind(genes);
        ind.chromosome = genesArray;
        int a = Random(0, population - 1);
        int b = Random(0, genes - 1);
        int c = Random(0, geneval);
        while (ind.chromosome[b] == (c - (geneval / 2))) {
            c = Random(0, geneval);
        }
        ind.chromosome[b] = c - (geneval / 2);

        GetFitness(&ind);
        individuals.push_back(ind);
        std::sort(individuals.begin(), individuals.end(), Comp);
    }
}

void GeneAlg::Selection() {
    int i = 0;
    for (auto it : individuals) {
        if (i == top)
            break;
        topIndividuals.push_back(it);
        i++;
    }
}

void GeneAlg::Crossover() {
    individuals.clear();
    for (int start = 0; start < population; start++) {
        Individual ind(genes);
        int x = Random(0, top - 1);
        int y = Random(0, top - 1);
        int rand = Random(0, 1);
        switch (rand) {
        case 0:
            ind = Single_Point_Crossover(&topIndividuals[x], &topIndividuals[y]);
            break;
        case 1:
            ind = Two_Point_Crossover(&topIndividuals[x], &topIndividuals[y]);
            break;
        default:
            break;
        }
        int m = Random(0, 30);
        if (m == 0) {
            int a = Random(0, population - 1);
            int b = Random(0, genes - 1);
            int c = Random(0, geneval);
            while (ind.chromosome[b] == (c - (geneval / 2))) {
                c = Random(0, geneval);
            }
            ind.chromosome[b] = c - (geneval / 2);
        }
        GetFitness(&ind);
        individuals.push_back(ind);
        std::sort(individuals.begin(), individuals.end(), Comp);
    }
}

void GeneAlg::Reset() {
    topIndividuals.clear();
}

void GeneAlg::Print() {
    int num = 0;
    std::vector<Individual>::iterator iter = individuals.begin();
    for (iter; iter != individuals.end(); iter++) {
        for (int i = 0; i < genes; i++) {
            std::cout << iter->chromosome[i] << "-";
        }
        printf(" [ %d ]  %d \t", iter->fitness, num);
        num++;
        std::cout << std::endl;
    }
}

void GeneAlg::Run(int num) {
    int age = 0;
    int clear = 0;
    //int m = 0;
    std::cout << " 遗传算法的代数" << age << std::endl;
    //std::cin >> m;
    while (true) {
        for (int i = 0; i < num; i++) {
            clear = 0;
            system("cls");
            Reset();
            Selection();
            Crossover();
            Print();
            for (auto it : individuals) {
                if (it.fitness > 50) {
                    clear++;
                }
                if (clear == population) {
                    age++;
                    std::cout << age << "代完成全部大于50" << std::endl;
                    return;
                }
            }
            age++;
            std::cout << age << "代" << std::endl;
        }
        return;
    }
}
