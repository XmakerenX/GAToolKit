//-------------------------------------------------------------------------------------------------
// GeneticsAlgorithem.cpp : This code is trying to find the best solution for the knapsack problem
// using genetic algorithm, The qulity of the answer is checked againts the optimal solution
//-------------------------------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
#include <limits>
#include <numeric>
#include <functional>
#include <list>  
#include <vector>
     
#include <cmath>
#include <boost/dynamic_bitset.hpp>
#include <boost/random/discrete_distribution.hpp>

typedef boost::dynamic_bitset<> DynamicBitSet;

//-----------------------------------------------------------------------------------------------
// Structures
//-----------------------------------------------------------------------------------------------
struct Item
{
    Item(int _val, int _weight)
    {
        val = _val;
        weight = _weight;
    }
    
    int val;
    int weight;
};

struct Pop
{
    Pop(DynamicBitSet _chromosome, int _fitness)
    {
        chromosome = _chromosome;
        fitness = _fitness;
    }
    
    bool operator== (const Pop& b)
    {
        return chromosome == b.chromosome;
    }
    
    DynamicBitSet chromosome;
    int fitness;
};

//-----------------------------------------------------------------------------------------------
// Functions
//-----------------------------------------------------------------------------------------------
int                         knapSackDynamic         (int W, const std::vector<Item> items);
std::vector<Item>           generateKnapsackItems   (int numberOfItems = 0);
std::vector<Pop>            generateFirstGeneration (int numberOfItems, const std::vector<Item>& items, int totalWeight);
void                        generatePopulation      (int popToGenerate, const std::vector<Item>& items, int totalWeight, std::vector<Pop>& population);
void                        pickOnlyTheFittests     (std::vector<Pop>& population, const std::vector<Item>& items, int totalWeight);
void                        createNewGeneration     (std::vector<Pop>& population);

void                        mateParents             (int parent1Index, int parent2Index, std::vector<Pop>& population, std::vector<Pop>& newPopulation);
void                        mutateChild             (DynamicBitSet& child);
int                         calcFitness             (DynamicBitSet a, const std::vector<Item>& items, int totalWeight);
bool                        cmp                     (const Pop& a, const Pop& b);
bool                        compareInterval         (const Pop& a, const Pop& b);

int                         AddFitness              (const int& left, const Pop& right);

const int populationNumber = 1000;
const int eliteSize = 100;

std::default_random_engine generator;

//-----------------------------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------------------------
int main()
{    
    typedef std::chrono::high_resolution_clock myclock;
    
    myclock::time_point beginning = myclock::now();
    generator.seed(beginning.time_since_epoch().count());
    
    std::uniform_int_distribution<int> weightDist(50, 300);
    
    std::vector<Item> items = generateKnapsackItems();
    int W = weightDist(generator);
    
    std::cout << "Number of items " << items.size() << "\n";
    std::cout << "Knapsack Weight " << W << "\n";
    
    std::cout << "the items are \n";
    for (const Item& item : items)
        std::cout << item.val << " " << item.weight << "\n";
    
    //Create first generation
    std::vector<Pop> population = generateFirstGeneration(items.size(), items, W);
              
    int lastFitnessAvg = -300000;
    int nGenOfNoChange = 0;
    int generation;
    
    //Apply GA for 1000 generations
    for (generation = 0; generation < 1000; generation++)
    {
        createNewGeneration(population);
        pickOnlyTheFittests(population, items, W);
        
        std::cout << "Gen " << generation << "'s fittest is:" << population[0].fitness << std::endl;
        std::cout << "Gen " << generation << "'s fittest chromosome is:" << population[0].chromosome << std::endl;
        
        for (int i = 0; i < populationNumber; i++) 
        {
            population[i].fitness = calcFitness(population[i].chromosome, items, W);
        }
        
        int fitnessAvg = std::accumulate(population.begin(), population.end(), 0, AddFitness) / populationNumber;
        std::cout << "fitness averge " << fitnessAvg << "\n";
        
//         if (lastFitnessAvg == fitnessAvg)
//         {
//             nGenOfNoChange++;
//             if (nGenOfNoChange == 10)
//             {
//                 // nuke the population but 1 from orbit and try again with new population
//                 population.erase(population.begin() + 1, population.end());
//                 generatePopulation(9, items, W,population);
//             }
//                 
//         }
//         else
//         {
//             nGenOfNoChange = 0;
//             lastFitnessAvg = fitnessAvg;
//         }
        
        
    }
    
    // make sure fitness was updated for the last gen
    for (Pop& p : population)
        p.fitness = calcFitness(p.chromosome, items, W);
            
    std::sort(population.begin(),population.end(), compareInterval);
    // print the last gen
    std::cout << "last generation was" << generation << std::endl;
    std::cout << "fittest is:" << population[0].fitness << std::endl;
    std::cout << "fittest chromosweightome is:" << population[0].chromosome << std::endl;
    std::cout << "dynamic knapsack solution : " << knapSackDynamic(W, items) << "\n";
    // print the knapsack problem that was solved
    std::cout << "Number of items " << items.size() << "\n";
    std::cout << "Knapsack Weight " << W << "\n";
    int weight  = 0;
    int i = population[0].chromosome.find_first();
    
    if (i == DynamicBitSet::npos)
    {
        return -20000;
    }
    
    weight += items[i].weight;
    
    while ( (i = population[0].chromosome.find_next(i)) != DynamicBitSet::npos)
    {
        weight += items[i].weight;
    }
    std::cout << "weight :" << weight << "\n";
    
    std::cout << "the items are \n";
    for (const Item& item : items)
        std::cout << item.val << " " << item.weight << "\n";
    
    return 0;
}

//-----------------------------------------------------------------------------------------------
// generateKnapsackItems
//-----------------------------------------------------------------------------------------------
std::vector<Item> generateKnapsackItems(int numberOfItems)
{
    std::uniform_int_distribution<int> itemNumDist(50, 64);
    std::uniform_int_distribution<int> itemDist(1, 100);
    
    std::vector<Item> items;
    
    if (numberOfItems == 0)
        numberOfItems = itemNumDist(generator);
    
    for (int i = 0; i < numberOfItems; i++)
        items.emplace_back(itemDist(generator), itemDist(generator) );
    
    return items;
}

//-----------------------------------------------------------------------------------------------
// knapSackDynamic
// Returns the maximum value that can be put in a knapsack of capacity W
//-----------------------------------------------------------------------------------------------
int knapSackDynamic(int W, const std::vector<Item> items)
{
    int i, w;
#ifdef _WIN32
	int ** K;
	K = new int*[items.size() + 1];
	for (int i = 0; i < items.size() + 1; i++)
		K[i] = new int[W + 1];
#else
	int K[items.size() + 1][W + 1];
#endif

    // Build table K[][] in bottom up manner
    for (i = 0; i <= items.size(); i++)
    {
        for (w = 0; w <= W; w++)
        {
            if (i == 0 || w == 0)
                K[i][w] = 0;
            else if (items[i - 1].weight <= w)
                    K[i][w] = std::max(items[i - 1].val + K[i - 1][w - items[i - 1].weight], K[i - 1][w]);
                 else
                    K[i][w] = K[i - 1][w];

        }
    }

	int ret = K[items.size()][W];
#ifdef _WIN32
	// free K 
	for (int i = 0; i < items.size() + 1; i++)
		delete[] K[i];
	delete[] K;
#endif
    return ret;
}

//-----------------------------------------------------------------------------------------------
// generateFirstGeneration
//-----------------------------------------------------------------------------------------------
std::vector<Pop> generateFirstGeneration(int numberOfItems, const std::vector<Item>& items, int totalWeight)
{
    unsigned long maxNum = std::pow(2, numberOfItems);
    std::uniform_int_distribution<long> longDistributaion(0, maxNum);
    int fitness;
    
    // the genes pool
    std::vector<Pop> population;
    
    for (int i = 0; i < populationNumber; i++)
    {
        DynamicBitSet chromosome = DynamicBitSet(items.size(), longDistributaion(generator));
        fitness = calcFitness(chromosome, items, totalWeight);
            
        //each solution is in the range of 0 to maxNum
        population.push_back( Pop(chromosome, fitness));                
    }
    
    return population;
}

//-----------------------------------------------------------------------------------------------
// generatePopulation
//-----------------------------------------------------------------------------------------------
void generatePopulation(int popToGenerate, const std::vector<Item>& items, int totalWeight, std::vector<Pop>& population)
{
    unsigned long maxNum = std::pow(2, items.size());
    std::uniform_int_distribution<long> longDistributaion(0, maxNum);
    int fitness;
    
    for (int i = 0; i < popToGenerate; i++)
    {
        DynamicBitSet chromosome = DynamicBitSet(items.size(), longDistributaion(generator));
        fitness = calcFitness(chromosome, items, totalWeight);
            
        //each solution is in the range of 0 to maxNum
        population.push_back( Pop(chromosome, fitness));                
    }
}

//-----------------------------------------------------------------------------------------------
// pickOnlyTheFittests
// input:      population
// algorithem: picks half of the population that has the best scores and "kills the rest"
// output:     the output array has the best s
//-----------------------------------------------------------------------------------------------
void pickOnlyTheFittests(std::vector<Pop>& population, const std::vector<Item>& items, int totalWeight)
{    
    for (Pop& p : population)
    {
        p.fitness = calcFitness(p.chromosome, items, totalWeight);
    }
            
    std::sort(population.begin(),population.end(), compareInterval);
    population.erase(population.begin() + populationNumber, population.end());
}
//-----------------------------------------------------------------------------------------------
// cmp
// compare from higest to lowest scores
//-----------------------------------------------------------------------------------------------
bool cmp(const Pop& a, const Pop& b)
{
    return (a.fitness < b.fitness);
}

//-----------------------------------------------------------------------------------------------
// compareInterval
// compare from lowest to higest scores
//-----------------------------------------------------------------------------------------------
bool compareInterval(const Pop& a, const Pop& b)
{
    return (a.fitness > b.fitness);
}

//-----------------------------------------------------------------------------------------------
// createNewGeneration
//-----------------------------------------------------------------------------------------------
void createNewGeneration(std::vector< Pop >& population)
{
    std::vector<Pop> newPopulation;
    
    int min = 10000;
    for (int i = 0; i < population.size(); i++)
    {
        if (population[i].fitness < min)
            min = population[i].fitness;
    }
    
    if (min < 0)
    {
        for (int i = 0; i < population.size(); i++)
            population[i].fitness -= (min + 1);
    }
    
    int fitnessSum = std::accumulate(population.begin(), population.end(), 0 ,AddFitness);
    
    std::vector<double> probabilities;
    for (int i = 0; i < population.size(); i++)
    {
        probabilities.push_back((double)population[i].fitness / fitnessSum);
    }
    
    boost::random::discrete_distribution<> dist(probabilities);
    for (int i = 0; i < populationNumber; i += 2)
    {
        int parent1Index, parent2Index;
        parent1Index = dist(generator);
        parent2Index = dist(generator);
        newPopulation.push_back(population[parent1Index]);
        newPopulation.push_back(population[parent2Index]);
        mateParents(parent1Index, parent2Index, population, newPopulation);
    }

    if (min < 0)
    {
        for (int i = 0; i < population.size(); i++)
            population[i].fitness += (min + 1);
    }
    
     std::sort(population.begin(),population.end(), compareInterval);
     for (int i = 0; i < 1; i++)
         newPopulation.push_back(population[i]);
        
    population = std::move(newPopulation);
}

//-----------------------------------------------------------------------------------------------
// mateParents
//-----------------------------------------------------------------------------------------------
void mateParents(int parent1Index, int parent2Index, std::vector<Pop>& population, std::vector<Pop>& newPopulation)
{
    DynamicBitSet& parent1 = population[parent1Index].chromosome;
    DynamicBitSet& parent2 = population[parent2Index].chromosome;
    
    DynamicBitSet child, child2;

    std::uniform_int_distribution<int> joinBitDist(1, parent1.size() - 1);
    //int joinBitLocation = (rand() % parent1.size());
    int joinBitLocation = joinBitDist(generator);
    joinBitLocation = parent1.size() - joinBitLocation;
        
    // mask starts as numberOfBits ones
    unsigned long temp = pow(2, parent1.size() + 1) - 1;
    DynamicBitSet mask = DynamicBitSet(parent1.size(), temp);
    // leave ones in the mask only from joinBitLocation to LSB
    mask = mask >> joinBitLocation;
    DynamicBitSet flippedMask = mask;
    flippedMask.flip();
            
    DynamicBitSet chromozomeX = parent1 & flippedMask; // upper  part
    DynamicBitSet chromozomeY = parent2 & mask;          // bottom part
    child = chromozomeX | chromozomeY;
    mutateChild(child);
        
    chromozomeX = parent2 & flippedMask; // upper  part
    chromozomeY = parent1 & mask;          // bottom part
    child2 = chromozomeX | chromozomeY;
    mutateChild(child2);
                
     newPopulation.push_back(Pop(child, 0));
     newPopulation.push_back(Pop(child2, 0));
}

//-----------------------------------------------------------------------------------------------
// mutateChild
//-----------------------------------------------------------------------------------------------
void mutateChild(DynamicBitSet& child)
{
    std::uniform_int_distribution<int> mutateDistributaion(0, 6);
    
    if (mutateDistributaion(generator) == 0)
    {
        std::uniform_int_distribution<int> bitDistributaion(0, child.size() - 1);
        // what bit to flip? when 0 is the LSB and 11 the MSB
        int randMutation = bitDistributaion(generator);
        child.flip(randMutation);
    }
}

//-----------------------------------------------------------------------------------------------
// calcFitness
//-----------------------------------------------------------------------------------------------
int calcFitness(DynamicBitSet a, const std::vector<Item>& items, int totalWeight)
{
    int fitness = 0;
    int weight = 0;
    int i = a.find_first();
    
    if (i == DynamicBitSet::npos)
    {
        return 0;
    }
    
    fitness += items[i].val;
    weight += items[i].weight;
    
    while ( (i = a.find_next(i)) != DynamicBitSet::npos)
    {
        fitness += items[i].val;
        weight += items[i].weight;
    }
    // penalized over the weight limit solutions
    if (weight > totalWeight)
    {
        fitness = totalWeight - weight;
    }
    
    return fitness;
}

//-----------------------------------------------------------------------------------------------
// AddFitness
// used for the std::accumulate
//-----------------------------------------------------------------------------------------------
int AddFitness(const int& left, const Pop& right)
{
    return left + right.fitness;
}
