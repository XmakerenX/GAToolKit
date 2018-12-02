// GeneticsAlgorithem.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
#include <limits>
#include <functional>
#include <list>  
#include <vector>
     
#include <cmath>
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<> DynamicBitSet;

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
    
    DynamicBitSet chromosome;
    int fitness;
};

int                         knapSackDynamic         (int W, const std::vector<Item> items);
std::vector<Item>           generateKnapsackItems   (int numberOfItems = 0);
std::vector<Pop>            generateFirstGeneration (int numberOfItems, const std::vector<Item>& items, int totalWeight);
void                        pickOnlyTheFittests     (std::vector<Pop>& population, const std::vector<Item>& items, int totalWeight);
void                        createNewGeneration     (std::vector<Pop>& population);

void                        mateParents             (int parent1Index, int parent2Index, std::vector< Pop >& population);
void                        mutateChild             (DynamicBitSet& child);
int                         calcFitness             (DynamicBitSet a, const std::vector<Item>& items, int totalWeight);
bool                        cmp                     (DynamicBitSet a, DynamicBitSet b);
bool                        compareInterval         (const Pop& a, const Pop& b);

//This code is trying to find the best solution
//to the following function: f(x) = x^2 where: 0 << x <<4096
//obviously, the solution is 4096 and the problem is easy, but its interesting to see how
//A genetic algorithem will behave, and that is why I am writing this code

const int populationNumber = 1000;
const int eliteSize = 100;
// const int numberOfBits = 12;
// const int searchSpace = pow(2,numberOfBits);

std::default_random_engine generator;

int main()
{    
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
//     
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
    
    int lastFittest = -300000;
    int con = 0;
    //Apply GA for 50 generations
    for (int gen = 0; gen < 1000; gen++)
    {
        pickOnlyTheFittests(population, items, W);
        std::cout << "Gen " << gen << "'s fittest is:" << population[0].fitness << std::endl;
        std::cout << "Gen " << gen << "'s fittest chromosome is:" << population[0].chromosome << std::endl;
        
//         if (lastFittest == population[0].fitness)
//         {
//             con++;
//             if (con == 10)
//                 break;
//         }
//         else
//         {
//             con = 0;
//             lastFittest = population[0].fitness;
//         }
        
        createNewGeneration(population);
        for (int i = 0; i < populationNumber; i++) 
            std::cout << population[i].fitness << std::endl;
    }
    
    for (Pop& p : population)
        p.fitness = calcFitness(p.chromosome, items, W);
            
    std::sort(population.begin(),population.end(), compareInterval);
    
    std::cout << "fittest is:" << population[0].fitness << std::endl;
    std::cout << "fittest chromosome is:" << population[0].chromosome << std::endl;
    std::cout << "dynamic knapsack solution : " << knapSackDynamic(W, items) << "\n";
    
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

std::vector<Item> generateKnapsackItems(int numberOfItems)
{
    std::uniform_int_distribution<int> itemNumDist(10, 64);
    std::uniform_int_distribution<int> itemDist(1, 100);
    
    std::vector<Item> items;
    
    if (numberOfItems == 0)
        numberOfItems = itemNumDist(generator);
    
    for (int i = 0; i < numberOfItems; i++)
        items.emplace_back(itemDist(generator), itemDist(generator) );
    
    return items;
}

// Returns the maximum value that can be put in a knapsack of capacity W
int knapSackDynamic(int W, const std::vector<Item> items)
{
    int i, w;
    int K[items.size() + 1][W + 1];
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

    return K[items.size()][W];
}

//--------------------------------------------------------------------------
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
        while  (calcFitness(chromosome, items, totalWeight) == 0)
        {
            chromosome = DynamicBitSet(items.size(), longDistributaion(generator));
            fitness = calcFitness(chromosome, items, totalWeight);
        }
            
        //each solution is in the range of 0 to maxNum
        population.push_back( Pop(chromosome, fitness));                
    }
    
    return population;
}
//--------------------------------------------------------------------------
//input:      population
//algorithem: picks half of the population that has the best scores and "kills the rest"
//output:     the output array has the best s
void pickOnlyTheFittests(std::vector<Pop>& population, const std::vector<Item>& items, int totalWeight)
{    
    for (Pop& p : population)
    {
        p.fitness = calcFitness(p.chromosome, items, totalWeight);
    }
    
    //sort:
    for (int i = 0; i < population.size(); i++) 
        std::cout << population[i].fitness << std::endl;
        
    std::sort(population.begin(),population.end(), compareInterval);

    // eraes half which are the least fittest
    //population.erase(population.begin() + (population.size() / 2) , population.end());
    population.erase(population.end() - eliteSize , population.end());
    
    for (int i = 0; i < eliteSize; i++)
        population.push_back(population[i]);
}
//--------------------------------------------------------------------------
//compare from higest to lowest scores
bool cmp(DynamicBitSet a, DynamicBitSet b)
{
	DynamicBitSet a_score = a&a;
	DynamicBitSet b_score = b&b;
	return (a_score >= b_score);
}
//--------------------------------------------------------------------------
//compare from lowest to higest scores
bool compareInterval(const Pop& a, const Pop& b)
{
	return (a.fitness > b.fitness);
}
//-------------------------------------------------------------------------
void createNewGeneration(std::vector< Pop >& population)
{
	std::random_shuffle(population.begin() + 2, population.end());

	//fill the "killed" individuals with offsprings of fittest parents
	// create 2 children from each pair
        int popSize = population.size();
        for (int i = eliteSize; i < popSize; i += 2)
        {
            mateParents(i, i + 1, population);
        }
}
//-------------------------------------------------------------------------
void mateParents(int parent1Index, int parent2Index, std::vector< Pop >& population)
{
    DynamicBitSet& parent1 = population[parent1Index].chromosome;
    DynamicBitSet& parent2 = population[parent2Index].chromosome;
    
    int joinBitLocation = (rand() % parent1.size());
    joinBitLocation = parent1.size() - joinBitLocation;
    //joinBitLocation = parent1.size() / 2;
    
    // mask starts as numberOfBits ones
    unsigned long temp = pow(2, parent1.size() + 1) - 1;
    DynamicBitSet mask = DynamicBitSet(parent1.size(), temp);
    // leave ones in the mask only from joinBitLocation to LSB
    mask = mask >> joinBitLocation;
        
    DynamicBitSet chromozomeX = parent1 & (mask.flip()); // upper  part
    DynamicBitSet chromozomeY = parent2 & mask;          // bottom part
    DynamicBitSet child = chromozomeX | chromozomeY;
    mutateChild(child);
    parent1 = child;
    //population.push_back(Pop(child, 0));
    
    chromozomeX = parent2 & (mask.flip()); // upper  part
    chromozomeY = parent1 & mask;          // bottom part
    child = chromozomeX | chromozomeY;
    mutateChild(child);
    parent2 = child;
    //population.push_back(Pop(child, 0));
}

void mutateChild(DynamicBitSet& child)
{
    std::uniform_int_distribution<int> mutateDistributaion(0, 6);
    
    if (mutateDistributaion(generator) == 0)
    {
        std::cout << child.size() << "\n";
        std::uniform_int_distribution<int> bitDistributaion(0, child.size() - 1);
        // what bit to flip? when 0 is the LSB and 11 the MSB
        int randMutation = bitDistributaion(generator);
        child.flip(randMutation);
    }
}

int calcFitness(DynamicBitSet a, const std::vector<Item>& items, int totalWeight)
{
    int fitness = 0;
    int weight = 0;
    int i = a.find_first();
    
    if (i == DynamicBitSet::npos)
    {
        return -20000;
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
        //fitness -= (100 + (totalWeight - weight));
        fitness = -10000 + (totalWeight - weight);
        //fitness = 0;
    
    return fitness;
}
