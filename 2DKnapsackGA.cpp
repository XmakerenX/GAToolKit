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
struct LegoItem
{
    LegoItem(int _width, int _height)
        :width(_width), height(_height)
    {}
    
    int width;
    int height;
};

struct Knapsack2DCase
{
    Knapsack2DCase(int _width, int _height, const std::vector<LegoItem>& _items)
        :width(_width), height(_height), items(_items)
    {
        coordinateBits = std::ceil(std::log2(std::max(width, height)));
    }
    
    Knapsack2DCase(int _width, int _height, const std::vector<LegoItem>&& _items)
        :width(_width), height(_height), items(_items)
    {
        coordinateBits = std::ceil(std::log2(std::max(width, height)));
    }
    
    int width;
    int height;
    long unsigned int coordinateBits;
    std::vector<LegoItem> items;
};

struct Rect
{
    Rect(int _left, int _top, int _right, int _bottom)
        :left(_left), top(_top), right(_right), bottom(_bottom)
    {}
    
    static Rect intersect(const Rect& a,const Rect& b)
    {
        return Rect(std::max(a.left, b.left), std::max(a.top, b.top), std::min(a.right, b.right), std::min(a.bottom, b.bottom));
    }
    
    int getWidth()
    {
        return right - left;
    }
    
    int getHeight()
    {
        return bottom - top;
    }
    
    int left;
    int top;
    int right;
    int bottom;
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

//-----------------------------------------------------------------------------------------------
// Functions
//-----------------------------------------------------------------------------------------------
std::vector<LegoItem>       generateKnapsackItems   (int width, int height, int numberOfItems = 0);
DynamicBitSet               generateChromosome      (long unsigned int totalBitsNum);
void                        repairChromosome        (const Knapsack2DCase& knapsackCase, DynamicBitSet& chromosome);
std::vector<Pop>            generateFirstGeneration (const Knapsack2DCase& knapsackCase);
void                        generatePopulation      (int popToGenerate, const Knapsack2DCase& knapsackCase, std::vector<Pop>& population);
void                        pickOnlyTheFittests     (const Knapsack2DCase knapsackCase,std::vector< Pop >& population);
void                        createNewGeneration     (const Knapsack2DCase& knapsackCase, std::vector<Pop>& population);
void                        mateParents             (int parent1Index, int parent2Index, std::vector< Pop >& population, std::vector<Pop>& newPopulation);
void                        mutateChild             (DynamicBitSet& child);
int                         calcFitness             (const Knapsack2DCase& knapsackCase, const DynamicBitSet& chromosome);
bool                        cmp                     (DynamicBitSet a, DynamicBitSet b);
bool                        compareInterval         (const Pop& a, const Pop& b);

int                         AddFitness              (const int& left, const Pop& right);
void                        PrintSolution           (const Knapsack2DCase& knapsackCase,const DynamicBitSet& chromosome);

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
    
    int paperWidth = 100;
    int paperHeight = 100;
    int sum = 0;
    
    Knapsack2DCase knapsackCase(paperWidth, paperHeight, generateKnapsackItems(paperWidth, paperHeight));
    
    std::cout << "Knapsack dimentions are " << paperWidth << "X" << paperHeight << "\n";
    std::cout << "Number of items " << knapsackCase.items.size() << "\n";
    std::cout << "the items are \n";
    for (const LegoItem& item : knapsackCase.items)
    {
        sum += item.width * item.height;
        std::cout << item.width << "X" << item.height << " " << item.width * item.height <<  "\n";
    }
    std::cout << "total items area " << sum << "\n";
    
    //Create first generation
    std::vector<Pop> population = generateFirstGeneration(knapsackCase);
    
    int lastFitnessAvg = 1050000;
    int nGenOfNoChange = 0;
    int generation;
    //Apply GA for 1000 generations
    for (generation = 0; generation < 1000; generation++)
    {
        createNewGeneration(knapsackCase, population);
        pickOnlyTheFittests(knapsackCase, population);
        
        std::cout << "Gen " << generation << "'s fittest is:" << population[0].fitness << std::endl;
        std::cout << "Gen " << generation << "'s fittest chromosome is:" << population[0].chromosome << std::endl;
        PrintSolution(knapsackCase, population[0].chromosome);
        
        for (int i = 0; i < populationNumber; i++) 
        {
            population[i].fitness = calcFitness(knapsackCase, population[i].chromosome);
            //std::cout << population[i].fitness << std::endl;
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
//                 generatePopulation(999, knapsackCase ,population);
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
        p.fitness = calcFitness(knapsackCase, p.chromosome);
            
    std::sort(population.begin(),population.end(), compareInterval);
    PrintSolution(knapsackCase ,population[0].chromosome);
    
    // print the last gen
    std::cout << "last generation was " << generation << std::endl;
    std::cout << "fittest is: " << population[0].fitness << std::endl;
    std::cout << "fittest chromosweightome is: " << population[0].chromosome << std::endl;
    // print the knapsack problem that was solved
    std::cout << "Knapsack dimentions are " << paperWidth << "X" << paperHeight << "\n";
    std::cout << "Number of items " << knapsackCase.items.size() << "\n";
    std::cout << "the items are \n";
    for (const LegoItem& item : knapsackCase.items)
        std::cout << item.width << "X" << item.height  << " " << item.width * item.height <<  "\n";
    std::cout << "total items area " << sum << "\n";
        
    return 0;
}

//-----------------------------------------------------------------------------------------------
// generateKnapsackItems
//-----------------------------------------------------------------------------------------------
std::vector<LegoItem> generateKnapsackItems(int width, int height, int numberOfItems)
{
    int minDimention = std::min(width,height);
    std::uniform_int_distribution<int> itemSizeDist(1, (minDimention / 5)*2);
    std::uniform_int_distribution<int> itemNumDist(20, 30);
    
    std::vector<LegoItem> legoItems;
    
    if (numberOfItems == 0)
        numberOfItems = itemNumDist(generator);
    
    for (int i = 0; i < numberOfItems; i++)
        legoItems.emplace_back(itemSizeDist(generator), itemSizeDist(generator));
    
    return legoItems;
}

//-----------------------------------------------------------------------------------------------
// generateChromosome
// generates a random chromosome with size of totalBitsNum bits
//-----------------------------------------------------------------------------------------------
DynamicBitSet generateChromosome(long unsigned int totalBitsNum)
{
    // get max value of long unsigned int
    long unsigned int  maxNum = 0;
    maxNum--;
    // generate values between 0 and the max value for long unsigned int
    std::uniform_int_distribution<long unsigned int> longDistributaion(0, maxNum);
    
    std::vector<long unsigned int> chromosomeInput;
    long unsigned int curTotalBitsNum = totalBitsNum;
    
    // DynamicBitSet bigger than 64bits need to be set by a vector of long unsigned int
    // so randomly generate such a vector with the last element being the reminder bits left
    while (curTotalBitsNum > 0)
    {
        if (curTotalBitsNum > 64)
        {
            chromosomeInput.push_back(longDistributaion(generator));
            curTotalBitsNum -= 64;
        }
        else
        {
            long unsigned int reminderMaxNum = std::pow(2, curTotalBitsNum);
            reminderMaxNum--;
            std::uniform_int_distribution<long unsigned int> reminderDistributaion(0, reminderMaxNum);
            chromosomeInput.push_back(reminderDistributaion(generator));
            curTotalBitsNum -= curTotalBitsNum;
        }
    }
    
    DynamicBitSet chromosome = DynamicBitSet(chromosomeInput.begin(), chromosomeInput.end());
    // the chromosome size will be in multiplication of 64 so shirnk it to totalBitsNum
    // to avoid usesless zeros at the binary string start
    chromosome.resize(totalBitsNum);
    
    return chromosome;
}

//-----------------------------------------------------------------------------------------------
// repairChromosome
// Move items that are in an invalid X,Y coordinates
// The Chromosomes are being randomly generated which can lead to X,Y positions for items
// that will put them out of the container so this function fixes thier coordinates
//-----------------------------------------------------------------------------------------------
void repairChromosome(const Knapsack2DCase& knapsackCase, DynamicBitSet& chromosome)
{
    // get how many bits are used to hold the coordinate number
    const long unsigned int& coordinateBits = knapsackCase.coordinateBits;
    long unsigned int bitsPerItem = 2 + 2 * coordinateBits;
    
    for (int i = 0; i < knapsackCase.items.size(); i++)
    {
        // craete mask to get current item info from the chromosome
        unsigned int itemMaskValue = std::pow(2, bitsPerItem);
        itemMaskValue--;
        DynamicBitSet itemMask =  DynamicBitSet(chromosome.size(), itemMaskValue);
        itemMask = itemMask << bitsPerItem * i;
        // get the item info from the chromosome
        DynamicBitSet itemValues = chromosome & itemMask;
        itemValues = itemValues >> bitsPerItem * i;
        // check if the item is in selected to be in the container
        if (itemValues[1 + 2 * coordinateBits] == 1)
        {
            // create mask to get the y value from the item bit string
            unsigned int yMaskValue = std::pow(2, coordinateBits);
            yMaskValue--;
            DynamicBitSet yMask = DynamicBitSet(chromosome.size(), yMaskValue);
            // get the item y coordinate
            unsigned int y = (itemValues & yMask).to_ulong();
            // create mask to get the x value from the item bit string
            DynamicBitSet xMask = yMask << coordinateBits;
            DynamicBitSet xValue = (itemValues & xMask) >> coordinateBits;
            // get the item x coordinate
            unsigned int x = xValue.to_ulong();
            unsigned int itemWidth = knapsackCase.items[i].width;
            unsigned int itemHeight = knapsackCase.items[i].height;
            // swap item width and height if item is vertical orientation
            if (itemValues[2 * coordinateBits] == 1)
            {
                std::swap(itemWidth, itemHeight);
            }
            // repair if needed the x and y coordinates
            if (x + itemWidth > knapsackCase.width)
                x = knapsackCase.width - itemWidth;
            if (y + itemHeight > knapsackCase.height)
                y = knapsackCase.height - itemHeight;
            
            x = x << coordinateBits;
            DynamicBitSet updatedXYValue = DynamicBitSet(chromosome.size(), x | y);
            // copy the updated value to chromosome
            for (int j = 0; j < (bitsPerItem - 2); j++)
            {
                chromosome[bitsPerItem * i + j] = updatedXYValue[j];
            }
        }
    }
}

//-----------------------------------------------------------------------------------------------
// generateFirstGeneration
//-----------------------------------------------------------------------------------------------
std::vector<Pop> generateFirstGeneration(const Knapsack2DCase& knapsackCase)
{
    long unsigned int coordinateBits = std::ceil(std::log2(std::max(knapsackCase.width, knapsackCase.height)));
    // one bit tell if item is in or out , second bit tell if item is horizontal or vertical
    // and the other bits are for the X and Y coordinate of the item
    long unsigned int totalBitsNum = knapsackCase.items.size()*2 + 2 * knapsackCase.items.size() * coordinateBits;
    
    int fitness;
    
    // the genes pool
    std::vector<Pop> population;
    
    for (int i = 0; i < populationNumber; i++)
    {
        DynamicBitSet chromosome = generateChromosome(totalBitsNum);
        // the random chromosome can have invalid values and might need to be repaired
        repairChromosome(knapsackCase, chromosome);
        fitness = calcFitness(knapsackCase, chromosome);
        population.push_back( Pop(chromosome, fitness));                
    }
    
    return population;
}

//-----------------------------------------------------------------------------------------------
// generatePopulation
// The same as generateFirstGeneration but generate population based popToGenerate
// used if the population become undiverse to attempt to move the algorithm from the current 
// optimium to a new one and better one
// ** This function should just replace generateFirstGeneration
//-----------------------------------------------------------------------------------------------
void generatePopulation(int popToGenerate, const Knapsack2DCase& knapsackCase, std::vector<Pop>& population)
{
    long unsigned int coordinateBits = std::ceil(std::log2(std::max(knapsackCase.width, knapsackCase.height)));
    // one bit tell if item is in or out , second bit tell if item is horizontal or vertical
    // and the other bits are for the X and Y coordinate of the item
    long unsigned int totalBitsNum = knapsackCase.items.size()*2 + 2 * knapsackCase.items.size() * coordinateBits;
    
    for (int i = 0; i < popToGenerate; i++)
    {
        DynamicBitSet chromosome = generateChromosome(totalBitsNum);
        repairChromosome(knapsackCase, chromosome);
        int fitness = calcFitness(knapsackCase, chromosome);
        population.push_back( Pop(chromosome, fitness));                
    }
}

//-----------------------------------------------------------------------------------------------
// pickOnlyTheFittests
// input:      population
// algorithem: clamps the population back to <populationNumber> and does so based on pop fitness
// output:     The <populationNumber> most fittest solutions 
//-----------------------------------------------------------------------------------------------
void pickOnlyTheFittests(const Knapsack2DCase knapsackCase,std::vector< Pop >& population)
{    
    for (Pop& p : population)
    {
        p.fitness = calcFitness(knapsackCase, p.chromosome);
    }
            
    std::sort(population.begin(),population.end(), compareInterval);
    population.erase(population.begin() + populationNumber, population.end());
}
//-----------------------------------------------------------------------------------------------
// cmp
// compare from higest to lowest scores
// !!!! unused !!!!
//-----------------------------------------------------------------------------------------------
bool cmp(DynamicBitSet a, DynamicBitSet b)
{
    DynamicBitSet a_score = a&a;
    DynamicBitSet b_score = b&b;
    return (a_score >= b_score);
}

//-----------------------------------------------------------------------------------------------
// compareInterval
// return who is fitter a or b 
//-----------------------------------------------------------------------------------------------
bool compareInterval(const Pop& a, const Pop& b)
{
    return (a.fitness > b.fitness);
}

//-----------------------------------------------------------------------------------------------
// createNewGeneration
// create new generation using roulette to choose parents 
//-----------------------------------------------------------------------------------------------
void createNewGeneration(const Knapsack2DCase& knapsackCase, std::vector<Pop>& population)
{
    std::vector<Pop> newPopulation;
    int fitnessSum = std::accumulate(population.begin(), population.end(), 0, AddFitness);
    
    // set the probabilities for each pop to be chosen
    std::vector<double> probabilities;
    for (int i = 0; i < population.size(); i++)
    {
        probabilities.push_back((double)population[i].fitness / fitnessSum);
    }
    
    // dist will choose random pop based on the probabilities array 
    boost::random::discrete_distribution<> dist(probabilities);
    for (int i = 0; i < populationNumber; i += 2)
    {
        // choose parents
        int parent1Index = dist(generator);
        int parent2Index = dist(generator);
        // add the parents to new population so if thier childrens are worse the parents will survive
        newPopulation.push_back(population[parent1Index]);
        newPopulation.push_back(population[parent2Index]);
        // create childrens from the parents
        mateParents(parent1Index, parent2Index, population, newPopulation);
        // check and repair if needed the newly created chromosomes
        repairChromosome(knapsackCase, newPopulation[newPopulation.size() -2].chromosome);
        repairChromosome(knapsackCase, newPopulation[newPopulation.size() -1].chromosome);
    }
    
    // apply elitizm so the 10 best solutions can never be not selected
    // TODO: might be good idea to not add best solutions that were anyways selected...
    std::sort(population.begin(),population.end(), compareInterval);
    for (int i = 0; i < 10; i++)
        newPopulation.push_back(population[i]);
    
    population = std::move(newPopulation);
}

//-----------------------------------------------------------------------------------------------
// mateParents
// create children using one point crossover from 2 parents and add them to newPopulation
//-----------------------------------------------------------------------------------------------
void mateParents(int parent1Index, int parent2Index, std::vector< Pop >& population, std::vector<Pop>& newPopulation)
{
    DynamicBitSet& parent1 = population[parent1Index].chromosome;
    DynamicBitSet& parent2 = population[parent2Index].chromosome;
    
    int joinBitLocation = (rand() % parent1.size());
    joinBitLocation = parent1.size() - joinBitLocation;
    
    // mask starts as numberOfBits ones
    unsigned long temp = pow(2, parent1.size() + 1) - 1;
    DynamicBitSet mask = DynamicBitSet(parent1.size(), temp);
    // leave ones in the mask only from joinBitLocation to LSB
    mask = mask >> joinBitLocation;
        
    DynamicBitSet chromozomeX = parent1 & (mask.flip()); // upper  part
    DynamicBitSet chromozomeY = parent2 & mask;          // bottom part
    DynamicBitSet child = chromozomeX | chromozomeY;
    mutateChild(child);
    
    chromozomeX = parent2 & (mask.flip()); // upper  part
    chromozomeY = parent1 & mask;          // bottom part
    DynamicBitSet child2 = chromozomeX | chromozomeY;
    mutateChild(child2);
    
    newPopulation.push_back(Pop(child, 0));
    newPopulation.push_back(Pop(child2, 0));
}

//-----------------------------------------------------------------------------------------------
// mutateChild
// mutate a created child with 0.1 probabilitie for each bit to flip
//-----------------------------------------------------------------------------------------------
void mutateChild(DynamicBitSet& child)
{
    std::uniform_int_distribution<int> mutateDistributaion(0, 10);
    
    for (int i = 0; i < child.size(); i++)
    {
        if (mutateDistributaion(generator) == 0)
            child.flip(i);
    }
}

//-----------------------------------------------------------------------------------------------
// calcFitness
// calculate the chromosome fitness
// The fitness is defined by how much area the items in the chromosome cover 
// every area that is being overlapped is subtracted from the fitness
// a chromosome with overlapping items has only 10% of fitness count to discourage chromosomes
// with overlapping items 
// TODO: overlapped area is still being added as it need to be subtracted twice
//-----------------------------------------------------------------------------------------------
int calcFitness(const Knapsack2DCase& knapsackCase, const DynamicBitSet& chromosome)
{
    std::vector<Rect> itemRects;
    //std::cout << chromosome << "\n";
    //int fitness = knapsackCase.width * knapsackCase.height;
    int fitness = 0;
    bool conflict = false;
    
    const long unsigned int& coordinateBits = knapsackCase.coordinateBits;
    long unsigned int bitsPerItem = 2 + 2 * coordinateBits;
    
    unsigned int itemMaskValue = std::pow(2, bitsPerItem);
    itemMaskValue--;
    DynamicBitSet itemMask =  DynamicBitSet(chromosome.size(), itemMaskValue);
    
    for (int i = 0; i < knapsackCase.items.size(); i++)
    {
        DynamicBitSet itemValues = chromosome & itemMask;
        itemValues = itemValues >> bitsPerItem * i;
        if (itemValues[1 + 2 * coordinateBits] == 1)
        {
            unsigned int yMaskValue = std::pow(2, coordinateBits);
            yMaskValue--;
            DynamicBitSet yMask = DynamicBitSet(chromosome.size(), yMaskValue);
            unsigned int y = (itemValues & yMask).to_ulong();
            DynamicBitSet xMask = yMask << coordinateBits;
            DynamicBitSet xValue = (itemValues & xMask) >> coordinateBits;
            unsigned int x = xValue.to_ulong();
            unsigned int itemWidth = knapsackCase.items[i].width;
            unsigned int itemHeight = knapsackCase.items[i].height;
            if (itemValues[2 * coordinateBits] == 1)
            {
                std::swap(itemWidth, itemHeight);
            }
            
            itemRects.emplace_back(x,y, x + itemWidth, y + itemHeight);
        }
        
        itemMask = itemMask << bitsPerItem;
    }
    
    for (int i = 0; i < itemRects.size(); i++)
    {
        int j;
        fitness += itemRects[i].getWidth() * itemRects[i].getHeight();
        for (j = i + 1; j < itemRects.size(); j++)
        {
            Rect c = Rect::intersect(itemRects[i], itemRects[j]);
            if (c.getWidth() > 0 && c.getHeight() > 0)
            {
                conflict = true;
                fitness -= c.getWidth() * c.getHeight();
            }
        }
    }
     
    if (conflict)
        return fitness * 0.1;
    else
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

//-----------------------------------------------------------------------------------------------
// PrintSolution
// prints what items are chromosome solution and where are they positioned and whether there are
// overlapping items
//-----------------------------------------------------------------------------------------------
void PrintSolution(const Knapsack2DCase& knapsackCase,const DynamicBitSet& chromosome)
{
    std::vector<Rect> itemRects;
    bool conflict = false;
    const long unsigned int& coordinateBits = knapsackCase.coordinateBits;
    long unsigned int bitsPerItem = 2 + 2 * coordinateBits;
    
    unsigned int itemMaskValue = std::pow(2, bitsPerItem);
    itemMaskValue--;
    DynamicBitSet itemMask =  DynamicBitSet(chromosome.size(), itemMaskValue);
    
    for (int i = 0; i < knapsackCase.items.size(); i++)
    {
        DynamicBitSet itemValues = chromosome & itemMask;
        itemValues = itemValues >> bitsPerItem * i;
        if (itemValues[1 + 2 * coordinateBits] == 1)
        {
            unsigned int yMaskValue = std::pow(2, coordinateBits);
            yMaskValue--;
            DynamicBitSet yMask = DynamicBitSet(chromosome.size(), yMaskValue);
            unsigned int y = (itemValues & yMask).to_ulong();
            DynamicBitSet xMask = yMask << coordinateBits;
            DynamicBitSet xValue = (itemValues & xMask) >> coordinateBits;
            unsigned int x = xValue.to_ulong();
            unsigned int itemWidth = knapsackCase.items[i].width;
            unsigned int itemHeight = knapsackCase.items[i].height;
            if (itemValues[2 * coordinateBits] == 1)
            {
                std::swap(itemWidth, itemHeight);
            }
            
            std::cout << "Item " << i << " size " << knapsackCase.items[i].width << "X" << knapsackCase.items[i].height << " Rect " << x << " " << y << " " << x + itemWidth << " " << y + itemHeight << "\n";  
            itemRects.emplace_back(x,y, x + itemWidth, y + itemHeight);
        }
        
        itemMask = itemMask << bitsPerItem;
    }
    
    for (int i = 0; i < itemRects.size(); i++)
    {
        int j;
        for (j = i + 1; j < itemRects.size(); j++)
        {
            Rect c = Rect::intersect(itemRects[i], itemRects[j]);
            if (c.getWidth() > 0 && c.getHeight() > 0)
            {
                conflict = true;
            }
        }
    }
    
    if (conflict)
        std::cout << "There were overlapping rects\n";
}
