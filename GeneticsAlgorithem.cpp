// GeneticsAlgorithem.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <algorithm>
#include <functional>
#include <stdlib.h>   
#include "time.h"
#include <list>  
#include <vector>

#include <stack>
#include <ctime>     
#include <cmath>
#include <cstdlib>  
#include <bitset>

std::vector<int> generateFirstGeneration();
void             pickOnlyTheFittests    (std::vector<int>& population);
void             createNewGeneration    (std::vector<int>& population);

int              mateParents            (int parent1, int parent2);
bool             cmp                    (int a, int b);
bool             compareInterval        (int a, int b);

//This code is trying to find the best solution
//to the following function: f(x) = x^2 where: 0 << x <<4096
//obviously, the solution is 4096 and the problem is easy, but its interesting to see how
//A genetic algorithem will behave, and that is why I am writing this code

const int populationNumber = 8;
const int numberOfBits = 12;
const int searchSpace = pow(2,numberOfBits);

int main()
{
	srand(time(0));

	//Create first generation
	std::vector<int> population = generateFirstGeneration();
	
	//Apply GA for 50 generations
	for (int gen = 0; gen < 50; gen++)
	{
		pickOnlyTheFittests(population);
		std::cout << "Gen " << gen << "'s fittest is:" << population[0] << std::endl;
		createNewGeneration(population);
		for (int i = 0; i < populationNumber; i++) std::cout << population[i] << std::endl;
	}

    return 0;
}
//--------------------------------------------------------------------------
std::vector<int> generateFirstGeneration()
{
	// the genes pool
	std::vector<int> population;
    
	for (int i = 0; i < populationNumber; i++)
		//each solution is in the range of 0 to 4095
		population.push_back( rand()%(searchSpace-1));
    
	return population;
}
//--------------------------------------------------------------------------
//input:      population
//algorithem: picks half of the population that has the best scores and "kills the rest"
//output:     the output array has the best s
void pickOnlyTheFittests(std::vector<int>& population)
{
	//sort:
	for (int i = 0; i < populationNumber; i++) std::cout << population[i] << std::endl;
		std::sort(population.begin(),population.end(), compareInterval);

	// eraes half which are the least fittest
	population.erase(population.begin() + (population.size() / 2) , population.end());
}
//--------------------------------------------------------------------------
//compare from higest to lowest scores
bool cmp(int a, int b)
{
	int a_score = a*a;
	int b_score = b*b;
	return (a_score >= b_score);
}
//--------------------------------------------------------------------------
//compare from lowest to higest scores
bool compareInterval(int a, int b)
{
	int a_score = -(a*a);
	int b_score = -(b*b);
	return (a_score > b_score);
}
//-------------------------------------------------------------------------
void createNewGeneration(std::vector<int>& population)
{
	std::random_shuffle(population.begin(), population.end());

	//fill the "killed" individuals with offsprings of fittest parents
	// create 2 children from each pair
	population.push_back(mateParents(population[0], population[1]));
	population.push_back(mateParents(population[0], population[1]));
	population.push_back(mateParents(population[2], population[3]));
	population.push_back(mateParents(population[2], population[3]));
}
//-------------------------------------------------------------------------
int mateParents(int parent1, int parent2)
{
	std::vector<int> out;
	int joinBitLocation = (rand() % numberOfBits);
	joinBitLocation = numberOfBits - joinBitLocation;
	
	// mask starts as numberOfBits ones
	int mask = pow(2, numberOfBits + 1) - 1;
	// leave ones in the mask only from joinBitLocation to LSB
    mask = mask >> joinBitLocation;
		
	int chromozomeX = parent1 & (!mask); // upper  part
	int chromozomeY = parent2 & mask;    // bottom part
	int child = chromozomeX | chromozomeY;

	if (rand() % 6 == 0)
	{
		// what bit to flip? when 0 is the LSB and 11 the MSB
		int randMutation = rand() % numberOfBits;
		int mask = 1 << randMutation;
		child ^= mask;
	}
    
	return child;
}
