//////////////////////
// Function description
//////////////////////

//Calculates the fruit quantity on a tree given the start and end date of fruiting, and how much food has been already eaten. Note that this food evolution is a parabole.

//////////////////////
// Libraries
//////////////////////

#include "Food_calculation.h"
#include <math.h> // Allows the use of mathematical functions
//Note: for fmod, because doing the modulo keeping sign (i.e. respecting mathematic def), I always add cycleLength, the maximum possible length, in order to avoid keeping the negative value instead of having the related positive one)
#include <cmath> // For abs and NA use
#include <ctgmath> // For NA use
#include <random> // Allows to have random distribution
#include <limits> //Allows to obtain max or min
#include <algorithm> // Allows to use sort

using namespace std; // For a given function, use the function whose name is in std, the default environment of C++

//////////////////////
// Function code
//////////////////////

double foodQuantity (double date, int startDate, int endDate, int maxFruitQuantity, double quantityEaten, double lengthOfCycle)
{
    //Translate so 0 matches start date
    date+=-startDate;
    endDate+=-startDate;
    startDate+=-startDate;

    //Do the modulo (actually only necessary for timer since there has been a translation)
    date=fmod(date + 10*lengthOfCycle, lengthOfCycle);//10* is in case, for the value to be not negative
    startDate=fmod(startDate + 10*lengthOfCycle, lengthOfCycle);//10* is in case, for the value to be not negative
    endDate=fmod(endDate + 10*lengthOfCycle, lengthOfCycle);//10* is in case, for the value to be not negative

    double quantityAvailable(0.0);
    if(date >= startDate && date <= endDate){
       quantityAvailable = maxFruitQuantity - quantityEaten;
    }
    else{
      quantityAvailable=0.0;
    }

    return quantityAvailable;
}



