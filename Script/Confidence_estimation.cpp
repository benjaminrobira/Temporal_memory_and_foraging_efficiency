//////////////////////
// Function description
//////////////////////

//Estimation of the fruiting confidence based on the longterm memory

//////////////////////
// Libraries
//////////////////////

#include <math.h> // Allows the use of mathematical functions
//Note: for fmod, because doing the modulo keeping sign (i.e. respecting mathematic def), I always add lengthOfCycle, the maximum possible length, in order to avoid keeping the negative value instead of having the related positive one)
#include <cmath> // For abs and NA use
#include <ctgmath> // For NA use
#include <random> // Allows to have random distribution
#include <limits> //Allows to obtain max or min
#include <algorithm> // Allows to use sort

//////////////////////
// Function code
//////////////////////

using namespace std; // For a given function, use the function whose name is in std, the default environment of C++

//Mimics the true distribution shape: triangular with slope varying as a linear function of IS.
double confidenceEstimation(
                            double time,
                            int lengthOfCycle,
                            int lengthOfFruiting,
                            int startDateSpecies,
                            int endDateSpecies
                            ){
//NOTE: the end date is the "small" triangle right point, that is, of the distribution of starting dates. The inclusion of the +fruitinglength extension is included in the following calculation. Do not use therefore the
//operational end date + fruitingLength directly

    //Translate so 0 matches start date
    time+=-startDateSpecies;
    endDateSpecies+=-startDateSpecies;
    startDateSpecies+=-startDateSpecies;

    //Do the modulo (actually only necessary for timer since there has been a translation)
    time=fmod(time + 10*lengthOfCycle, lengthOfCycle);//10* is in case, for the value to be not negative
    startDateSpecies=fmod(startDateSpecies + 10*lengthOfCycle, lengthOfCycle);//10* is in case, for the value to be not negative

    if(endDateSpecies==lengthOfCycle || endDateSpecies == lengthOfCycle*2 || endDateSpecies==-lengthOfCycle){//Bc otherwise will be set to 0 ...
        endDateSpecies=lengthOfCycle;
    }
    else{
        endDateSpecies=fmod(endDateSpecies + 10*lengthOfCycle, lengthOfCycle);//10* is in case, for the value to be not negative
    }
    double maxDateSpecies((startDateSpecies + endDateSpecies)/2.0);

    double confidence(0.0);
    double slope(4.0/((endDateSpecies-startDateSpecies)*(endDateSpecies-startDateSpecies))); // This slope is calculated by integrating the triangular shaped curve to have a surface of 1.


    //Confidence should be equal to the area under the curve from t-30 to t which will therefore be described as:
    double valueFRight(0.0);
    double valueFLeft(0.0);

    double timeLeft(time - lengthOfFruiting);

    //Associate left value: depends on position on the triangular curve
    if(timeLeft < 0){
            if(fmod(timeLeft + 10*lengthOfCycle, lengthOfCycle) < endDateSpecies){
                valueFLeft= 1.0-
                            (-slope*(fmod(timeLeft + 10*lengthOfCycle, lengthOfCycle)-endDateSpecies)*(fmod(timeLeft + 10*lengthOfCycle, lengthOfCycle)-endDateSpecies)/2.0  + slope*(endDateSpecies + startDateSpecies)*(endDateSpecies + startDateSpecies)/8  - slope*endDateSpecies*startDateSpecies/2 + 1/2);
            }
            else{
                valueFLeft=0.0;
            }

    }
    else if(timeLeft < maxDateSpecies){
        valueFLeft= slope*(timeLeft-startDateSpecies)*(timeLeft-startDateSpecies)/2.0;
    }
    else if(timeLeft < endDateSpecies){
    //In the paper, a simpler formula is written that this coded here, but it has been showed with the following R code that you obtain the same value than this formula:

//startDateSpecies=0
//endDateSpecies=50
//currentDate=30
//slope=4/((startDateSpecies-endDateSpecies)**2)
//
//value=4*(currentDate - startDateSpecies)/(endDateSpecies-startDateSpecies) - 2*(currentDate - startDateSpecies)**2/((startDateSpecies - endDateSpecies)**2) - 1
//value
//
//value=-slope*(currentDate-endDateSpecies)*(currentDate-endDateSpecies)/2.0  +
//slope*(endDateSpecies + startDateSpecies)*(endDateSpecies + startDateSpecies)/8.0  -
//slope*endDateSpecies*startDateSpecies/2.0 + 1.0/2.0;
//value
         valueFLeft=-slope*(timeLeft-endDateSpecies)*(timeLeft-endDateSpecies)/2.0  + slope*(endDateSpecies + startDateSpecies)*(endDateSpecies + startDateSpecies)/8.0  - slope*endDateSpecies*startDateSpecies/2.0 + 1.0/2.0;
    }
    else{
        valueFLeft=1.0;
    }

    //Associate right value: depends on position on the triangular curve
    //if(time < 0){
    //can't be negative
    //}
    if(time < maxDateSpecies){
        valueFRight= slope*(time-startDateSpecies)*(time-startDateSpecies)/2.0;
    }
    else if(time < endDateSpecies){
        valueFRight=-slope*(time-endDateSpecies)*(time-endDateSpecies)/2.0  + slope*(endDateSpecies + startDateSpecies)*(endDateSpecies + startDateSpecies)/8.0  - slope*endDateSpecies*startDateSpecies/2.0 + 1.0/2.0;

    }
    else{
        valueFRight=1.0;
    }

    confidence=valueFRight-valueFLeft;
return confidence;
}


