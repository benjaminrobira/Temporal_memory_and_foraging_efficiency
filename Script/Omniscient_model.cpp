//////////////////////#
// Function description
//////////////////////#

//This function simulates the behaviour of a singular agent in a given environment given omniscience.


//////////////////////#
// Libraries
//////////////////////#

#include <iostream> // Allows to display input and output in the console
#include <fstream> // Allows to add/extract input/output to a file
#include <stdio.h>
#include <string.h>// To use memcpy
#include <math.h> // Allows the use of mathematical functions
//Note: for fmod, because doing the modulo keeping sign (i.e. respecting mathematic def), I always add cycleLength, the maximum possible length, in order to avoid keeping the negative value instead of having the related positive one)
#include <cmath> // For abs and NA use
#include <ctgmath> // For NA use
#include <random> // Allows to have random distribution
#include <limits> //Allows to obtain max or min
#include <algorithm> // Allows to use sort
#include <array> // Allows to use array
#include <vector> // Allows to use dynamic array and size
#include <string> // Allow to use string
#include <stdlib.h> // To convert string to double

using namespace std; // For a given function, use the function whose name is in std, the default environment of C++

/*Import own function */
#include "Food_calculation.h"
#include "Confidence_estimation.h"

//////////////////////
// Function code
//////////////////////

void omniscientModel(

        //Run ID
        int step,

        //Output variable
        std::string pathOfTheFileOutput,

        //Time
        double timer,
        double previousTimer,
        double cycleLength,

        //Agent ability
        double visionRadius,
        double coefficientDistanceToTime,
        std::string memoryType,
        bool extendibleArms,

        //Agent location and success
        double xCoordinateAgent,
        double yCoordinateAgent,

        //Initial environmental condition

        //On Species
        int numberSpecies,
        int numberIndividualsPerSpecies,
        vector<int>quantityFoodMaxTreeForSpecies,
		int numberTreesRareSpecies,
        double fruitingLength,
        int earliestStartFruitingDate,
        int latestEndFruitingDate,
        double averageMinDistanceStart,
        double sdMinDistanceStart,
        double averageMinDistanceEnd,
        double sdMinDistanceEnd,
        std::vector<double> intraSpeciesSynchrony,
        double percentSpeciesWithLowSynchrony,

        //On Tree
        std::vector<double> xCoordinateTree,
        std::vector<double> yCoordinateTree,
        std::vector<double> foodQuantityTree,
        std::vector<int> speciesTree,
        std::vector<int> endFruitingDateTree,
        std::vector<int> startFruitingDateTree,

        //Noise term=competition
        double probaTreeAlreadyDepleted

){
    //--------------
    //Opening output file
    //--------------

    std::ofstream outputFlux;
    outputFlux.open(pathOfTheFileOutput.c_str(), std::ios_base::app);//Re-open to append results
    if(outputFlux)
    {
        // Write the file only if correctly opened
    }
    else //If impossible to write in the file, say it
    {
        cout << "ERROR: Impossible to open the output file" << endl;
    }

    //--------------
    //Running model
    //--------------

    /* Random generator */
    std::random_device rd;  //Will be used to obtain a seed for the random number engine: ISSUE, IS EVERYTIME THE SAME
    std::mt19937_64 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> disU(0, 1); //Uniform distrib

    /* Variables/Parameters */
    int boutNumber(0);
    int daysElapsedWithoutFoodDuringRegistredPart(0);

    //Agent's output
    double oneStepForagingSuccessAgentTarget(0.0);
    double oneStepForagingSuccessAgentArms(0.0);
    double totFoodCollectedTarget(0.0);
    double totFoodCollectedArms(0.0);
    double totDistanceTravelled(0.0);

    //Other variables to collect for testing the trade-off exploration-exploitation hypothesis
    vector<double> foodQuantityVisitedTree;
    vector<double> targetFood;
    vector<int> totNumberVisitsTree(foodQuantityTree.size(), 0.0);
    vector<int> visitedPreviousBout(foodQuantityTree.size(), 0);

    //Initial priors on trees
    int numberTreeVisited(0);
    std::vector<double> totalFoodEatenAtTree(foodQuantityTree.size(),0);
    std::vector <int> treeVisitedDuringSimulation(foodQuantityTree.size(),0);
    std::vector <int> alreadyAteOnTreeInCycle(foodQuantityTree.size(),0);

    //Rare tree count (for rare tree experiment)
    int numberTreesRareSpeciesVisitedWithFood(0);

    int whichTarget(-1);
    int whichTargetPrevious(-1);

    while(timer<cycleLength){

        //Update vector assessing if visited during cycle for counting rare tree
        if (previousTimer <= floor(timer/ (double) cycleLength) * cycleLength && timer >= floor(timer/ (double) cycleLength) * cycleLength){
           alreadyAteOnTreeInCycle.assign(foodQuantityTree.size(),0);
        }

        vector<int> hasFoodTree(foodQuantityTree.size(), 0);
        whichTargetPrevious=whichTarget;
        whichTarget=-1;
        double interestInTarget(-1);
        //double directionTarget(0);
        double distanceToTarget(0);
        double xCoordinateTarget(-1);
        double yCoordinateTarget(-1);


        /*Target selection and move (+ food collection en route)*/
        //Associate food quantity/confidence value in it, and distance to calculate interest
        for (unsigned int t=0; t<foodQuantityTree.size(); t++){

        //Reinitialize the food quantity eaten if date is passed its ending date
        double timerTranslated_r(timer);
        double startDateTranslated_r(startFruitingDateTree[t]);
        double endDateTranslated_r(endFruitingDateTree[t]);

        //Translate so 0 matches start date
        timerTranslated_r+=-startDateTranslated_r;
        endDateTranslated_r+=-startDateTranslated_r;
        //startDateTranslated+=-startDateTranslated;

        //Do the modulo (actually only necessary for timer since there has been a translation)
        timerTranslated_r=fmod(timerTranslated_r + 10*cycleLength, cycleLength);
        //startDateTranslated=fmod(startDateTranslated + 10*cycleLength, cycleLength);
        endDateTranslated_r=fmod(endDateTranslated_r + 10*cycleLength, cycleLength);

        if(
           timerTranslated_r > endDateTranslated_r
           ){
            totalFoodEatenAtTree[t]=0;
        }

        //Update the food quantity at the tree
        foodQuantityTree[t]=foodQuantity(timer, startFruitingDateTree[t], endFruitingDateTree[t],quantityFoodMaxTreeForSpecies[speciesTree[t]], totalFoodEatenAtTree[t], cycleLength);

        if(foodQuantityTree[t]>0.0001){
            hasFoodTree[t]=1;
        }
        if((int) floor(previousTimer)!=(int) floor(timer)){
            if(foodQuantityTree[t]>0.0001){
                if(disU(gen) < probaTreeAlreadyDepleted){//If competition, food, if present, can be depleted
                    foodQuantityTree[t]=0.0;
                    totalFoodEatenAtTree[t]=1.0;//To mimic depletion, bc this variable is only used for estimating food quantity, consider that the agent already depleted it
                    hasFoodTree[t]=0;
                }
           }
        }

        //Then proceed to association confidence/quantity estimation
        double distanceToTree(
                              sqrt(
                              (xCoordinateAgent - xCoordinateTree[t])*(xCoordinateAgent - xCoordinateTree[t]) +
                              (yCoordinateAgent - yCoordinateTree[t])*(yCoordinateAgent - yCoordinateTree[t])
                              )
                              );

        double interestInCurrentTree(foodQuantityTree[t]/distanceToTree);

        //Added to avoid targeting trees in the starting vision field
        if(distanceToTree <= visionRadius){
            //Not done since now, it is eaten at the start with the telescopic arms, so not targeted anymore
            //foodQuantityEstimated=foodQuantityTree[t];
            //So put negative value not to be targeted anymore
            interestInCurrentTree=-1;

        }

        //Avoiding retargeting the target; actually useless with omniscient but was done when resource could "regrow" on earlier versions
        if(whichTargetPrevious==t){
            interestInCurrentTree=-1;
        }
        else if (distanceToTree > 0){
            if(interestInCurrentTree > interestInTarget){//Change target to that tree if of higher interest
                whichTarget=t;
                interestInTarget=interestInCurrentTree;
                //directionTarget=atan2((yCoordinateTree[t] - yCoordinateAgent),(xCoordinateTree[t] - xCoordinateAgent));
                distanceToTarget=distanceToTree;
                xCoordinateTarget=xCoordinateTree[t];
                yCoordinateTarget=yCoordinateTree[t];
            }
            else if(interestInCurrentTree == interestInTarget){
                    //cout << "Interest is equal" << endl;
                    if (distanceToTarget > distanceToTree){
                        //Take the one with the higher estimation of food (=> greater distance)
                        whichTarget=t;
                        interestInTarget=interestInCurrentTree;
                        //directionTarget=atan2((yCoordinateTree[t] - yCoordinateAgent),(xCoordinateTree[t] - xCoordinateAgent));
                        distanceToTarget=distanceToTree;
                        xCoordinateTarget=xCoordinateTree[t];
                        yCoordinateTarget=yCoordinateTree[t];
                    }
                    else{
                        //Keep the first one
                    }
            }
            else{
                //Do nothing
                }
            }
        }

  double foragingTimeTarget(0);
  double foodTarget(foodQuantityTree[whichTarget]);
  double foodEatenTarget(totalFoodEatenAtTree[whichTarget]);
//********************************************************************************************************************************************************************************************
//**MEMORY
        if((int) accumulate(hasFoodTree.begin(), hasFoodTree.end(), (int) 0)== (int) 0){//Do not run simulation if no tree with food
            previousTimer=timer;
            timer=floor(timer +1);
            if(timer >= earliestStartFruitingDate){
                daysElapsedWithoutFoodDuringRegistredPart+=1;
            }
            whichTarget=whichTargetPrevious;//Reput the correct previous target since did not move
        }
//********************************************************************************************************************************************************************************************
        else{
            //Register food at target
            targetFood.push_back(foodQuantityTree[whichTarget]);

            double foodCollectedTarget(0);
            double foodCollectedArms(0);


            //Determine all trees on the way to the target and then update their food quantity and foraging success accordingly
            for (unsigned int t=0; t<foodQuantityTree.size(); t++){
                double xCoordinateTreeSystemAgentOriginTarget((xCoordinateTree[t]-xCoordinateAgent)*(xCoordinateTarget-xCoordinateAgent)/distanceToTarget+ (yCoordinateTree[t]-yCoordinateAgent)*(yCoordinateTarget-yCoordinateAgent)/distanceToTarget);
                double yCoordinateTreeSystemAgentOriginTarget(-(xCoordinateTree[t]-xCoordinateAgent)*(yCoordinateTarget-yCoordinateAgent)/distanceToTarget+ (yCoordinateTree[t]-yCoordinateAgent)*(xCoordinateTarget-xCoordinateAgent)/distanceToTarget);
                if(t==(unsigned) whichTarget){//To be sure select target
                                totalFoodEatenAtTree[t]+=foodQuantityTree[t];
                                foodCollectedTarget+= foodQuantityTree[t];
                                //foodQuantityTree[t]=0; //Useless bc accounting for eaten quantity before

                                if(timer >=0 &&
                                //Do not count trees that were just visited: therefore in this case, those included in the initial vision circle. For NM/NWM/chrono/asso, using the "seenDuringBout" variable. Here this var doesn't exist.
                                visitedPreviousBout[t]==0){
                                    //Update if tree visited (==SEEN)
                                    treeVisitedDuringSimulation[t]=1;
                                    numberTreeVisited+=1;
                                    totNumberVisitsTree[t]+=1;
                                    //Add quantity food during visits
                                    foodQuantityVisitedTree.push_back(foodQuantityTree[t]);
                                }

                                //Rare tree count
                                if(speciesTree[t]==numberSpecies-1 && foodQuantityTree[t] > 0.0001 && timer > 0 && timer < cycleLength && alreadyAteOnTreeInCycle[t]==0){
                                    numberTreesRareSpeciesVisitedWithFood+=1;
                                }
                                if(foodQuantityTree[t] > 0.0001){
                                    alreadyAteOnTreeInCycle[t]=1;
                                }
                                visitedPreviousBout[t]=1;
                }
                else if(abs(yCoordinateTreeSystemAgentOriginTarget)<=visionRadius&&//The tree encountered en route is no more than visionRadius away from the linear path to the target
                   xCoordinateTreeSystemAgentOriginTarget >= 0 && //This tree is in the same direction than the target
                   xCoordinateTreeSystemAgentOriginTarget <= distanceToTarget // This tree is not farther than the target
                   ){
                        if(extendibleArms==true){
                            //Update food quantity collected
                            totalFoodEatenAtTree[t]+=foodQuantityTree[t];
                            foodCollectedArms+= foodQuantityTree[t];
                            //Rare tree count
                            if(speciesTree[t]==numberSpecies-1 && foodQuantityTree[t] > 0.0001 && timer > 0 && timer < cycleLength && alreadyAteOnTreeInCycle[t]==0){
                                    numberTreesRareSpeciesVisitedWithFood+=1;
                            }
                            if(foodQuantityTree[t] > 0.0001){
                                alreadyAteOnTreeInCycle[t]=1;
                            }
                        }
                        else{
                            //Do nothing
                        }
                        if(timer >=0 &&
                           //Do not count trees that were just visited: therefore in this case, those included in the initial vision circle. For NM/NWM/chrono/asso, using the "seenDuringBout" variable. Here this var doesn't exist.
                           visitedPreviousBout[t]==0
                           ){
                            //Update if tree visited (==SEEN)
                            totNumberVisitsTree[t]+=1;
                            treeVisitedDuringSimulation[t]=1;
                            numberTreeVisited+=1;
                            //Add quantity food during visits
                            foodQuantityVisitedTree.push_back(foodQuantityTree[t]);
                        }
                        visitedPreviousBout[t]=1;
                }
                else if(//Do the same for all food centered on the target location within the vision circle
                        (xCoordinateTree[t]-xCoordinateTarget)*(xCoordinateTree[t]-xCoordinateTarget) +
                        (yCoordinateTree[t]-yCoordinateTarget)*(yCoordinateTree[t]-yCoordinateTarget)<=visionRadius*visionRadius ||
                        //Do the same for all food centered on the start location
                        (xCoordinateTree[t]-xCoordinateAgent)*(xCoordinateTree[t]-xCoordinateAgent) +
                        (yCoordinateTree[t]-yCoordinateAgent)*(yCoordinateTree[t]-yCoordinateAgent)<=visionRadius*visionRadius
                        ){

                        if(extendibleArms==true){
                            //Update food quantity collected
                            totalFoodEatenAtTree[t]+=foodQuantityTree[t];
                            foodCollectedArms+= foodQuantityTree[t];
                            //Rare tree count
                            if(speciesTree[t]==numberSpecies-1 && foodQuantityTree[t] > 0.0001 && timer > 0 && timer < cycleLength && alreadyAteOnTreeInCycle[t]==0){
                                    numberTreesRareSpeciesVisitedWithFood+=1;
                            }
                            if(foodQuantityTree[t] > 0.0001){
                                alreadyAteOnTreeInCycle[t]=1;
                            }
                        }
                        else{
                        }
                            //foodQuantityTree[t]=0; //Useless bc accounting for eaten quantity before

                        if(timer>=0 &&
                           //Do not count trees that were just visited: therefore in this case, those included in the initial vision circle. For NM/NWM/chrono/asso, using the "seenDuringBout" variable. Here this var doesn't exist.
                           visitedPreviousBout[t]==0){
                            //Update if tree visited (==SEEN)
                            totNumberVisitsTree[t]+=1;
                            treeVisitedDuringSimulation[t]=1;
                            numberTreeVisited+=1;
                            //Add quantity food during visits
                            foodQuantityVisitedTree.push_back(foodQuantityTree[t]);
                        }
                        visitedPreviousBout[t]=1;
                }
                else{
                    visitedPreviousBout[t]=0;
                }
            }
            /*Update agents' features*/
            if(timer >= 0){
                totDistanceTravelled+=distanceToTarget;
                totFoodCollectedArms+=foodCollectedArms;
                totFoodCollectedTarget+=foodCollectedTarget;
                oneStepForagingSuccessAgentArms+=(foodCollectedArms + foodCollectedTarget)/distanceToTarget;
                oneStepForagingSuccessAgentTarget+=foodCollectedTarget/distanceToTarget;
                boutNumber+=1;
            }
            xCoordinateAgent=xCoordinateTree[whichTarget];
            yCoordinateAgent=yCoordinateTree[whichTarget];
            previousTimer=timer;
            timer+=coefficientDistanceToTime*distanceToTarget;
        }
    }

    cout << "Ending omniscient" << endl;
    //--------------
    //Ending procedure
    //--------------

    /*Save results in a final array*/
    outputFlux <<
    step << " " <<
    memoryType << " " <<
    timer<< " " <<
    daysElapsedWithoutFoodDuringRegistredPart << " " <<
    oneStepForagingSuccessAgentArms/boutNumber << " " <<
    oneStepForagingSuccessAgentTarget/boutNumber << " " <<
    (totFoodCollectedArms+totFoodCollectedTarget)/totDistanceTravelled << " " <<
    totFoodCollectedTarget/totDistanceTravelled << " " <<
    std::accumulate(treeVisitedDuringSimulation.begin(), treeVisitedDuringSimulation.end(), 0.0)/(double) treeVisitedDuringSimulation.size() << " " <<
    std::accumulate(targetFood.begin(), targetFood.end(), 0.0)/(double) targetFood.size() << " " <<
    numberTreesRareSpecies << " " <<
    numberTreesRareSpeciesVisitedWithFood/ (double) numberTreesRareSpecies << " " <<
    probaTreeAlreadyDepleted << " " <<
    "NA" << " " <<
    "NA" << " " <<
    visionRadius << " " <<
    cycleLength << " " <<
    fruitingLength << " " <<
    numberSpecies << " " <<
    numberIndividualsPerSpecies << " " <<
    std::accumulate(intraSpeciesSynchrony.begin(), intraSpeciesSynchrony.end(), 0.0)/(double) intraSpeciesSynchrony.size() << " " <<
    percentSpeciesWithLowSynchrony << " " <<
    averageMinDistanceStart << " " <<
    sdMinDistanceStart << " " <<
    averageMinDistanceEnd << " " <<
    sdMinDistanceEnd << " " <<
    "NA" <<
    endl;

    //Close output
    outputFlux.close();
}






