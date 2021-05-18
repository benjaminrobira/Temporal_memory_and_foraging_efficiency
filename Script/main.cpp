//------------------------------------------
// LOADING NECESSARY LIBRARIES AND FUNCTIONS
//------------------------------------------


/*Import own libraries*/
#include <iostream> // Allows to display input and output in the console
#include <fstream> // Allows to add/extract input/output to a file
#include <stdio.h>
#include <string.h>// To use memcpy
#include <sstream> // To transform double into string
#include <math.h> // Allows the use of mathematical functions
//Note: for fmod, because doing the modulo keeping sign (i.e. respecting mathematic def), I always add cycleLengthInit, the maximum possible length,
//in order to avoid keeping the negative value instead of having the related positive one)
#include <cmath> // For abs and NA use
#include <ctgmath> // For NA use
#include <random> // Allows to have random distribution
#include <limits> //Allows to obtain max or min
#include <algorithm> // Allows to use sort, fill
#include <array> // Allows to use array
#include <vector> // Allows to use dynamic array and size
#include <string> // Allow to use string
#include <stdlib.h> // To convert string to double
#include<omp.h> //For using OpenMP to parallelize independent for loop
#include <thread>
//I have -fopenmp in my compiler "Compiler->Compiler Settings->Other Options"
//I have -lgomp  in "Compiler->Linker Settings->Other linker options" (-lpthread and -lmthread might be necessary too)
//-mthread  in "Compiler->Linker Settings->Other linker options" if windows and exceptions for parallelization (i.e. private ...)
//I have added the libgomp-1.dll file and other with the location of the exe file in case


/*Import own function*/
#include "Food_calculation.h"
//#include "Interest_estimation.h"
#include "Confidence_estimation.h"
#include "Null_model.h"
#include "Null_working_memory_model.h"
#include "Chronological_model.h"
#include "Chronological_model_for_delay.h"
#include "Associative_model.h"
#include "Omniscient_model.h"

//NOTE:
//For associative and chronological memories, an earlier version included a test if they were able to compute a sampling strategy and combine it with long-term knowledge
// To do so a weight, varying between -1 (no long-term knwoledge used) to 1(only longterm knowledge used) was implemented. The code is left with this but only a weight of 1 is used in the article simulations.

using namespace std; // For a given function, use the function whose name is in std, the default environment of C++

//-------------------
// Main function: divided for each scenario to test namely Homogeneous scenario, Heterogeneous scenario, Working memory span and environment composition, Delay Experiment
// Each scenario will start with setting environment. Then each memory will be tested!
//-------------------

int main()
{
    /*Opening output and preparing column header*/

    //Output path
    string pathOfTheFileOutputInit("output_test");

//    outputFlux.open(pathOfTheFileOutputInit.c_str(), std::ios_base::app);//Re-open to append results
//    if(outputFlux)  // Write the file only if correctly opened
//    {
//    }
//    else //If impossible to write in the file, say it
//    {
//        cout << "ERROR: Impossible to open the output file" << endl;
//    }


    /*Seeding random sequence */
    std::random_device rd;  //Will be used to obtain a seed for the random number engine: ISSUE, IS EVERYTIME THE SAME
    std::mt19937_64 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    //--------------
    //FIXED VARIABLES
    //--------------


    int numberSpeciesInit(30);
    int numberIndividualsPerSpeciesInit(100/*Import own function */);
    vector<int> quantityFoodMaxTreeForSpeciesInit(numberSpeciesInit, 1);
    int numberTreesRareSpeciesInit(0);
	int mapSizeInit(1000);
    double visionRadiusInit((double) 1/(double) 2/sqrt((double) numberSpeciesInit* (double) numberIndividualsPerSpeciesInit/ (double) mapSizeInit/(double) mapSizeInit));
    //The lowest limit for the radius must be 1/(2sqrt(d)) with d the density. This value represents the average closest neighbour distance expected.
    //However, it must be chosen to also be reduced compared to the distance travelled.
    int cycleLengthInit(360);
    int fruitingLengthInit(30);
    double timeLapseNotToReturnInit(15.0);

    int numberPatch(10);
    double dispersion(20);
///////////////////////////////////////////////////////////////////////
// WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ///////////////////////////////
    // IMPORTANT NOTE /////
    //////////////////////////////

    //If cycleLength=365
    //The model cannot work properly for the associative/confidence part when the triangular-shaped distribution has a base of more than base=335 (i.e. cycleLength - fruitingLength), because of the confidence function, which takes
    //a new triangular-shaped curve for confidence estimation, of base'=base + fruitinglength, because 1) the ending date is at a distance greater than cycleLength, which is taken as a modulo for the date and timer and
    //2) because if another modulo is taken, actually trees ending their fruiting would in the end occur at the same time than some just starting, hence we should consider a different confidence as the sum of the end
    //of the triangle and the beginning of the other. To avoid such useless calculations, please constrain values of IS (intraspecies synchrony) equal or above 1- (cycleLength - fruitingLength)/cycleLength=1-335/cycleLengthInit=0.08

///////////////////////////////////////////////////////////////////////

    std::vector<double> intraSpeciesSynchronyValuesInit={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90};

    //--------------
    //FIXED FUNCTIONS
    //--------------

    std::uniform_real_distribution<double> disULocTree(0, mapSizeInit); //Uniform distrib
    std::uniform_real_distribution<double> disULocPatch(sqrt(6)*mapSizeInit/dispersion, mapSizeInit - sqrt(6)*mapSizeInit/dispersion); //Uniform distrib: the sqrt(6) define where 95% will fall for a bivariate gaussian.
    //Hence, as I use mapSize/dispersion as the sd for the gaussian, it will say that for the patches, tree fall at 95% in a radius of sqrt(6)*mapSizeInit/dispersion from the center with coordinates the means of the bivariate gaussian,
    // I therefore substract this radius from the borders to avoid border effects
    std::uniform_real_distribution<double> disUDateSpecies(0,cycleLengthInit);//Uniform distrib, getting rid of border effects

    string isHeterogeneous("NO");

    ////////////////////////////////////////////////////////////////////
    // SCENARIO 1: Temporally/Spatially homogeneous environment with changing synchrony
    ////////////////////////////////////////////////////////////////////

    #pragma omp parallel for
    for (int s=0; s<300; s++){

        std::string pathOfTheFileOutputInitMdf(pathOfTheFileOutputInit);
        pathOfTheFileOutputInitMdf.append(to_string(s));
        pathOfTheFileOutputInitMdf.append("_IS_test.txt");
        std::ofstream outputFlux;
        outputFlux.open(pathOfTheFileOutputInitMdf.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
        if(outputFlux)  // Write the file only if correctly opened
        {
        }
        else //If impossible to write in the file, say it
        {
            cout << "ERROR: Impossible to open the output file" << endl;
        }

        outputFlux << "Run_simulation" << " " << "Memory_type" << " " << "Time_length_run" << " " << "Days_wo_food" << " " << "Foraging_success_cum" << " " << "Foraging_success_cum_target_only" << " " << "Foraging_success_tot" << " " << "Foraging_success_tot_target_only" << " " << "Percent_tree_visited" << " " <<
        "Percent_targets_with_food" << " " << "Number_trees_rare_species" << " " << "Rate_trees_rare_eaten" << " " << "Proba_competition" << " " << "Working_memory_length_time_unit" << " " << "Time_not_to_return" << " " << "Vision_radius" << " " << "Length_cycle" << " " << "Length_fruiting" << " " << "Number_species" << " " << "Number_individuals_per_species" << " " <<
        "Intra_species_synchrony" << " " << "Percentage_low_synchrony" << " " << "Distance_min_start" << " " << "Sd_min_start" << " " << "Distance_min_end" << " " << "Sd_min_end" << " " << "Weight_for_longterm_memory" << endl;

        outputFlux.close();
        for(int iSS=0; iSS<intraSpeciesSynchronyValuesInit.size(); iSS++){

            //--------------
            //INITIALISATION
            //--------------

            /*Time-associated variables*/
            double timerInit(-cycleLengthInit);
            double previousTimerInit(-cycleLengthInit);
            double coefficientDistanceToTimeInit((double) 1/((double) 1/cycleLengthInit*(double) mapSizeInit*(double) mapSizeInit/(double) visionRadiusInit/2.0));
            //cout << 1/coefficientDistanceToTimeInit << endl;

            // This is set up by making the hypothesis that during the cycle all of the trees can be visited. We assume that, during the cycle the individual hence travels a distance d
            // which makes it cover d*visionRadius surface (approximation not considering overlap, hence overestimation). Applying the "produit en croix", we find that d can't be more than 1/2*mapSize**2/visionRadius. Hence this speed if defining working memory as FruitingLength/2

            /*Species-level variables*/

            vector<int> startFruitingDateSpeciesOperationalInit;
            vector<int> endFruitingDateSpeciesOperationalInit;

            vector<int> startFruitingDateSpeciesAverageInitResult;
            vector<int> startFruitingDateSpeciesStdevInitResult;

            vector<double> intraSpeciesSynchronyInit;

            /*Tree-level variables*/
            vector<double> xCoordinateTreeInit;
            vector<double> yCoordinateTreeInit;
            vector<int> startFruitingDateTreeInit;
            vector<int> endFruitingDateTreeInit;
            vector<int> speciesTreeInit;
            vector<double> foodQuantityTreeInit;

            /*Agent-level variables*/
            double xCoordinateAgentInit(disULocTree(gen));
            double yCoordinateAgentInit(disULocTree(gen));
            double workingMemoryTimeInit(15.0);
            double weightInFavourOfLongtermKnowledgeInit(1.0); //Only for associative/chronological memory, should be comprised between [-1;1]. Used for averaging between observed (=estimated) fruiting, and Longterm knowledge on fruiting dates

            /* Environment */
            double earliestStartFruitingDateInit(cycleLengthInit);
            double latestEndFruitingDateInit(0.0);

            std::uniform_real_distribution<double> disU01(0, 1); //Uniform distrib
            std::bernoulli_distribution bernouilliDisInit(0.5);

            /*Setting the environment: resource distribution and tree features*/
            for (int i=0; i<numberSpeciesInit; i++){

                    intraSpeciesSynchronyInit.push_back(intraSpeciesSynchronyValuesInit[iSS]);
                    //int speciesToSplit(0);
                    vector<int> dateStartVector;
                    vector<int> dateEndVector;
                        int speciesStartDate(disUDateSpecies(gen));

                        //ADDED ON THE 01/02/2021
                        //If heterogeneous distribution
                        std::vector<int> patchesCoordinatesX(numberPatch,0);
                        std::vector<int> patchesCoordinatesY(numberPatch,0);
                        for(int p=0; p<patchesCoordinatesX.size(); p++){
                            patchesCoordinatesX[p]=disULocPatch(gen);
                            patchesCoordinatesY[p]=disULocPatch(gen);
                        }
                        int counterTree=0;
                        int refPatch=0;
                        std::normal_distribution<double> disNormalPatchX(patchesCoordinatesX[0],mapSizeInit/dispersion);//
                        std::normal_distribution<double> disNormalPatchY(patchesCoordinatesY[0],mapSizeInit/dispersion);//
                        //////////////////////////

                        for (int j=0; j<numberIndividualsPerSpeciesInit; j++){
                            //New method triangular: sum of uniform distributions (NOTE: in the article for formula, we are using the Tm (the max point), and not the TS, therefore you can see that it sums the uniform - 1 !!!!!!!!
                            int dateTransitory(speciesStartDate +  (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit/2*(disU01(gen)+disU01(gen)));

                            startFruitingDateTreeInit.push_back(dateTransitory);
                            endFruitingDateTreeInit.push_back(dateTransitory + fruitingLengthInit);

                            //Food quantity
                            foodQuantityTreeInit.push_back(0);//No food at first (will be updated afterwards)

                            //Coordinates

                            /////////// ADDED ON THE 01/02/2021
                            if(isHeterogeneous=="YES"){
                                counterTree+=1;
                                if(counterTree>=floor(numberIndividualsPerSpeciesInit/(double) numberPatch)){
                                    if(refPatch!=patchesCoordinatesX.size()-1){
                                        refPatch+=1;
                                        counterTree=0;
                                        disNormalPatchX.param(std::normal_distribution<double>::param_type(patchesCoordinatesX[refPatch],mapSizeInit/dispersion));//
                                        disNormalPatchY.param(std::normal_distribution<double>::param_type(patchesCoordinatesY[refPatch],mapSizeInit/dispersion));//
                                    }
                                }
                                double locX(disNormalPatchX(gen));
                                while(locX < 0 || locX > mapSizeInit){
                                    locX=disNormalPatchX(gen);
                                }
                                double locY(disNormalPatchY(gen));
                                while(locY < 0 || locY > mapSizeInit){
                                    locY=disNormalPatchY(gen);
                                }
                                xCoordinateTreeInit.push_back(locX);
                                yCoordinateTreeInit.push_back(locY);
                            }
                            else{
                            ////////////////////////////////////
                                xCoordinateTreeInit.push_back(disULocTree(gen));
                                yCoordinateTreeInit.push_back(disULocTree(gen));
                            }

                            //Update vector for mean and var calculation
                            dateStartVector.push_back(dateTransitory);
                            dateEndVector.push_back(dateTransitory+fruitingLengthInit);

                            speciesTreeInit.push_back(i);
                        }

                    /*Calculating the operational starting date*/
                    //std::sort (dateStartVector.begin(), dateStartVector.end());//Ranging in increasing order starting dates
                    startFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate + 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));//Correction adding bc subsampling triggers that we never reach this extreme value. Hence Simon conducted simulation to see what correction we have to apply when sampling in triangular distrib, the extreme values we obtained on average.

                    /*Calculating the operational ending date*/
                    //std::sort (dateEndVector.rbegin(), dateEndVector.rend());//Ranging in decreasing order ending dates
                    endFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate  + (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit - 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));

                    if(fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit) < earliestStartFruitingDateInit
                       ){
                            earliestStartFruitingDateInit=fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                    }
                    if(fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit) > latestEndFruitingDateInit
                       ){
                            latestEndFruitingDateInit=fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                    }

                    /*Extract the finally obtained mean and var for all species*/
                    //Calculate sum, mean, square sum and stdv
                    double sum = std::accumulate(dateStartVector.begin(), dateStartVector.end(), 0.0);//arguments are: first value, last value, initialisation value
                    double mean = sum / dateStartVector.size();
                    double sq_sum = std::inner_product(dateStartVector.begin(), dateStartVector.end(), dateStartVector.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
                    double stdev = std::sqrt(sq_sum / dateStartVector.size() - mean * mean);

                    if(intraSpeciesSynchronyInit[i] > 0){
                    }
                    else{
                        mean=cycleLengthInit/2.0;
                        stdev=sqrt(cycleLengthInit*cycleLengthInit/12.0);
                    }
                    startFruitingDateSpeciesAverageInitResult.push_back(mean);
                    startFruitingDateSpeciesStdevInitResult.push_back(stdev);
                }

                earliestStartFruitingDateInit=max(0.0, earliestStartFruitingDateInit);
                latestEndFruitingDateInit=min((double) cycleLengthInit, latestEndFruitingDateInit);

                /*Determinate the associative strength between start-start events and start-end events*/
                //and /*Create the matrix for associative model: first col is the species, second is the distance to the date, and third col characterises which date is to use 0=start, 1=end*/

                //Create the result matrix
                int distanceStartStartInit[numberSpeciesInit][numberSpeciesInit];
                int distanceStartEndInit[numberSpeciesInit][numberSpeciesInit];
                int distanceEndStartInit[numberSpeciesInit][numberSpeciesInit];
                int distanceEndEndInit[numberSpeciesInit][numberSpeciesInit];
                int matrixInferenceAssociationInit[numberSpeciesInit][3];
                std::fill(*matrixInferenceAssociationInit, *matrixInferenceAssociationInit + 3*numberSpeciesInit, cycleLengthInit); //Initialise matrix with only cycleLengthInit

                //Create the min vectors of distance
                vector<double> vectorMinDistanceStart;
                vector<double> vectorMinDistanceEnd;

                //Loop to fill the result matrices and calculate alongside the minimum: first rows are always start-start then start-end
                for (int j=0; j<numberSpeciesInit; j++){

                    for (int k=0; k<numberSpeciesInit; k++){
                    if(j!=k){
                        //Updating the results matrices
                            distanceStartStartInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));
                            distanceStartEndInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]-fruitingLengthInit),cycleLengthInit));
                            distanceEndStartInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]+fruitingLengthInit),cycleLengthInit));
                            distanceEndEndInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));

                            vectorMinDistanceStart.push_back(min(distanceStartStartInit[j][k],
                                                                 distanceStartEndInit[j][k]
                                                                 )
                                                             );
                            vectorMinDistanceEnd.push_back(min(distanceEndStartInit[j][k],
                                                               distanceEndEndInit[j][k]
                                                               )
                                                           );
                        if(distanceStartStartInit[j][k]==min(distanceStartStartInit[j][k],distanceStartEndInit[j][k])){
                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                               (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] >=
                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                            ){
                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                       (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] ==
                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]

                                    ){
                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                            ){
                                               matrixInferenceAssociationInit[j][0]=k;
                                               matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                               matrixInferenceAssociationInit[j][2]=0;
                                            }
                                            else{
                                                //Do nothing
                                            }
                                    }
                                    else{
                                         matrixInferenceAssociationInit[j][0]=k;
                                         matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                         matrixInferenceAssociationInit[j][2]=0;
                                    }
                            }
                            else{
                                //Do nothing
                            }
                        }
                        else{
                           if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                               (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] >=
                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                            ){
                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                       (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] ==
                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]

                                    ){
                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                            ){
                                               matrixInferenceAssociationInit[j][0]=k;
                                               matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                               matrixInferenceAssociationInit[j][2]=1;
                                            }
                                            else{
                                                //Do nothing
                                            }
                                    }
                                    else{
                                         matrixInferenceAssociationInit[j][0]=k;
                                         matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                         matrixInferenceAssociationInit[j][2]=1;
                                    }
                            }
                        }
                    }
                }
            }

            //Calculate sum, mean, square sum and stdv
            double sum = std::accumulate(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), 0.0);//arguments are: first value, last value, initialisation value
            double mean = sum / vectorMinDistanceStart.size();
            double sq_sum = std::inner_product(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), vectorMinDistanceStart.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
            double stdev = std::sqrt(sq_sum / vectorMinDistanceStart.size() - mean * mean);

            double averageMinDistanceStartInit(mean);
            double sdMinDistanceStartInit(stdev);

            sum = std::accumulate(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), 0.0);//arguments are: first value, last value, initialisation value
            mean = sum / vectorMinDistanceEnd.size();
            sq_sum = std::inner_product(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), vectorMinDistanceEnd.begin(), 0.0);//arguments are: first value to End with, last value, first value to match to the other one, initialisation value
            stdev = std::sqrt(sq_sum / vectorMinDistanceEnd.size() - mean * mean);

            double averageMinDistanceEndInit(mean);
            double sdMinDistanceEndInit(stdev);


        //--------------
        //Running models
        //--------------

       /* Run model: NULL*/

        cout << "I work until line setting loop for model: NULL" << endl;

        nullModel(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double workingMemoryTime=
                workingMemoryTimeInit, //max(2.0, (double)fruitingLengthInit/20),
                //double timeLapseNotToReturn=
                timeLapseNotToReturnInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Null_model_telescopic",
                // bool extendibleArms,
                true,

                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0

        );

        /*Run model: NULL WORKING MEMORY*/
        cout << "I work until line setting loop for model: NULL WORKING MEMORY" << endl;

        nullModelWorkingMemory(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double workingMemoryTime=
                workingMemoryTimeInit,
                //double timeLapseNotToReturn=
                timeLapseNotToReturnInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Null_working_memory_model_telescopic",
                // bool extendibleArms,
                true,

                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0
        );

        /*Run model: CHRONOLOGICAL*/
        cout << "I work until line setting loop for model: CHRONOLOGICAL" << endl;

        chronologicalMemoryModel(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double workingMemoryTime=
                workingMemoryTimeInit,
                //double timeLapseNotToReturn=
                timeLapseNotToReturnInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Chronological_model_telescopic",
                // bool extendibleArms,
                true,
                //double weightInFavourOfLongtermKnowledge,
                weightInFavourOfLongtermKnowledgeInit,

                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,
                //std::vector<int> startFruitingDateSpeciesOperational=
                startFruitingDateSpeciesOperationalInit,
                //std::vector<int> endFruitingDateSpeciesOperational=
                endFruitingDateSpeciesOperationalInit,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0
        );


        /*Run model: ASSOCIATIVE*/
        cout << "I work until line setting loop for model: ASSOCIATIVE" << endl;

        associativeMemoryModel(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //-cycleLengthInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double workingMemoryTime=
                workingMemoryTimeInit,
                //double timeLapseNotToReturn=
                timeLapseNotToReturnInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Associative_model_telescopic",
                // bool extendibleArms,
                true,
                //double weightInFavourOfLongtermKnowledge,
                weightInFavourOfLongtermKnowledgeInit,


                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,
                //int matrixInferenceAssociation,
                *matrixInferenceAssociationInit,
                //std::vector<int> startFruitingDateSpeciesOperational=
                startFruitingDateSpeciesOperationalInit,
                //std::vector<int> endFruitingDateSpeciesOperational=
                endFruitingDateSpeciesOperationalInit,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0
        );

        /*Run model: OMNISCIENT*/
        cout << "I work until line setting loop for model: OMNISCIENT" << endl;

        omniscientModel(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Omniscient_model_telescopic",
                // bool extendibleArms,
                true,

                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0
        );
    }
}

////////////////////////////////////////////////////////////////
// SCENARIO 2: Temporally heterogeneous/spatially homogeneous environment with change in proportion of little/highly synchronous species
////////////////////////////////////////////////////////////////

    #pragma omp parallel for
    for (int s=0; s<300; s++){
        cout << s << endl;
        std::string pathOfTheFileOutputInitMdf(pathOfTheFileOutputInit);
        pathOfTheFileOutputInitMdf.append(to_string(s));
        pathOfTheFileOutputInitMdf.append("_HE.txt");
        std::ofstream outputFlux;
        outputFlux.open(pathOfTheFileOutputInitMdf.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
        if(outputFlux)  // Write the file only if correctly opened
        {
        }
        else //If impossible to write in the file, say it
        {
            cout << "ERROR: Impossible to open the output file" << endl;
        }

        outputFlux << "Run_simulation" << " " << "Memory_type" << " " << "Time_length_run" << " " << "Days_wo_food" << " " << "Foraging_success_cum" << " " << "Foraging_success_cum_target_only" << " " << "Foraging_success_tot" << " " << "Foraging_success_tot_target_only" << " " << "Percent_tree_visited" << " " <<
        "Percent_targets_with_food" << " " << "Number_trees_rare_species" << " " << "Rate_trees_rare_eaten" << " " << "Proba_competition" << " " << "Working_memory_length_time_unit" << " " << "Time_not_to_return" << " " << "Vision_radius" << " " << "Length_cycle" << " " << "Length_fruiting" << " " << "Number_species" << " " << "Number_individuals_per_species" << " " <<
        "Intra_species_synchrony" << " " << "Percentage_low_synchrony" << " " << "Distance_min_start" << " " << "Sd_min_start" << " " << "Distance_min_end" << " " << "Sd_min_end" << " " << "Weight_for_longterm_memory" << endl;

        outputFlux.close();

        for(int p=0; p<(intraSpeciesSynchronyValuesInit.size()); p++){

            double proportionLowSynchrony((intraSpeciesSynchronyValuesInit[p] - 0.1)/0.8);

                //--------------
                //INITIALISATION
                //--------------

                /*Time-associated variables*/
                double timerInit(-cycleLengthInit);
                double previousTimerInit(-cycleLengthInit/4);
                double coefficientDistanceToTimeInit((double) 1/((double) 1/cycleLengthInit*(double) mapSizeInit*(double) mapSizeInit/(double) visionRadiusInit/2.0));

                /*Species-level variables*/

                vector<int> startFruitingDateSpeciesOperationalInit;
                vector<int> endFruitingDateSpeciesOperationalInit;

                vector<int> startFruitingDateSpeciesAverageInitResult;
                vector<int> startFruitingDateSpeciesStdevInitResult;

                vector<double> intraSpeciesSynchronyInit;

                /*Tree-level variables*/
                vector<double> xCoordinateTreeInit;
                vector<double> yCoordinateTreeInit;
                vector<int> startFruitingDateTreeInit;
                vector<int> endFruitingDateTreeInit;
                vector<int> speciesTreeInit;
                vector<double> foodQuantityTreeInit;

                /*Agent-level variables*/
                double xCoordinateAgentInit(disULocTree(gen));
                double yCoordinateAgentInit(disULocTree(gen));
                double workingMemoryTimeInit(fruitingLengthInit/2);
                double weightInFavourOfLongtermKnowledgeInit(1.0); //Only for associative/chronological memory, should be comprised between [-1;1]. Used for averaging between observed (=estimated) fruiting, and Longterm knowledge on fruiting dates

                /* Environment */
                double earliestStartFruitingDateInit(cycleLengthInit);
                double latestEndFruitingDateInit(0.0);

                std::uniform_real_distribution<double> disU01(0, 1); //Uniform distrib
                std::bernoulli_distribution bernouilliDisInit(0.5);

                int speciesNumberWithLowSynchrony(0);

                /*Setting the environment: resource distribution and tree features*/

                for (int i=0; i<numberSpeciesInit; i++){

                        if(speciesNumberWithLowSynchrony/ (double) numberSpeciesInit < proportionLowSynchrony){
                            intraSpeciesSynchronyInit.push_back(0.1);
                            speciesNumberWithLowSynchrony+=1;
                        }
                        else{
                            intraSpeciesSynchronyInit.push_back(0.9);
                        }

                        //ADDED ON THE 01/02/2021
                        //If heterogeneous distribution
                        std::vector<int> patchesCoordinatesX(numberPatch,0);
                        std::vector<int> patchesCoordinatesY(numberPatch,0);
                        for(int p=0; p<patchesCoordinatesX.size(); p++){
                            patchesCoordinatesX[p]=disULocPatch(gen);
                            patchesCoordinatesY[p]=disULocPatch(gen);
                        }
                        int counterTree=0;
                        int refPatch=0;
                        std::normal_distribution<double> disNormalPatchX(patchesCoordinatesX[0],mapSizeInit/dispersion);//
                        std::normal_distribution<double> disNormalPatchY(patchesCoordinatesY[0],mapSizeInit/dispersion);//
                        //////////////////////////

                        //int speciesToSplit(0);
                        vector<int> dateStartVector;
                        vector<int> dateEndVector;
                        int speciesStartDate(disUDateSpecies(gen));

                        for (int j=0; j<numberIndividualsPerSpeciesInit; j++){

                            //New method triangular: sum of uniform distributions
                            int dateTransitory(speciesStartDate +  (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit/2*(disU01(gen)+disU01(gen)));

                            startFruitingDateTreeInit.push_back(dateTransitory);
                            endFruitingDateTreeInit.push_back(dateTransitory + fruitingLengthInit);

                            //Food quantity
                            foodQuantityTreeInit.push_back(0);//No food at first (will be updated afterwards)

                            //Coordinates

                            /////////// ADDED ON THE 01/02/2021
                            if(isHeterogeneous=="YES"){
                                counterTree+=1;
                                if(counterTree>=floor(numberIndividualsPerSpeciesInit/(double) numberPatch)){
                                    if(refPatch!=patchesCoordinatesX.size()-1){
                                        refPatch+=1;
                                        counterTree=0;
                                        disNormalPatchX.param(std::normal_distribution<double>::param_type(patchesCoordinatesX[refPatch],mapSizeInit/dispersion));//
                                        disNormalPatchY.param(std::normal_distribution<double>::param_type(patchesCoordinatesY[refPatch],mapSizeInit/dispersion));//
                                    }
                                }
                                double locX(disNormalPatchX(gen));
                                while(locX < 0 || locX > mapSizeInit){
                                    locX=disNormalPatchX(gen);
                                }
                                double locY(disNormalPatchY(gen));
                                while(locY < 0 || locY > mapSizeInit){
                                    locY=disNormalPatchY(gen);
                                }
                                xCoordinateTreeInit.push_back(locX);
                                yCoordinateTreeInit.push_back(locY);
                            }
                            else{
                            ////////////////////////////////////
                                xCoordinateTreeInit.push_back(disULocTree(gen));
                                yCoordinateTreeInit.push_back(disULocTree(gen));
                            }

                            //Update the sum for mean and var calculation
                            dateStartVector.push_back(dateTransitory);
                            dateEndVector.push_back(dateTransitory+fruitingLengthInit);

                            speciesTreeInit.push_back(i);
                        }

                        /*Calculating the operational starting date*/
                        //std::sort (dateStartVector.begin(), dateStartVector.end());//Ranging in increasing order starting dates
                        startFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate + 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));//Correction adding bc subsampling triggers that we never reach this extreme value. Hence Simon conducted simulation to see what correction we have to apply when sampling in triangular distrib, the extreme values we obtained on average.//startFruitingDateSpeciesOperationalInit.push_back(dateStartVector[0]);

                        /*Calculating the operational ending date*/
                        //std::sort (dateEndVector.rbegin(), dateEndVector.rend());//Ranging in decreasing order ending dates
                        endFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate  + (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit - 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));//endFruitingDateSpeciesOperationalInit.push_back(dateEndVector[0]);

                        if(fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit) < earliestStartFruitingDateInit
                           //&& i > 0
                           ){//For the associative scenario it is necessary to get rid of the period with no fruits or no fruit predictable
                            //NOT ANYMORE BECAUSE SIMULATING MANY CYCLES: therefore eliminating when only sp1
                                earliestStartFruitingDateInit=fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                        }
                        if(fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit) > latestEndFruitingDateInit){
                                latestEndFruitingDateInit=fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                        }

                        /*Extract the finally obtained mean and var for all species*/
                        //Calculate sum, mean, square sum and stdv
                        double sum = std::accumulate(dateStartVector.begin(), dateStartVector.end(), 0.0);//arguments are: first value, last value, initialisation value
                        double mean = sum / dateStartVector.size();
                        double sq_sum = std::inner_product(dateStartVector.begin(), dateStartVector.end(), dateStartVector.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
                        double stdev = std::sqrt(sq_sum / dateStartVector.size() - mean * mean);

                        startFruitingDateSpeciesAverageInitResult.push_back(mean);
                        startFruitingDateSpeciesStdevInitResult.push_back(stdev);
                    }

                    earliestStartFruitingDateInit=max(0.0, earliestStartFruitingDateInit);
                    latestEndFruitingDateInit=min((double) cycleLengthInit, latestEndFruitingDateInit);

                /*Determiner the associative strength between start-start events and start-end events*/
                //and /*Create the matrix for associative model: first col is the species, second is the distance to the date, and third col characterises which date is to use 0=start, 1=end*/

                //Create the result matrix
                int distanceStartStartInit[numberSpeciesInit][numberSpeciesInit];
                int distanceStartEndInit[numberSpeciesInit][numberSpeciesInit];
                int distanceEndStartInit[numberSpeciesInit][numberSpeciesInit];
                int distanceEndEndInit[numberSpeciesInit][numberSpeciesInit];
                int matrixInferenceAssociationInit[numberSpeciesInit][3];
                std::fill(*matrixInferenceAssociationInit, *matrixInferenceAssociationInit + 3*numberSpeciesInit, cycleLengthInit); //Initialise matrix with only cycleLengthInit

                //Create the min vectors of distance
                vector<double> vectorMinDistanceStart;
                vector<double> vectorMinDistanceEnd;

                //Loop to fill the result matrices and calculate alongside the minimum: first rows are always start-start then start-end
                for (int j=0; j<numberSpeciesInit; j++){

                    for (int k=0; k<numberSpeciesInit; k++){
                    if(j!=k){
                        //Updating the results matrices
                            distanceStartStartInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));
                            distanceStartEndInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]-fruitingLengthInit),cycleLengthInit));
                            distanceEndStartInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]+fruitingLengthInit),cycleLengthInit));
                            distanceEndEndInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));

                            vectorMinDistanceStart.push_back(min(distanceStartStartInit[j][k],
                                                                 distanceStartEndInit[j][k]
                                                                 )
                                                             );
                            vectorMinDistanceEnd.push_back(min(distanceEndStartInit[j][k],
                                                               distanceEndEndInit[j][k]
                                                               )
                                                           );
                        if(distanceStartStartInit[j][k]==min(distanceStartStartInit[j][k],distanceStartEndInit[j][k])){
                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                               (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] >=
                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                            ){
                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                       (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] ==
                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                                    ){
                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                            ){
                                               matrixInferenceAssociationInit[j][0]=k;
                                               matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                               matrixInferenceAssociationInit[j][2]=0;
                                            }
                                            else{
                                                //Do nothing
                                            }
                                    }
                                    else{
                                         matrixInferenceAssociationInit[j][0]=k;
                                         matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                         matrixInferenceAssociationInit[j][2]=0;
                                    }
                            }
                            else{
                                //Do nothing
                            }
                        }
                        else{
                           if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                               (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] >=
                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                            ){
                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                       (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] ==
                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                                    ){
                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                            ){
                                               matrixInferenceAssociationInit[j][0]=k;
                                               matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                               matrixInferenceAssociationInit[j][2]=1;
                                            }
                                            else{
                                                //Do nothing
                                            }
                                    }
                                    else{
                                         matrixInferenceAssociationInit[j][0]=k;
                                         matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                         matrixInferenceAssociationInit[j][2]=1;
                                    }
                            }
                        }
                    }
                }

                ////Checking what are the chosen referential species
                //cout << matrixInferenceAssociationInit[j][0] << " " << intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]] << endl;
            }

            //Calculate sum, mean, square sum and stdv
            double sum = std::accumulate(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), 0.0);//arguments are: first value, last value, initialisation value
            double mean = sum / vectorMinDistanceStart.size();
            double sq_sum = std::inner_product(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), vectorMinDistanceStart.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
            double stdev = std::sqrt(sq_sum / vectorMinDistanceStart.size() - mean * mean);

            double averageMinDistanceStartInit(mean);
            double sdMinDistanceStartInit(stdev);

            sum = std::accumulate(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), 0.0);//arguments are: first value, last value, initialisation value
            mean = sum / vectorMinDistanceEnd.size();
            sq_sum = std::inner_product(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), vectorMinDistanceEnd.begin(), 0.0);//arguments are: first value to End with, last value, first value to match to the other one, initialisation value
            stdev = std::sqrt(sq_sum / vectorMinDistanceEnd.size() - mean * mean);

            double averageMinDistanceEndInit(mean);
            double sdMinDistanceEndInit(stdev);

            //--------------
            //Running models
            //--------------

            /*Run model: NULL*/

            cout << "I work until line setting loop for model: NULL" << endl;

            nullModel(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double workingMemoryTime=
                    workingMemoryTimeInit, //max(2.0, (double)fruitingLengthInit/20),
                    //double timeLapseNotToReturn=
                    timeLapseNotToReturnInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Null_model_telescopic",
                    // bool extendibleArms,
                    true,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,

                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,

                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0
            );

            /*Run model: NULL WORKING MEMORY*/
            cout << "I work until line setting loop for model: NULL WORKING MEMORY" << endl;

            nullModelWorkingMemory(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double workingMemoryTime=
                    workingMemoryTimeInit,
                    //double timeLapseNotToReturn=
                    timeLapseNotToReturnInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Null_working_memory_model_telescopic",
                    // bool extendibleArms,
                    true,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,

                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,

                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0

            );

            /*Run model: CHRONOLOGICAL*/
            cout << "I work until line setting loop for model: CHRONOLOGICAL" << endl;

            chronologicalMemoryModel(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double workingMemoryTime=
                    workingMemoryTimeInit,
                    //double timeLapseNotToReturn=
                    timeLapseNotToReturnInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Chronological_model_telescopic",
                    // bool extendibleArms,
                    true,
                    //double weightInFavourOfLongtermKnowledge,
                    weightInFavourOfLongtermKnowledgeInit,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,
                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,
                    //std::vector<int> startFruitingDateSpeciesOperational=
                    startFruitingDateSpeciesOperationalInit,
                    //std::vector<int> endFruitingDateSpeciesOperational=
                    endFruitingDateSpeciesOperationalInit,
                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0

            );

            /*Run model: ASSOCIATIVE*/
            cout << "I work until line setting loop for model: ASSOCIATIVE" << endl;

            associativeMemoryModel(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double workingMemoryTime=
                    workingMemoryTimeInit,
                    //double timeLapseNotToReturn=
                    timeLapseNotToReturnInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Associative_model_telescopic",
                    // bool extendibleArms,
                    true,
                    //double weightInFavourOfLongtermKnowledge,
                    weightInFavourOfLongtermKnowledgeInit,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,

                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,
                    //int matrixInferenceAssociation,
                    *matrixInferenceAssociationInit,

                    //std::vector<int> startFruitingDateSpeciesOperational=
                    startFruitingDateSpeciesOperationalInit,
                    //std::vector<int> endFruitingDateSpeciesOperational=
                    endFruitingDateSpeciesOperationalInit,

                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0

            );


            /*Run model: OMNISCIENT*/
            cout << "I work until line setting loop for model: OMNISCIENT" << endl;

            omniscientModel(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Omniscient_model_telescopic",
                    // bool extendibleArms,
                    true,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,

                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,

                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0
            );

        }
    }


// Set back spatial environment to heterogeneous
isHeterogeneous="YES";


    ////////////////////////////////////////////////////////////////////
    // SCENARIO 1: Temporally homogeneous/Spatially heterogeneous environment with changing synchrony
    ////////////////////////////////////////////////////////////////////

    #pragma omp parallel for
    for (int s=0; s<300; s++){

        std::string pathOfTheFileOutputInitMdf(pathOfTheFileOutputInit);
        pathOfTheFileOutputInitMdf.append(to_string(s));
        pathOfTheFileOutputInitMdf.append("_IS_spatialHeterogeneous.txt");
        std::ofstream outputFlux;
        outputFlux.open(pathOfTheFileOutputInitMdf.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
        if(outputFlux)  // Write the file only if correctly opened
        {
        }
        else //If impossible to write in the file, say it
        {
            cout << "ERROR: Impossible to open the output file" << endl;
        }

        outputFlux << "Run_simulation" << " " << "Memory_type" << " " << "Time_length_run" << " " << "Days_wo_food" << " " << "Foraging_success_cum" << " " << "Foraging_success_cum_target_only" << " " << "Foraging_success_tot" << " " << "Foraging_success_tot_target_only" << " " << "Percent_tree_visited" << " " <<
        "Percent_targets_with_food" << " " << "Number_trees_rare_species" << " " << "Rate_trees_rare_eaten" << " " << "Proba_competition" << " " << "Working_memory_length_time_unit" << " " << "Time_not_to_return" << " " << "Vision_radius" << " " << "Length_cycle" << " " << "Length_fruiting" << " " << "Number_species" << " " << "Number_individuals_per_species" << " " <<
        "Intra_species_synchrony" << " " << "Percentage_low_synchrony" << " " << "Distance_min_start" << " " << "Sd_min_start" << " " << "Distance_min_end" << " " << "Sd_min_end" << " " << "Weight_for_longterm_memory" << endl;

        outputFlux.close();
        for(int iSS=0; iSS<intraSpeciesSynchronyValuesInit.size(); iSS++){

            //--------------
            //INITIALISATION
            //--------------

            /*Time-associated variables*/
            double timerInit(-cycleLengthInit);
            double previousTimerInit(-cycleLengthInit);
            double coefficientDistanceToTimeInit((double) 1/((double) 1/cycleLengthInit*(double) mapSizeInit*(double) mapSizeInit/(double) visionRadiusInit/2.0));
            //cout << 1/coefficientDistanceToTimeInit << endl;

            // This is set up by making the hypothesis that during the cycle all of the trees can be visited. We assume that, during the cycle the individual hence travels a distance d
            // which makes it cover d*visionRadius surface (approximation not considering overlap, hence overestimation). Applying the "produit en croix", we find that d can't be more than 1/2*mapSize**2/visionRadius. Hence this speed if defining working memory as FruitingLength/2

            /*Species-level variables*/

            vector<int> startFruitingDateSpeciesOperationalInit;
            vector<int> endFruitingDateSpeciesOperationalInit;

            vector<int> startFruitingDateSpeciesAverageInitResult;
            vector<int> startFruitingDateSpeciesStdevInitResult;

            vector<double> intraSpeciesSynchronyInit;

            /*Tree-level variables*/
            vector<double> xCoordinateTreeInit;
            vector<double> yCoordinateTreeInit;
            vector<int> startFruitingDateTreeInit;
            vector<int> endFruitingDateTreeInit;
            vector<int> speciesTreeInit;
            vector<double> foodQuantityTreeInit;

            /*Agent-level variables*/
            double xCoordinateAgentInit(disULocTree(gen));
            double yCoordinateAgentInit(disULocTree(gen));
            double workingMemoryTimeInit(15.0);
            double weightInFavourOfLongtermKnowledgeInit(1.0); //Only for associative/chronological memory, should be comprised between [-1;1]. Used for averaging between observed (=estimated) fruiting, and Longterm knowledge on fruiting dates

            /* Environment */
            double earliestStartFruitingDateInit(cycleLengthInit);
            double latestEndFruitingDateInit(0.0);

            std::uniform_real_distribution<double> disU01(0, 1); //Uniform distrib
            std::bernoulli_distribution bernouilliDisInit(0.5);


            //if(s==0){
                ////CHECK DISTRIBUTION
//                std::string pathOfTheFileOutputInitMdf(pathOfTheFileOutputInit);
//                pathOfTheFileOutputInitMdf.append(to_string(s));
//                pathOfTheFileOutputInitMdf.append("_check_distribution.txt");
//                std::ofstream outputFluxCheck;
//                outputFluxCheck.open(pathOfTheFileOutputInitMdf.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
//                if(outputFluxCheck)  // Write the file only if correctly opened
//                {
//                }
//                else //If impossible to write in the file, say it
//                {
//                    cout << "ERROR: Impossible to open the output file" << endl;
//                }
//
//                outputFluxCheck << "Species " << "x " << "y" << endl;
            //}

            /*Setting the environment: resource distribution and tree features*/
            for (int i=0; i<numberSpeciesInit; i++){

                    intraSpeciesSynchronyInit.push_back(intraSpeciesSynchronyValuesInit[iSS]);
                    //int speciesToSplit(0);
                    vector<int> dateStartVector;
                    vector<int> dateEndVector;
                        int speciesStartDate(disUDateSpecies(gen));

                        //ADDED ON THE 01/02/2021
                        //If heterogeneous distribution
                        std::vector<int> patchesCoordinatesX(numberPatch,0);
                        std::vector<int> patchesCoordinatesY(numberPatch,0);
                        for(int p=0; p<patchesCoordinatesX.size(); p++){
                            patchesCoordinatesX[p]=disULocPatch(gen);
                            //cout << patchesCoordinatesX[p] << endl;
                            patchesCoordinatesY[p]=disULocPatch(gen);
                            //cout << patchesCoordinatesY[p] << endl;
                        }
                        int counterTree=0;
                        int refPatch=0;
                        std::normal_distribution<double> disNormalPatchX(patchesCoordinatesX[0],mapSizeInit/dispersion);//
                        std::normal_distribution<double> disNormalPatchY(patchesCoordinatesY[0],mapSizeInit/dispersion);//
                        //////////////////////////

                        for (int j=0; j<numberIndividualsPerSpeciesInit; j++){
                            //New method triangular: sum of uniform distributions (NOTE: in the article for formula, we are using the Tm (the max point), and not the TS, therefore you can see that it sums the uniform - 1 !!!!!!!!
                            int dateTransitory(speciesStartDate +  (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit/2*(disU01(gen)+disU01(gen)));

                            startFruitingDateTreeInit.push_back(dateTransitory);
                            endFruitingDateTreeInit.push_back(dateTransitory + fruitingLengthInit);

                            //Food quantity
                            foodQuantityTreeInit.push_back(0);//No food at first (will be updated afterwards)

                            //Coordinates
                            //cout << floor(numberIndividualsPerSpeciesInit/(double) numberPatch) << endl;
                            /////////// ADDED ON THE 01/02/2021
                            if(isHeterogeneous=="YES"){
                                //cout << counterTree << endl;
                                counterTree+=1;
                                //cout << patchesCoordinatesX.size()-1 << endl;
                                if(counterTree>=floor(numberIndividualsPerSpeciesInit/(double) numberPatch)){
                                    if(refPatch!=patchesCoordinatesX.size()-1){
                                        refPatch+=1;
                                        //cout << refPatch << endl;
                                        counterTree=0;
                                        //cout << patchesCoordinatesX[refPatch] << " " << patchesCoordinatesY[refPatch] << endl;
                                        disNormalPatchX.param(std::normal_distribution<double>::param_type(patchesCoordinatesX[refPatch],mapSizeInit/dispersion));//
                                        disNormalPatchY.param(std::normal_distribution<double>::param_type(patchesCoordinatesY[refPatch],mapSizeInit/dispersion));//
                                    }
                                }
                                double locX(disNormalPatchX(gen));
                                while(locX < 0 || locX > mapSizeInit){
                                    locX=disNormalPatchX(gen);
                                    //cout << locX << endl;
                                }
                                double locY(disNormalPatchY(gen));
                                while(locY < 0 || locY > mapSizeInit){
                                    locY=disNormalPatchY(gen);
                                    //cout << locY << endl;
                                }
                                xCoordinateTreeInit.push_back(locX);
                                yCoordinateTreeInit.push_back(locY);
                            }
                            else{
                            ////////////////////////////////////
                                xCoordinateTreeInit.push_back(disULocTree(gen));
                                yCoordinateTreeInit.push_back(disULocTree(gen));
                            }
//                            if(s==0){
//                                outputFluxCheck << i << " " << xCoordinateTreeInit[xCoordinateTreeInit.size()-1] << " " << yCoordinateTreeInit[yCoordinateTreeInit.size()-1] << endl;
//                            }
                            //Update vector for mean and var calculation
                            dateStartVector.push_back(dateTransitory);
                            dateEndVector.push_back(dateTransitory+fruitingLengthInit);

                            speciesTreeInit.push_back(i);
                        }

                    /*Calculating the operational starting date*/
                    //std::sort (dateStartVector.begin(), dateStartVector.end());//Ranging in increasing order starting dates
                    startFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate + 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));//Correction adding bc subsampling triggers that we never reach this extreme value. Hence Simon conducted simulation to see what correction we have to apply when sampling in triangular distrib, the extreme values we obtained on average.

                    /*Calculating the operational ending date*/
                    //std::sort (dateEndVector.rbegin(), dateEndVector.rend());//Ranging in decreasing order ending dates
                    endFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate  + (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit - 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));

                    if(fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit) < earliestStartFruitingDateInit
                       ){
                            earliestStartFruitingDateInit=fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                    }
                    if(fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit) > latestEndFruitingDateInit
                       ){
                            latestEndFruitingDateInit=fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                    }

                    /*Extract the finally obtained mean and var for all species*/
                    //Calculate sum, mean, square sum and stdv
                    double sum = std::accumulate(dateStartVector.begin(), dateStartVector.end(), 0.0);//arguments are: first value, last value, initialisation value
                    double mean = sum / dateStartVector.size();
                    double sq_sum = std::inner_product(dateStartVector.begin(), dateStartVector.end(), dateStartVector.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
                    double stdev = std::sqrt(sq_sum / dateStartVector.size() - mean * mean);

                    if(intraSpeciesSynchronyInit[i] > 0){
                    }
                    else{
                        mean=cycleLengthInit/2.0;
                        stdev=sqrt(cycleLengthInit*cycleLengthInit/12.0);
                    }
                    startFruitingDateSpeciesAverageInitResult.push_back(mean);
                    startFruitingDateSpeciesStdevInitResult.push_back(stdev);
                }

                earliestStartFruitingDateInit=max(0.0, earliestStartFruitingDateInit);
                latestEndFruitingDateInit=min((double) cycleLengthInit, latestEndFruitingDateInit);

                /*Determinate the associative strength between start-start events and start-end events*/
                //and /*Create the matrix for associative model: first col is the species, second is the distance to the date, and third col characterises which date is to use 0=start, 1=end*/

                //Create the result matrix
                int distanceStartStartInit[numberSpeciesInit][numberSpeciesInit];
                int distanceStartEndInit[numberSpeciesInit][numberSpeciesInit];
                int distanceEndStartInit[numberSpeciesInit][numberSpeciesInit];
                int distanceEndEndInit[numberSpeciesInit][numberSpeciesInit];
                int matrixInferenceAssociationInit[numberSpeciesInit][3];
                std::fill(*matrixInferenceAssociationInit, *matrixInferenceAssociationInit + 3*numberSpeciesInit, cycleLengthInit); //Initialise matrix with only cycleLengthInit

                //Create the min vectors of distance
                vector<double> vectorMinDistanceStart;
                vector<double> vectorMinDistanceEnd;

                //Loop to fill the result matrices and calculate alongside the minimum: first rows are always start-start then start-end
                for (int j=0; j<numberSpeciesInit; j++){

                    for (int k=0; k<numberSpeciesInit; k++){
                    if(j!=k){
                        //Updating the results matrices
                            distanceStartStartInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));
                            distanceStartEndInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]-fruitingLengthInit),cycleLengthInit));
                            distanceEndStartInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]+fruitingLengthInit),cycleLengthInit));
                            distanceEndEndInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));

                            vectorMinDistanceStart.push_back(min(distanceStartStartInit[j][k],
                                                                 distanceStartEndInit[j][k]
                                                                 )
                                                             );
                            vectorMinDistanceEnd.push_back(min(distanceEndStartInit[j][k],
                                                               distanceEndEndInit[j][k]
                                                               )
                                                           );
                        if(distanceStartStartInit[j][k]==min(distanceStartStartInit[j][k],distanceStartEndInit[j][k])){
                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                               (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] >=
                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                            ){
                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                       (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] ==
                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]

                                    ){
                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                            ){
                                               matrixInferenceAssociationInit[j][0]=k;
                                               matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                               matrixInferenceAssociationInit[j][2]=0;
                                            }
                                            else{
                                                //Do nothing
                                            }
                                    }
                                    else{
                                         matrixInferenceAssociationInit[j][0]=k;
                                         matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                         matrixInferenceAssociationInit[j][2]=0;
                                    }
                            }
                            else{
                                //Do nothing
                            }
                        }
                        else{
                           if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                               (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] >=
                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                            ){
                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                       (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] ==
                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]

                                    ){
                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                            ){
                                               matrixInferenceAssociationInit[j][0]=k;
                                               matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                               matrixInferenceAssociationInit[j][2]=1;
                                            }
                                            else{
                                                //Do nothing
                                            }
                                    }
                                    else{
                                         matrixInferenceAssociationInit[j][0]=k;
                                         matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                         matrixInferenceAssociationInit[j][2]=1;
                                    }
                            }
                        }
                    }
                }
            }
//            outputFluxCheck.close();

            //Calculate sum, mean, square sum and stdv
            double sum = std::accumulate(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), 0.0);//arguments are: first value, last value, initialisation value
            double mean = sum / vectorMinDistanceStart.size();
            double sq_sum = std::inner_product(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), vectorMinDistanceStart.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
            double stdev = std::sqrt(sq_sum / vectorMinDistanceStart.size() - mean * mean);

            double averageMinDistanceStartInit(mean);
            double sdMinDistanceStartInit(stdev);

            sum = std::accumulate(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), 0.0);//arguments are: first value, last value, initialisation value
            mean = sum / vectorMinDistanceEnd.size();
            sq_sum = std::inner_product(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), vectorMinDistanceEnd.begin(), 0.0);//arguments are: first value to End with, last value, first value to match to the other one, initialisation value
            stdev = std::sqrt(sq_sum / vectorMinDistanceEnd.size() - mean * mean);

            double averageMinDistanceEndInit(mean);
            double sdMinDistanceEndInit(stdev);


        //--------------
        //Running models
        //--------------

       /* Run model: NULL*/

        cout << "I work until line setting loop for model: NULL" << endl;

        nullModel(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double workingMemoryTime=
                workingMemoryTimeInit, //max(2.0, (double)fruitingLengthInit/20),
                //double timeLapseNotToReturn=
                timeLapseNotToReturnInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Null_model_telescopic",
                // bool extendibleArms,
                true,

                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0

        );

        /*Run model: NULL WORKING MEMORY*/
        cout << "I work until line setting loop for model: NULL WORKING MEMORY" << endl;

        nullModelWorkingMemory(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double workingMemoryTime=
                workingMemoryTimeInit,
                //double timeLapseNotToReturn=
                timeLapseNotToReturnInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Null_working_memory_model_telescopic",
                // bool extendibleArms,
                true,

                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0
        );

        /*Run model: CHRONOLOGICAL*/
        cout << "I work until line setting loop for model: CHRONOLOGICAL" << endl;

        chronologicalMemoryModel(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double workingMemoryTime=
                workingMemoryTimeInit,
                //double timeLapseNotToReturn=
                timeLapseNotToReturnInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Chronological_model_telescopic",
                // bool extendibleArms,
                true,
                //double weightInFavourOfLongtermKnowledge,
                weightInFavourOfLongtermKnowledgeInit,

                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,
                //std::vector<int> startFruitingDateSpeciesOperational=
                startFruitingDateSpeciesOperationalInit,
                //std::vector<int> endFruitingDateSpeciesOperational=
                endFruitingDateSpeciesOperationalInit,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0
        );


        /*Run model: ASSOCIATIVE*/
        cout << "I work until line setting loop for model: ASSOCIATIVE" << endl;

        associativeMemoryModel(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //-cycleLengthInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double workingMemoryTime=
                workingMemoryTimeInit,
                //double timeLapseNotToReturn=
                timeLapseNotToReturnInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Associative_model_telescopic",
                // bool extendibleArms,
                true,
                //double weightInFavourOfLongtermKnowledge,
                weightInFavourOfLongtermKnowledgeInit,


                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,
                //int matrixInferenceAssociation,
                *matrixInferenceAssociationInit,
                //std::vector<int> startFruitingDateSpeciesOperational=
                startFruitingDateSpeciesOperationalInit,
                //std::vector<int> endFruitingDateSpeciesOperational=
                endFruitingDateSpeciesOperationalInit,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0
        );

        /*Run model: OMNISCIENT*/
        cout << "I work until line setting loop for model: OMNISCIENT" << endl;

        omniscientModel(
                //Run ID
                //int step=
                s,

                //Output variable
                //string pathOfTheFileOutput=
                pathOfTheFileOutputInitMdf,

                //Time
                //double timer=
                timerInit,
                //double previousTimer=
                previousTimerInit,
                //double cycleLength=
                cycleLengthInit,

                //Agent ability
                //double visionRadius=
                visionRadiusInit,
                //double coefficientDistanceToTime=
                coefficientDistanceToTimeInit,
                //string memoryType=
                "Omniscient_model_telescopic",
                // bool extendibleArms,
                true,

                //Agent location and success
                //double xCoordinateAgent=
                xCoordinateAgentInit,
                //double yCoordinateAgent=
                yCoordinateAgentInit,

                //Initial environmental condition

                //On Species
                //int numberSpecies=
                numberSpeciesInit,
                //int numberIndividualsPerSpecies=
                numberIndividualsPerSpeciesInit,
                //vector<int>quantityFoodMaxTreeForSpecies,
                quantityFoodMaxTreeForSpeciesInit,
                //int numberTreesRareSpeciesInit,
                numberTreesRareSpeciesInit,
                //double fruitingLength=
                fruitingLengthInit,
                //int earliestStartFruitingDate=
                earliestStartFruitingDateInit,
                //int latestEndFruitingDate=
                latestEndFruitingDateInit,
                //double averageMinDistanceStart=
                averageMinDistanceStartInit,
                //double sdMinDistanceStart,
                sdMinDistanceStartInit,
                //double averageMinDistanceEnd=
                averageMinDistanceEndInit,
                //double sdMinDistanceEnd,
                sdMinDistanceEndInit,
                //std::vector<double> intraSpeciesSynchrony=
                intraSpeciesSynchronyInit,
                //double percentSpeciesWithLowSynchrony=
                1.0,

                //On Tree
                //std::vector<double> xCoordinateTree=
                xCoordinateTreeInit,
                //std::vector<double> yCoordinateTree=
                yCoordinateTreeInit,
                //std::vector<double> foodQuantityTree=
                foodQuantityTreeInit,
                //std::vector<int> speciesTree=
                speciesTreeInit,
                //std::vector<int> endFruitingDateTree=
                endFruitingDateTreeInit,
                //std::vector<int> startFruitingDateTree=
                startFruitingDateTreeInit,
                //double probaTreeAlreadyDepleted,
                0.0
        );
    }
}

////////////////////////////////////////////////////////////////
// SCENARIO 2b: Temporally heterogeneous/spatially heterogeneous environment with change in proportion of little/highly synchronous species
////////////////////////////////////////////////////////////////

    #pragma omp parallel for
    for (int s=0; s<300; s++){
        cout << s << endl;
        std::string pathOfTheFileOutputInitMdf(pathOfTheFileOutputInit);
        pathOfTheFileOutputInitMdf.append(to_string(s));
        pathOfTheFileOutputInitMdf.append("_HE_spatialHeterogeneous.txt");
        std::ofstream outputFlux;
        outputFlux.open(pathOfTheFileOutputInitMdf.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
        if(outputFlux)  // Write the file only if correctly opened
        {
        }
        else //If impossible to write in the file, say it
        {
            cout << "ERROR: Impossible to open the output file" << endl;
        }

        outputFlux << "Run_simulation" << " " << "Memory_type" << " " << "Time_length_run" << " " << "Days_wo_food" << " " << "Foraging_success_cum" << " " << "Foraging_success_cum_target_only" << " " << "Foraging_success_tot" << " " << "Foraging_success_tot_target_only" << " " << "Percent_tree_visited" << " " <<
        "Percent_targets_with_food" << " " << "Number_trees_rare_species" << " " << "Rate_trees_rare_eaten" << " " << "Proba_competition" << " " << "Working_memory_length_time_unit" << " " << "Time_not_to_return" << " " << "Vision_radius" << " " << "Length_cycle" << " " << "Length_fruiting" << " " << "Number_species" << " " << "Number_individuals_per_species" << " " <<
        "Intra_species_synchrony" << " " << "Percentage_low_synchrony" << " " << "Distance_min_start" << " " << "Sd_min_start" << " " << "Distance_min_end" << " " << "Sd_min_end" << " " << "Weight_for_longterm_memory" << endl;

        outputFlux.close();

        for(int p=0; p<(intraSpeciesSynchronyValuesInit.size()); p++){

            double proportionLowSynchrony((intraSpeciesSynchronyValuesInit[p] - 0.1)/0.8);

                //--------------
                //INITIALISATION
                //--------------

                /*Time-associated variables*/
                double timerInit(-cycleLengthInit);
                double previousTimerInit(-cycleLengthInit/4);
                double coefficientDistanceToTimeInit((double) 1/((double) 1/cycleLengthInit*(double) mapSizeInit*(double) mapSizeInit/(double) visionRadiusInit/2.0));

                /*Species-level variables*/

                vector<int> startFruitingDateSpeciesOperationalInit;
                vector<int> endFruitingDateSpeciesOperationalInit;

                vector<int> startFruitingDateSpeciesAverageInitResult;
                vector<int> startFruitingDateSpeciesStdevInitResult;

                vector<double> intraSpeciesSynchronyInit;

                /*Tree-level variables*/
                vector<double> xCoordinateTreeInit;
                vector<double> yCoordinateTreeInit;
                vector<int> startFruitingDateTreeInit;
                vector<int> endFruitingDateTreeInit;
                vector<int> speciesTreeInit;
                vector<double> foodQuantityTreeInit;

                /*Agent-level variables*/
                double xCoordinateAgentInit(disULocTree(gen));
                double yCoordinateAgentInit(disULocTree(gen));
                double workingMemoryTimeInit(fruitingLengthInit/2);
                double weightInFavourOfLongtermKnowledgeInit(1.0); //Only for associative/chronological memory, should be comprised between [-1;1]. Used for averaging between observed (=estimated) fruiting, and Longterm knowledge on fruiting dates

                /* Environment */
                double earliestStartFruitingDateInit(cycleLengthInit);
                double latestEndFruitingDateInit(0.0);

                std::uniform_real_distribution<double> disU01(0, 1); //Uniform distrib
                std::bernoulli_distribution bernouilliDisInit(0.5);

                int speciesNumberWithLowSynchrony(0);

                /*Setting the environment: resource distribution and tree features*/

                for (int i=0; i<numberSpeciesInit; i++){

                        if(speciesNumberWithLowSynchrony/ (double) numberSpeciesInit < proportionLowSynchrony){
                            intraSpeciesSynchronyInit.push_back(0.1);
                            speciesNumberWithLowSynchrony+=1;
                        }
                        else{
                            intraSpeciesSynchronyInit.push_back(0.9);
                        }

                        //ADDED ON THE 01/02/2021
                        //If heterogeneous distribution
                        std::vector<int> patchesCoordinatesX(numberPatch,0);
                        std::vector<int> patchesCoordinatesY(numberPatch,0);
                        for(int p=0; p<patchesCoordinatesX.size(); p++){
                            patchesCoordinatesX[p]=disULocPatch(gen);
                            patchesCoordinatesY[p]=disULocPatch(gen);
                        }
                        int counterTree=0;
                        int refPatch=0;
                        std::normal_distribution<double> disNormalPatchX(patchesCoordinatesX[0],mapSizeInit/dispersion);//
                        std::normal_distribution<double> disNormalPatchY(patchesCoordinatesY[0],mapSizeInit/dispersion);//
                        //////////////////////////

                        //int speciesToSplit(0);
                        vector<int> dateStartVector;
                        vector<int> dateEndVector;
                        int speciesStartDate(disUDateSpecies(gen));

                        for (int j=0; j<numberIndividualsPerSpeciesInit; j++){

                            //New method triangular: sum of uniform distributions
                            int dateTransitory(speciesStartDate +  (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit/2*(disU01(gen)+disU01(gen)));

                            startFruitingDateTreeInit.push_back(dateTransitory);
                            endFruitingDateTreeInit.push_back(dateTransitory + fruitingLengthInit);

                            //Food quantity
                            foodQuantityTreeInit.push_back(0);//No food at first (will be updated afterwards)

                            //Coordinates

                            /////////// ADDED ON THE 01/02/2021
                            if(isHeterogeneous=="YES"){
                                counterTree+=1;
                                if(counterTree>=floor(numberIndividualsPerSpeciesInit/(double) numberPatch)){
                                    if(refPatch!=patchesCoordinatesX.size()-1){
                                        refPatch+=1;
                                        counterTree=0;
                                        disNormalPatchX.param(std::normal_distribution<double>::param_type(patchesCoordinatesX[refPatch],mapSizeInit/dispersion));//
                                        disNormalPatchY.param(std::normal_distribution<double>::param_type(patchesCoordinatesY[refPatch],mapSizeInit/dispersion));//
                                    }
                                }
                                double locX(disNormalPatchX(gen));
                                while(locX < 0 || locX > mapSizeInit){
                                    locX=disNormalPatchX(gen);
                                }
                                double locY(disNormalPatchY(gen));
                                while(locY < 0 || locY > mapSizeInit){
                                    locY=disNormalPatchY(gen);
                                }
                                xCoordinateTreeInit.push_back(locX);
                                yCoordinateTreeInit.push_back(locY);
                            }
                            else{
                            ////////////////////////////////////
                                xCoordinateTreeInit.push_back(disULocTree(gen));
                                yCoordinateTreeInit.push_back(disULocTree(gen));
                            }

                            //Update the sum for mean and var calculation
                            dateStartVector.push_back(dateTransitory);
                            dateEndVector.push_back(dateTransitory+fruitingLengthInit);

                            speciesTreeInit.push_back(i);
                        }

                        /*Calculating the operational starting date*/
                        //std::sort (dateStartVector.begin(), dateStartVector.end());//Ranging in increasing order starting dates
                        startFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate + 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));//Correction adding bc subsampling triggers that we never reach this extreme value. Hence Simon conducted simulation to see what correction we have to apply when sampling in triangular distrib, the extreme values we obtained on average.//startFruitingDateSpeciesOperationalInit.push_back(dateStartVector[0]);

                        /*Calculating the operational ending date*/
                        //std::sort (dateEndVector.rbegin(), dateEndVector.rend());//Ranging in decreasing order ending dates
                        endFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate  + (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit - 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));//endFruitingDateSpeciesOperationalInit.push_back(dateEndVector[0]);

                        if(fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit) < earliestStartFruitingDateInit
                           //&& i > 0
                           ){//For the associative scenario it is necessary to get rid of the period with no fruits or no fruit predictable
                            //NOT ANYMORE BECAUSE SIMULATING MANY CYCLES: therefore eliminating when only sp1
                                earliestStartFruitingDateInit=fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                        }
                        if(fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit) > latestEndFruitingDateInit){
                                latestEndFruitingDateInit=fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                        }

                        /*Extract the finally obtained mean and var for all species*/
                        //Calculate sum, mean, square sum and stdv
                        double sum = std::accumulate(dateStartVector.begin(), dateStartVector.end(), 0.0);//arguments are: first value, last value, initialisation value
                        double mean = sum / dateStartVector.size();
                        double sq_sum = std::inner_product(dateStartVector.begin(), dateStartVector.end(), dateStartVector.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
                        double stdev = std::sqrt(sq_sum / dateStartVector.size() - mean * mean);

                        startFruitingDateSpeciesAverageInitResult.push_back(mean);
                        startFruitingDateSpeciesStdevInitResult.push_back(stdev);
                    }

                    earliestStartFruitingDateInit=max(0.0, earliestStartFruitingDateInit);
                    latestEndFruitingDateInit=min((double) cycleLengthInit, latestEndFruitingDateInit);

                /*Determiner the associative strength between start-start events and start-end events*/
                //and /*Create the matrix for associative model: first col is the species, second is the distance to the date, and third col characterises which date is to use 0=start, 1=end*/

                //Create the result matrix
                int distanceStartStartInit[numberSpeciesInit][numberSpeciesInit];
                int distanceStartEndInit[numberSpeciesInit][numberSpeciesInit];
                int distanceEndStartInit[numberSpeciesInit][numberSpeciesInit];
                int distanceEndEndInit[numberSpeciesInit][numberSpeciesInit];
                int matrixInferenceAssociationInit[numberSpeciesInit][3];
                std::fill(*matrixInferenceAssociationInit, *matrixInferenceAssociationInit + 3*numberSpeciesInit, cycleLengthInit); //Initialise matrix with only cycleLengthInit

                //Create the min vectors of distance
                vector<double> vectorMinDistanceStart;
                vector<double> vectorMinDistanceEnd;

                //Loop to fill the result matrices and calculate alongside the minimum: first rows are always start-start then start-end
                for (int j=0; j<numberSpeciesInit; j++){

                    for (int k=0; k<numberSpeciesInit; k++){
                    if(j!=k){
                        //Updating the results matrices
                            distanceStartStartInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));
                            distanceStartEndInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]-fruitingLengthInit),cycleLengthInit));
                            distanceEndStartInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]+fruitingLengthInit),cycleLengthInit));
                            distanceEndEndInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));

                            vectorMinDistanceStart.push_back(min(distanceStartStartInit[j][k],
                                                                 distanceStartEndInit[j][k]
                                                                 )
                                                             );
                            vectorMinDistanceEnd.push_back(min(distanceEndStartInit[j][k],
                                                               distanceEndEndInit[j][k]
                                                               )
                                                           );
                        if(distanceStartStartInit[j][k]==min(distanceStartStartInit[j][k],distanceStartEndInit[j][k])){
                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                               (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] >=
                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                            ){
                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                       (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] ==
                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                                    ){
                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                            ){
                                               matrixInferenceAssociationInit[j][0]=k;
                                               matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                               matrixInferenceAssociationInit[j][2]=0;
                                            }
                                            else{
                                                //Do nothing
                                            }
                                    }
                                    else{
                                         matrixInferenceAssociationInit[j][0]=k;
                                         matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                         matrixInferenceAssociationInit[j][2]=0;
                                    }
                            }
                            else{
                                //Do nothing
                            }
                        }
                        else{
                           if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                               (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] >=
                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                            ){
                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                       (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] ==
                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                                    ){
                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                            ){
                                               matrixInferenceAssociationInit[j][0]=k;
                                               matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                               matrixInferenceAssociationInit[j][2]=1;
                                            }
                                            else{
                                                //Do nothing
                                            }
                                    }
                                    else{
                                         matrixInferenceAssociationInit[j][0]=k;
                                         matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                         matrixInferenceAssociationInit[j][2]=1;
                                    }
                            }
                        }
                    }
                }

                ////Checking what are the chosen referential species
                //cout << matrixInferenceAssociationInit[j][0] << " " << intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]] << endl;
            }

            //Calculate sum, mean, square sum and stdv
            double sum = std::accumulate(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), 0.0);//arguments are: first value, last value, initialisation value
            double mean = sum / vectorMinDistanceStart.size();
            double sq_sum = std::inner_product(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), vectorMinDistanceStart.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
            double stdev = std::sqrt(sq_sum / vectorMinDistanceStart.size() - mean * mean);

            double averageMinDistanceStartInit(mean);
            double sdMinDistanceStartInit(stdev);

            sum = std::accumulate(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), 0.0);//arguments are: first value, last value, initialisation value
            mean = sum / vectorMinDistanceEnd.size();
            sq_sum = std::inner_product(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), vectorMinDistanceEnd.begin(), 0.0);//arguments are: first value to End with, last value, first value to match to the other one, initialisation value
            stdev = std::sqrt(sq_sum / vectorMinDistanceEnd.size() - mean * mean);

            double averageMinDistanceEndInit(mean);
            double sdMinDistanceEndInit(stdev);

            //--------------
            //Running models
            //--------------

            /*Run model: NULL*/

            cout << "I work until line setting loop for model: NULL" << endl;

            nullModel(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double workingMemoryTime=
                    workingMemoryTimeInit, //max(2.0, (double)fruitingLengthInit/20),
                    //double timeLapseNotToReturn=
                    timeLapseNotToReturnInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Null_model_telescopic",
                    // bool extendibleArms,
                    true,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,

                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,

                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0
            );

            /*Run model: NULL WORKING MEMORY*/
            cout << "I work until line setting loop for model: NULL WORKING MEMORY" << endl;

            nullModelWorkingMemory(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double workingMemoryTime=
                    workingMemoryTimeInit,
                    //double timeLapseNotToReturn=
                    timeLapseNotToReturnInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Null_working_memory_model_telescopic",
                    // bool extendibleArms,
                    true,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,

                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,

                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0

            );

            /*Run model: CHRONOLOGICAL*/
            cout << "I work until line setting loop for model: CHRONOLOGICAL" << endl;

            chronologicalMemoryModel(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double workingMemoryTime=
                    workingMemoryTimeInit,
                    //double timeLapseNotToReturn=
                    timeLapseNotToReturnInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Chronological_model_telescopic",
                    // bool extendibleArms,
                    true,
                    //double weightInFavourOfLongtermKnowledge,
                    weightInFavourOfLongtermKnowledgeInit,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,
                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,
                    //std::vector<int> startFruitingDateSpeciesOperational=
                    startFruitingDateSpeciesOperationalInit,
                    //std::vector<int> endFruitingDateSpeciesOperational=
                    endFruitingDateSpeciesOperationalInit,
                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0

            );

            /*Run model: ASSOCIATIVE*/
            cout << "I work until line setting loop for model: ASSOCIATIVE" << endl;

            associativeMemoryModel(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double workingMemoryTime=
                    workingMemoryTimeInit,
                    //double timeLapseNotToReturn=
                    timeLapseNotToReturnInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Associative_model_telescopic",
                    // bool extendibleArms,
                    true,
                    //double weightInFavourOfLongtermKnowledge,
                    weightInFavourOfLongtermKnowledgeInit,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,

                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,
                    //int matrixInferenceAssociation,
                    *matrixInferenceAssociationInit,

                    //std::vector<int> startFruitingDateSpeciesOperational=
                    startFruitingDateSpeciesOperationalInit,
                    //std::vector<int> endFruitingDateSpeciesOperational=
                    endFruitingDateSpeciesOperationalInit,

                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0

            );


            /*Run model: OMNISCIENT*/
            cout << "I work until line setting loop for model: OMNISCIENT" << endl;

            omniscientModel(
                    //Run ID
                    //int step=
                    s,

                    //Output variable
                    //string pathOfTheFileOutput=
                    pathOfTheFileOutputInitMdf,

                    //Time
                    //double timer=
                    timerInit,
                    //double previousTimer=
                    previousTimerInit,


                    //double cycleLength=
                    cycleLengthInit,

                    //Agent ability
                    //double visionRadius=
                    visionRadiusInit,
                    //double coefficientDistanceToTime=
                    coefficientDistanceToTimeInit,
                    //string memoryType=
                    "Omniscient_model_telescopic",
                    // bool extendibleArms,
                    true,

                    //Agent location and success
                    //double xCoordinateAgent=
                    xCoordinateAgentInit,
                    //double yCoordinateAgent=
                    yCoordinateAgentInit,

                    //Initial environmental condition

                    //On Species
                    //int numberSpecies=
                    numberSpeciesInit,
                    //int numberIndividualsPerSpecies=
                    numberIndividualsPerSpeciesInit,
                    //vector<int>quantityFoodMaxTreeForSpecies,
                    quantityFoodMaxTreeForSpeciesInit,
                    //int numberTreesRareSpeciesInit,
                    numberTreesRareSpeciesInit,
                    //double fruitingLength=
                    fruitingLengthInit,
                    //int earliestStartFruitingDate=
                    earliestStartFruitingDateInit,
                    //int latestEndFruitingDate=
                    latestEndFruitingDateInit,
                    //double averageMinDistanceStart=
                    averageMinDistanceStartInit,
                    //double sdMinDistanceStart,
                    sdMinDistanceStartInit,
                    //double averageMinDistanceEnd=
                    averageMinDistanceEndInit,
                    //double sdMinDistanceEnd,
                    sdMinDistanceEndInit,

                    //std::vector<double> intraSpeciesSynchrony=
                    intraSpeciesSynchronyInit,
                    //double percentSpeciesWithLowSynchrony=
                    proportionLowSynchrony,

                    //On Tree
                    //std::vector<double> xCoordinateTree=
                    xCoordinateTreeInit,
                    //std::vector<double> yCoordinateTree=
                    yCoordinateTreeInit,
                    //std::vector<double> foodQuantityTree=
                    foodQuantityTreeInit,
                    //std::vector<int> speciesTree=
                    speciesTreeInit,
                    //std::vector<int> endFruitingDateTree=
                    endFruitingDateTreeInit,
                    //std::vector<int> startFruitingDateTree=
                    startFruitingDateTreeInit,
                    //double probaTreeAlreadyDepleted,
                    0.0
            );

        }
    }


//
//
//// Set back spatial environment to homogeneous
//isHeterogeneous="NO";
//
///////////////////////////////////////
//// SCENARIO 3: Working memory span and environment composition variation
///////////////////////////////////////
//
//    std::vector<int> workingMemoryLengthVector = {1,5,10,15,20,25,30,40,50,60};
//
//    vector<double> valueForSynchronyInWML={0.9};
//    for (int iSS=0; iSS < valueForSynchronyInWML.size(); iSS++){
//    //Fixing synchrony
//    double valueToUseSynchrony(valueForSynchronyInWML[iSS]);//Synchrony
//
//    //Testing for different repartition Nspe and Nind
//
//    std::vector<int> vectorNumberSpecies = {10,15,30,50,100};
//    std::vector<int> vectorNumberTrees = {300,200,100,60,30};
//
//    for(int r=0; r<vectorNumberSpecies.size(); r++){
//    vector<int> quantityFoodMaxTreeForSpeciesInitWM(vectorNumberSpecies[r]*vectorNumberTrees[r], 1);
//    unsigned int s;
//    #pragma omp parallel for
//    for (unsigned s=0; s<300; s++){
//        std::string pathOfTheFileOutputInitMdf(pathOfTheFileOutputInit);
//        pathOfTheFileOutputInitMdf.append(to_string(s));
//        pathOfTheFileOutputInitMdf.append("_");
//        pathOfTheFileOutputInitMdf.append(to_string(valueToUseSynchrony));
//        pathOfTheFileOutputInitMdf.append("_");
//        pathOfTheFileOutputInitMdf.append(to_string(vectorNumberSpecies[r]));
//        pathOfTheFileOutputInitMdf.append("_WM.txt");
//        std::ofstream outputFlux;
//        outputFlux.open(pathOfTheFileOutputInitMdf.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
//        if(outputFlux)  // Write the file only if correctly opened
//        {
//        }
//        else //If impossible to write in the file, say it
//        {
//            cout << "ERROR: Impossible to open the output file" << endl;
//        }
//
//        outputFlux << "Run_simulation" << " " << "Memory_type" << " " << "Time_length_run" << " " << "Days_wo_food" << " " << "Foraging_success_cum" << " " << "Foraging_success_cum_target_only" << " " << "Foraging_success_tot" << " " << "Foraging_success_tot_target_only" << " " << "Percent_tree_visited" << " " <<
//        "Percent_targets_with_food" << " " << "Number_trees_rare_species" << " " << "Rate_trees_rare_eaten" << " " << "Proba_competition" << " " << "Working_memory_length_time_unit" << " " << "Time_not_to_return" << " " << "Vision_radius" << " " << "Length_cycle" << " " << "Length_fruiting" << " " << "Number_species" << " " << "Number_individuals_per_species" << " " <<
//        "Intra_species_synchrony" << " " << "Percentage_low_synchrony" << " " << " " << "Distance_min_start" << " " << "Sd_min_start" << " " << "Distance_min_end" << " " << "Sd_min_end" << " " << "Weight_for_longterm_memory" << endl;
//
//        outputFlux.close();
//
//            //--------------
//            //INITIALISATION
//            //--------------
//
//            /*Time-associated variables*/
//            double timerInit(-cycleLengthInit);
//            double previousTimerInit(-cycleLengthInit/4);
//            double coefficientDistanceToTimeInit((double) 1/((double) 1/cycleLengthInit*(double) mapSizeInit*(double) mapSizeInit/(double) visionRadiusInit/2.0)); // This is set up by making the hypothesis that during the working memory 1/5 of the trees can be visited. We assume that, during the working memory, the individual hence travels a distance d
//            // which makes it cover d*visionRadius surface (approximation not considering overlap, hence overestimation). Applying the "produit en croix", we find that d can't be more than 1/5*mapSize**2/visionRadius. Hence this speed if defining working memory as FruitingLength/2
//
//            /*Species-level variables*/
//
//            vector<int> startFruitingDateSpeciesOperationalInit;
//            vector<int> endFruitingDateSpeciesOperationalInit;
//
//            vector<int> startFruitingDateSpeciesAverageInitResult;
//            vector<int> startFruitingDateSpeciesStdevInitResult;
//
//            vector<double> intraSpeciesSynchronyInit;
//
//            /*Tree-level variables*/
//            vector<double> xCoordinateTreeInit;
//            vector<double> yCoordinateTreeInit;
//            vector<int> startFruitingDateTreeInit;
//            vector<int> endFruitingDateTreeInit;
//            vector<int> speciesTreeInit;
//            vector<double> foodQuantityTreeInit;
//
//            /*Agent-level variables*/
//            double xCoordinateAgentInit(disULocTree(gen));
//            double yCoordinateAgentInit(disULocTree(gen));
//            double workingMemoryTimeInit(fruitingLengthInit/2.0);
//            double weightInFavourOfLongtermKnowledgeInit(0.0); //Only for associative/chronological memory, should be comprised between [-1;1]. Used for averaging between observed (=estimated) fruiting, and Longterm knowledge on fruiting dates
//
//            /* Environment */
//            double earliestStartFruitingDateInit(cycleLengthInit);
//            double latestEndFruitingDateInit(0.0);
//
//            std::uniform_real_distribution<double> disU01(0, 1); //Uniform distrib
//            std::bernoulli_distribution bernouilliDisInit(0.5);
//
//           /*Setting the environment: resource distribution and tree features*/
//            for (int i=0; i<numberSpeciesInit; i++){
//
//                    intraSpeciesSynchronyInit.push_back(intraSpeciesSynchronyValuesInit[iSS]);
//                    //int speciesToSplit(0);
//                    vector<int> dateStartVector;
//                    vector<int> dateEndVector;
//                        int speciesStartDate(disUDateSpecies(gen));
//
//                        //ADDED ON THE 01/02/2021
//                        //If heterogeneous distribution
//                        std::vector<int> patchesCoordinatesX(numberPatch,0);
//                        std::vector<int> patchesCoordinatesY(numberPatch,0);
//                        for(int p=0; p<patchesCoordinatesX.size(); p++){
//                            patchesCoordinatesX[p]=disULocPatch(gen);
//                            patchesCoordinatesY[p]=disULocPatch(gen);
//                        }
//                        int counterTree=0;
//                        int refPatch=0;
//                        std::normal_distribution<double> disNormalPatchX(patchesCoordinatesX[0],mapSizeInit/dispersion);//
//                        std::normal_distribution<double> disNormalPatchY(patchesCoordinatesY[0],mapSizeInit/dispersion);//
//                        //////////////////////////
//
//                        for (int j=0; j<numberIndividualsPerSpeciesInit; j++){
//                            //New method triangular: sum of uniform distributions (NOTE: in the article for formula, we are using the Tm (the max point), and not the TS, therefore you can see that it sums the uniform - 1 !!!!!!!!
//                            int dateTransitory(speciesStartDate +  (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit/2*(disU01(gen)+disU01(gen)));
//
//                            startFruitingDateTreeInit.push_back(dateTransitory);
//                            endFruitingDateTreeInit.push_back(dateTransitory + fruitingLengthInit);
//
//                            //Food quantity
//                            foodQuantityTreeInit.push_back(0);//No food at first (will be updated afterwards)
//
//                            //Coordinates
//
//                            /////////// ADDED ON THE 01/02/2021
//                            if(isHeterogeneous=="YES"){
//                                counterTree+=1;
//                                if(counterTree>=floor(numberIndividualsPerSpeciesInit/(double) numberPatch)){
//                                    if(refPatch!=patchesCoordinatesX.size()-1){
//                                        refPatch+=1;
//                                        counterTree=0;
//                                        disNormalPatchX.param(std::normal_distribution<double>::param_type(patchesCoordinatesX[refPatch],mapSizeInit/dispersion));//
//                                        disNormalPatchY.param(std::normal_distribution<double>::param_type(patchesCoordinatesY[refPatch],mapSizeInit/dispersion));//
//                                    }
//                                }
//                                double locX(disNormalPatchX(gen));
//                                while(locX < 0 || locX > mapSizeInit){
//                                    locX=disNormalPatchX(gen);
//                                }
//                                double locY(disNormalPatchY(gen));
//                                while(locY < 0 || locY > mapSizeInit){
//                                    locY=disNormalPatchY(gen);
//                                }
//                                xCoordinateTreeInit.push_back(locX);
//                                yCoordinateTreeInit.push_back(locY);
//                            }
//                            else{
//                            ////////////////////////////////////
//                                xCoordinateTreeInit.push_back(disULocTree(gen));
//                                yCoordinateTreeInit.push_back(disULocTree(gen));
//                            }
//
//                            //Update vector for mean and var calculation
//                            dateStartVector.push_back(dateTransitory);
//                            dateEndVector.push_back(dateTransitory+fruitingLengthInit);
//
//                            speciesTreeInit.push_back(i);
//                        }
//
//                    /*Calculating the operational starting date*/
//                    //std::sort (dateStartVector.begin(), dateStartVector.end());//Ranging in increasing order starting dates
//                    startFruitingDateSpeciesOperationalInit.push_back(speciesStartDate + 0.625/sqrt(vectorNumberTrees[r])*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit);//Correction adding bc subsampling triggers that we never reach this extreme value. Hence Simon conducted simulation to see what correction we have to apply when sampling in triangular distrib, the extreme values we obtained on average.//startFruitingDateSpeciesOperationalInit.push_back(dateStartVector[0]);
//
//                    /*Calculating the operational ending date*/
//                    //std::sort (dateEndVector.rbegin(), dateEndVector.rend());//Ranging in decreasing order ending dates
//                    endFruitingDateSpeciesOperationalInit.push_back(speciesStartDate  + (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit - 0.625/sqrt(vectorNumberTrees[r])*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit);//endFruitingDateSpeciesOperationalInit.push_back(dateEndVector[0]);
//
//                    if(fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit) < earliestStartFruitingDateInit
//                       ){
//                            earliestStartFruitingDateInit=fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
//                    }
//                    if(fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit) > latestEndFruitingDateInit
//                       ){
//                            latestEndFruitingDateInit=fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
//                    }
//
//                    /*Extract the finally obtained mean and var for all species*/
//                    //Calculate sum, mean, square sum and stdv
//                    double sum = std::accumulate(dateStartVector.begin(), dateStartVector.end(), 0.0);//arguments are: first value, last value, initialisation value
//                    double mean = sum / dateStartVector.size();
//                    double sq_sum = std::inner_product(dateStartVector.begin(), dateStartVector.end(), dateStartVector.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
//                    double stdev = std::sqrt(sq_sum / dateStartVector.size() - mean * mean);
//
//                    if(intraSpeciesSynchronyInit[i] > 0){
//                    }
//                    else{
//                        mean=cycleLengthInit/2.0;
//                        stdev=sqrt(cycleLengthInit*cycleLengthInit/12.0);
//                    }
//                    startFruitingDateSpeciesAverageInitResult.push_back(mean);
//                    startFruitingDateSpeciesStdevInitResult.push_back(stdev);
//                }
//
//                earliestStartFruitingDateInit=max(0.0, earliestStartFruitingDateInit);
//                latestEndFruitingDateInit=min((double) cycleLengthInit, latestEndFruitingDateInit);
//
//                /*Determiner the associative strength between start-start events and start-end events*/
//                //and /*Create the matrix for associative model: first col is the species, second is the distance to the date, and third col characterises which date is to use 0=start, 1=end*/
//
//                //Create the result matrix
//                int distanceStartStartInit[vectorNumberSpecies[r]][vectorNumberSpecies[r]];
//                int distanceStartEndInit[vectorNumberSpecies[r]][vectorNumberSpecies[r]];
//                int distanceEndStartInit[vectorNumberSpecies[r]][vectorNumberSpecies[r]];
//                int distanceEndEndInit[vectorNumberSpecies[r]][vectorNumberSpecies[r]];
//                int matrixInferenceAssociationInit[vectorNumberSpecies[r]][3];
//                std::fill(*matrixInferenceAssociationInit, *matrixInferenceAssociationInit + 3*vectorNumberSpecies[r], cycleLengthInit); //Initialise matrix with only cycleLengthInit
//
//                //Create the min vectors of distance
//                vector<double> vectorMinDistanceStart;
//                vector<double> vectorMinDistanceEnd;
//
//                //Loop to fill the result matrices and calculate alongside the minimum: first rows are always start-start then start-end
//                for (int j=0; j<vectorNumberSpecies[r]; j++){
//
//                    for (int k=0; k<vectorNumberSpecies[r]; k++){
//                    if(j!=k){
//                        //Updating the results matrices
//                            distanceStartStartInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));
//                            distanceStartEndInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]-fruitingLengthInit),cycleLengthInit));
//                            distanceEndStartInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]+fruitingLengthInit),cycleLengthInit));
//                            distanceEndEndInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));
//
//                            vectorMinDistanceStart.push_back(min(distanceStartStartInit[j][k],
//                                                                 distanceStartEndInit[j][k]
//                                                                 )
//                                                             );
//                            vectorMinDistanceEnd.push_back(min(distanceEndStartInit[j][k],
//                                                               distanceEndEndInit[j][k]
//                                                               )
//                                                           );
//                        if(distanceStartStartInit[j][k]==min(distanceStartStartInit[j][k],distanceStartEndInit[j][k])){
//                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
//                               (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] >=
//                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
//                            ){
//                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
//                                       (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] ==
//                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
//
//                                    ){
//                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
//                                            ){
//                                               matrixInferenceAssociationInit[j][0]=k;
//                                               matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
//                                               matrixInferenceAssociationInit[j][2]=0;
//                                            }
//                                            else{
//                                                //Do nothing
//                                            }
//                                    }
//                                    else{
//                                         matrixInferenceAssociationInit[j][0]=k;
//                                         matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
//                                         matrixInferenceAssociationInit[j][2]=0;
//                                    }
//                            }
//                            else{
//                                //Do nothing
//                            }
//                        }
//                        else{
//                           if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
//                               (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] >=
//                               (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
//                            ){
//                                    if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
//                                       (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] ==
//                                       (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
//
//                                    ){
//                                            if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
//                                            ){
//                                               matrixInferenceAssociationInit[j][0]=k;
//                                               matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
//                                               matrixInferenceAssociationInit[j][2]=1;
//                                            }
//                                            else{
//                                                //Do nothing
//                                            }
//                                    }
//                                    else{
//                                         matrixInferenceAssociationInit[j][0]=k;
//                                         matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
//                                         matrixInferenceAssociationInit[j][2]=1;
//                                    }
//                            }
//                        }
//                    }
//                }
//            }
//
//            //Calculate sum, mean, square sum and stdv
//            double sum = std::accumulate(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), 0.0);//arguments are: first value, last value, initialisation value
//            double mean = sum / vectorMinDistanceStart.size();
//            double sq_sum = std::inner_product(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), vectorMinDistanceStart.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
//            double stdev = std::sqrt(sq_sum / vectorMinDistanceStart.size() - mean * mean);
//
//            double averageMinDistanceStartInit(mean);
//            double sdMinDistanceStartInit(stdev);
//
//            sum = std::accumulate(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), 0.0);//arguments are: first value, last value, initialisation value
//            mean = sum / vectorMinDistanceEnd.size();
//            sq_sum = std::inner_product(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), vectorMinDistanceEnd.begin(), 0.0);//arguments are: first value to End with, last value, first value to match to the other one, initialisation value
//            stdev = std::sqrt(sq_sum / vectorMinDistanceEnd.size() - mean * mean);
//
//            double averageMinDistanceEndInit(mean);
//            double sdMinDistanceEndInit(stdev);
//
//        for(int wml=0; wml<workingMemoryLengthVector.size(); wml++){
//
//        //--------------
//        //Running models
//        //--------------
//
//        /*Run model: NULL WORKING MEMORY*/
//        cout << "I work until line setting loop for model: NULL WORKING MEMORY" << endl;
//
//        nullModelWorkingMemory(
//                //Run ID
//                //int step=
//                s,
//
//                //Output variable
//                //string pathOfTheFileOutput=
//                pathOfTheFileOutputInitMdf,
//
//                //Time
//                //double timer=
//                timerInit,
//                //double previousTimer=
//                previousTimerInit,
//                //double cycleLength=
//                cycleLengthInit,
//
//                //Agent ability
//                //double visionRadius=
//                visionRadiusInit,
//                //double workingMemoryTime=
//                workingMemoryLengthVector[wml],
//                //double timeLapseNotToReturn=
//                timeLapseNotToReturnInit,
//                //double coefficientDistanceToTime=
//                coefficientDistanceToTimeInit,
//                //string memoryType=
//                "Null_working_memory_model_telescopic",
//
//                // bool extendibleArms,
//                true,
//
//                //Agent location and success
//                //double xCoordinateAgent=
//                xCoordinateAgentInit,
//                //double yCoordinateAgent=
//                yCoordinateAgentInit,
//
//                //Initial environmental condition
//
//                //On Species
//                //int numberSpecies=
//                vectorNumberSpecies[r],
//                //int numberIndividualsPerSpecies=
//                vectorNumberTrees[r],
//				//vector<int>quantityFoodMaxTreeForSpecies,
//				quantityFoodMaxTreeForSpeciesInitWM,
//				//int numberTreesRareSpeciesInit,
//				numberTreesRareSpeciesInit,
//                //double fruitingLength=
//                fruitingLengthInit,
//                //int earliestStartFruitingDate=
//                earliestStartFruitingDateInit,
//                //int latestEndFruitingDate=
//                latestEndFruitingDateInit,
//                //double averageMinDistanceStart=
//                averageMinDistanceStartInit,
//                //double sdMinDistanceStart,
//                sdMinDistanceStartInit,
//                //double averageMinDistanceEnd=
//                averageMinDistanceEndInit,
//                //double sdMinDistanceEnd,
//                sdMinDistanceEndInit,
//
//                //std::vector<double> intraSpeciesSynchrony=
//                intraSpeciesSynchronyInit,
//                //double percentSpeciesWithLowSynchrony=
//                1.0,
//
//                //On Tree
//                //std::vector<double> xCoordinateTree=
//                xCoordinateTreeInit,
//                //std::vector<double> yCoordinateTree=
//                yCoordinateTreeInit,
//                //std::vector<double> foodQuantityTree=
//                foodQuantityTreeInit,
//                //std::vector<int> speciesTree=
//                speciesTreeInit,
//                //std::vector<int> endFruitingDateTree=
//                endFruitingDateTreeInit,
//                //std::vector<int> startFruitingDateTree=
//                startFruitingDateTreeInit,
//                //double probaTreeAlreadyDepleted,
//                0.0
//            );
//         }
//      }
//   }
//}

    ////////////////////////////////////////////////
    // SCENARIO 4: "THE DELAY EXPERIMENT" Delaying the fruiting time from known timing in chronological memory
    ////////////////////////////////////////////////

    // NOTE: Contrary to what was done previously, here, simulations for different parameterization won't have the same environment (it was convenient to avoid time-consumming useless modifications
    vector<int> vectorTranslation={-60,-45,-30,-15,-5,-1,1,5,15,30,45,60};

    for(int x=0; x < vectorTranslation.size(); x++){
        int translation(vectorTranslation[x]);

        //unsigned int s;
        #pragma omp parallel for
        for (unsigned int s=0; s<1; s++){

            std::string pathOfTheFileOutputInitMdf(pathOfTheFileOutputInit);
            pathOfTheFileOutputInitMdf.append(to_string(s));
            pathOfTheFileOutputInitMdf.append("_");
            pathOfTheFileOutputInitMdf.append(to_string(translation));
            pathOfTheFileOutputInitMdf.append("_IS_delay.txt");
            std::ofstream outputFlux;
            outputFlux.open(pathOfTheFileOutputInitMdf.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
            if(outputFlux)  // Write the file only if correctly opened
            {
            }
            else //If impossible to write in the file, say it
            {
                cout << "ERROR: Impossible to open the output file" << endl;
            }

            outputFlux << "Run_simulation" << " " << "Memory_type" << " " << "Time_length_run" << " " << "Days_wo_food" << " " << "Foraging_success_cum" << " " << "Foraging_success_cum_target_only" << " " << "Foraging_success_tot" << " " << "Foraging_success_tot_target_only" << " " << "Percent_tree_visited" << " " <<
            "Percent_targets_with_food" << " " << "Number_trees_rare_species" << " " << "Rate_trees_rare_eaten" << " " << "Proba_competition" << " " << "Working_memory_length_time_unit" << " " << "Time_not_to_return" << " " << "Vision_radius" << " " << "Length_cycle" << " " << "Length_fruiting" << " " << "Number_species" << " " << "Number_individuals_per_species" << " " <<
            "Intra_species_synchrony" << " " << "Percentage_low_synchrony" << " " << "Distance_min_start" << " " << "Sd_min_start" << " " << "Distance_min_end" << " " << "Sd_min_end" << " " << "Weight_for_longterm_memory" << endl;

            outputFlux.close();

            for(int iSS=0; iSS<intraSpeciesSynchronyValuesInit.size(); iSS++){

                //--------------
                //INITIALISATION
                //--------------

                /*Time-associated variables*/
                double timerInit(-cycleLengthInit);//NO NEED EARLIER BC ONLY CHRONOLOGICAL
                double previousTimerInit(-cycleLengthInit);
                double coefficientDistanceToTimeInit((double) 1/((double) 1/cycleLengthInit*(double) mapSizeInit*(double) mapSizeInit/(double) visionRadiusInit/2.0));

                // This is set up by making the hypothesis that during the cycle all of the trees can be visited. We assume that, during the cycle the individual hence travels a distance d
                // which makes it cover d*visionRadius surface (approximation not considering overlap, hence overestimation). Applying the "produit en croix", we find that d can't be more than 1/2*mapSize**2/visionRadius. Hence this speed if defining working memory as FruitingLength/2
                /*Species-level variables*/

                vector<int> startFruitingDateSpeciesOperationalInit;
                vector<int> endFruitingDateSpeciesOperationalInit;

                vector<int> startFruitingDateSpeciesAverageInitResult;
                vector<int> startFruitingDateSpeciesStdevInitResult;

                vector<double> intraSpeciesSynchronyInit;

                /*Tree-level variables*/
                vector<double> xCoordinateTreeInit;
                vector<double> yCoordinateTreeInit;
                vector<int> startFruitingDateTreeInit;
                vector<int> endFruitingDateTreeInit;
                vector<int> speciesTreeInit;
                vector<double> foodQuantityTreeInit;

                /*Agent-level variables*/
                double xCoordinateAgentInit(disULocTree(gen));
                double yCoordinateAgentInit(disULocTree(gen));
                double workingMemoryTimeInit(15.0);
                double weightInFavourOfLongtermKnowledgeInit(1.0); //Only for associative/chronological memory, should be comprised between [-1;1]. Used for averaging between observed (=estimated) fruiting, and Longterm knowledge on fruiting dates

                /* Environment */
                double earliestStartFruitingDateInit(cycleLengthInit);
                double latestEndFruitingDateInit(0.0);

                std::uniform_real_distribution<double> disU01(0, 1); //Uniform distrib
                std::bernoulli_distribution bernouilliDisInit(0.5);

                /*Setting the environment: resource distribution and tree features*/
                for (int i=0; i<numberSpeciesInit; i++){

                    intraSpeciesSynchronyInit.push_back(intraSpeciesSynchronyValuesInit[iSS]);
                    //int speciesToSplit(0);
                    vector<int> dateStartVector;
                    vector<int> dateEndVector;
                        int speciesStartDate(disUDateSpecies(gen));

                        //ADDED ON THE 01/02/2021
                        //If heterogeneous distribution
                        std::vector<int> patchesCoordinatesX(numberPatch,0);
                        std::vector<int> patchesCoordinatesY(numberPatch,0);
                        for(int p=0; p<patchesCoordinatesX.size(); p++){
                            patchesCoordinatesX[p]=disULocPatch(gen);
                            patchesCoordinatesY[p]=disULocPatch(gen);
                        }
                        int counterTree=0;
                        int refPatch=0;
                        std::normal_distribution<double> disNormalPatchX(patchesCoordinatesX[0],mapSizeInit/dispersion);//
                        std::normal_distribution<double> disNormalPatchY(patchesCoordinatesY[0],mapSizeInit/dispersion);//
                        //////////////////////////

                        for (int j=0; j<numberIndividualsPerSpeciesInit; j++){
                            //New method triangular: sum of uniform distributions (NOTE: in the article for formula, we are using the Tm (the max point), and not the TS, therefore you can see that it sums the uniform - 1 !!!!!!!!
                            int dateTransitory(speciesStartDate +  (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit/2*(disU01(gen)+disU01(gen)));

                            startFruitingDateTreeInit.push_back(dateTransitory);
                            endFruitingDateTreeInit.push_back(dateTransitory + fruitingLengthInit);

                            //Food quantity
                            foodQuantityTreeInit.push_back(0);//No food at first (will be updated afterwards)

                            //Coordinates

                            /////////// ADDED ON THE 01/02/2021
                            if(isHeterogeneous=="YES"){
                                counterTree+=1;
                                if(counterTree>=floor(numberIndividualsPerSpeciesInit/(double) numberPatch)){
                                    if(refPatch!=patchesCoordinatesX.size()-1){
                                        refPatch+=1;
                                        counterTree=0;
                                        disNormalPatchX.param(std::normal_distribution<double>::param_type(patchesCoordinatesX[refPatch],mapSizeInit/dispersion));//
                                        disNormalPatchY.param(std::normal_distribution<double>::param_type(patchesCoordinatesY[refPatch],mapSizeInit/dispersion));//
                                    }
                                }
                                double locX(disNormalPatchX(gen));
                                while(locX < 0 || locX > mapSizeInit){
                                    locX=disNormalPatchX(gen);
                                }
                                double locY(disNormalPatchY(gen));
                                while(locY < 0 || locY > mapSizeInit){
                                    locY=disNormalPatchY(gen);
                                }
                                xCoordinateTreeInit.push_back(locX);
                                yCoordinateTreeInit.push_back(locY);
                            }
                            else{
                            ////////////////////////////////////
                                xCoordinateTreeInit.push_back(disULocTree(gen));
                                yCoordinateTreeInit.push_back(disULocTree(gen));
                            }

                            //Update vector for mean and var calculation
                            dateStartVector.push_back(dateTransitory);
                            dateEndVector.push_back(dateTransitory+fruitingLengthInit);

                            speciesTreeInit.push_back(i);
                        }

                        /*Calculating the operational starting date*/
                        //std::sort (dateStartVector.begin(), dateStartVector.end());//Ranging in increasing order starting dates
                        startFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate + 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));//Correction adding bc subsampling triggers that we never reach this extreme value. Hence Simon conducted simulation to see what correction we have to apply when sampling in triangular distrib, the extreme values we obtained on average.//startFruitingDateSpeciesOperationalInit.push_back(dateStartVector[0]);
                        /*Calculating the operational ending date*/
                        //std::sort (dateEndVector.rbegin(), dateEndVector.rend());//Ranging in decreasing order ending dates
                        endFruitingDateSpeciesOperationalInit.push_back(floor(speciesStartDate  + (1 - intraSpeciesSynchronyInit[i])*cycleLengthInit - 0.625/sqrt(numberIndividualsPerSpeciesInit)*(1 - intraSpeciesSynchronyInit[i])*cycleLengthInit));//endFruitingDateSpeciesOperationalInit.push_back(dateEndVector[0]);

                        if(fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit) < earliestStartFruitingDateInit
                           ){
                                earliestStartFruitingDateInit=fmod(10*cycleLengthInit + startFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                        }
                        if(fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit) > latestEndFruitingDateInit
                           ){
                                latestEndFruitingDateInit=fmod(10*cycleLengthInit + endFruitingDateSpeciesOperationalInit[i], cycleLengthInit);
                        }


                        /*Extract the finally obtained mean and var for all species*/
                        //Calculate sum, mean, square sum and stdv
                        double sum = std::accumulate(dateStartVector.begin(), dateStartVector.end(), 0.0);//arguments are: first value, last value, initialisation value
                        double mean = sum / dateStartVector.size();
                        double sq_sum = std::inner_product(dateStartVector.begin(), dateStartVector.end(), dateStartVector.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
                        double stdev = std::sqrt(sq_sum / dateStartVector.size() - mean * mean);

                        if(intraSpeciesSynchronyInit[i] > 0){
                        }
                        else{
                            mean=cycleLengthInit/2.0;
                            stdev=sqrt(cycleLengthInit*cycleLengthInit/12.0);
                        }
                        startFruitingDateSpeciesAverageInitResult.push_back(mean);
                        startFruitingDateSpeciesStdevInitResult.push_back(stdev);
                    }

                    earliestStartFruitingDateInit=max(0.0, earliestStartFruitingDateInit);
                    latestEndFruitingDateInit=min((double) cycleLengthInit, latestEndFruitingDateInit);

                    /*Determiner the associative strength between start-start events and start-end events*/
                    //and /*Create the matrix for associative model: first col is the species, second is the distance to the date, and third col characterises which date is to use 0=start, 1=end*/

                    //Create the result matrix
                    int distanceStartStartInit[numberSpeciesInit][numberSpeciesInit];
                    int distanceStartEndInit[numberSpeciesInit][numberSpeciesInit];
                    int distanceEndStartInit[numberSpeciesInit][numberSpeciesInit];
                    int distanceEndEndInit[numberSpeciesInit][numberSpeciesInit];
                    int matrixInferenceAssociationInit[numberSpeciesInit][3];
                    std::fill(*matrixInferenceAssociationInit, *matrixInferenceAssociationInit + 3*numberSpeciesInit, cycleLengthInit); //Initialise matrix with only cycleLengthInit

                    //Create the min vectors of distance
                    vector<double> vectorMinDistanceStart;
                    vector<double> vectorMinDistanceEnd;

                    //Loop to fill the result matrices and calculate alongside the minimum: first rows are always start-start then start-end
                    for (int j=0; j<numberSpeciesInit; j++){

                        for (int k=0; k<numberSpeciesInit; k++){
                        if(j!=k){
                            //Updating the results matrices
                                distanceStartStartInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));
                                distanceStartEndInit[j][k]= floor(fmod(10*cycleLengthInit+(startFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]-fruitingLengthInit),cycleLengthInit));
                                distanceEndStartInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-startFruitingDateSpeciesOperationalInit[k]+fruitingLengthInit),cycleLengthInit));
                                distanceEndEndInit[j][k]= floor(fmod(10*cycleLengthInit+(endFruitingDateSpeciesOperationalInit[j]-endFruitingDateSpeciesOperationalInit[k]),cycleLengthInit));

                                vectorMinDistanceStart.push_back(min(distanceStartStartInit[j][k],
                                                                     distanceStartEndInit[j][k]
                                                                     )
                                                                 );
                                vectorMinDistanceEnd.push_back(min(distanceEndStartInit[j][k],
                                                                   distanceEndEndInit[j][k]
                                                                   )
                                                               );
                            if(distanceStartStartInit[j][k]==min(distanceStartStartInit[j][k],distanceStartEndInit[j][k])){
                                if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                   (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] >=
                                   (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                                ){
                                        if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                           (double) intraSpeciesSynchronyInit[k]/distanceStartStartInit[j][k] ==
                                           (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]

                                        ){
                                                if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                                ){
                                                   matrixInferenceAssociationInit[j][0]=k;
                                                   matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                                   matrixInferenceAssociationInit[j][2]=0;
                                                }
                                                else{
                                                    //Do nothing
                                                }
                                        }
                                        else{
                                             matrixInferenceAssociationInit[j][0]=k;
                                             matrixInferenceAssociationInit[j][1]=distanceStartStartInit[j][k];
                                             matrixInferenceAssociationInit[j][2]=0;
                                        }
                                }
                                else{
                                    //Do nothing
                                }
                            }
                            else{
                               if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                   (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] >=
                                   (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]
                                ){
                                        if(matrixInferenceAssociationInit[j][0]==cycleLengthInit ||
                                           (double) intraSpeciesSynchronyInit[k]/distanceStartEndInit[j][k] ==
                                           (double) intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]/matrixInferenceAssociationInit[j][1]

                                        ){
                                                if(matrixInferenceAssociationInit[j][0]==cycleLengthInit || intraSpeciesSynchronyInit[k] > intraSpeciesSynchronyInit[matrixInferenceAssociationInit[j][0]]
                                                ){
                                                   matrixInferenceAssociationInit[j][0]=k;
                                                   matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                                   matrixInferenceAssociationInit[j][2]=1;
                                                }
                                                else{
                                                    //Do nothing
                                                }
                                        }
                                        else{
                                             matrixInferenceAssociationInit[j][0]=k;
                                             matrixInferenceAssociationInit[j][1]=distanceStartEndInit[j][k];
                                             matrixInferenceAssociationInit[j][2]=1;
                                        }
                                }
                            }
                        }
                    }
                }

                //Calculate sum, mean, square sum and stdv
                double sum = std::accumulate(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), 0.0);//arguments are: first value, last value, initialisation value
                double mean = sum / vectorMinDistanceStart.size();
                double sq_sum = std::inner_product(vectorMinDistanceStart.begin(), vectorMinDistanceStart.end(), vectorMinDistanceStart.begin(), 0.0);//arguments are: first value to start with, last value, first value to match to the other one, initialisation value
                double stdev = std::sqrt(sq_sum / vectorMinDistanceStart.size() - mean * mean);

                double averageMinDistanceStartInit(mean);
                double sdMinDistanceStartInit(stdev);

                sum = std::accumulate(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), 0.0);//arguments are: first value, last value, initialisation value
                mean = sum / vectorMinDistanceEnd.size();
                sq_sum = std::inner_product(vectorMinDistanceEnd.begin(), vectorMinDistanceEnd.end(), vectorMinDistanceEnd.begin(), 0.0);//arguments are: first value to End with, last value, first value to match to the other one, initialisation value
                stdev = std::sqrt(sq_sum / vectorMinDistanceEnd.size() - mean * mean);

                double averageMinDistanceEndInit(mean);
                double sdMinDistanceEndInit(stdev);

            //--------------
            //Running models
            //--------------

            //cout << startFruitingDateSpeciesOperationalInit[0] << " " << endFruitingDateSpeciesOperationalInit[0] << endl;
            /*Run model: CHRONOLOGICAL*/
            cout << "I work until line setting loop for model: CHRONOLOGICAL" << endl;


            /* Create the erroneous species dates (i.e positive/negative translation) */

            vector<int> startFruitingDateSpeciesOperationalInitTranslated=startFruitingDateSpeciesOperationalInit;
            vector<int> endFruitingDateSpeciesOperationalInitTranslated=endFruitingDateSpeciesOperationalInit;

                for(int a=0; a<startFruitingDateSpeciesOperationalInitTranslated.size(); a++){
                    startFruitingDateSpeciesOperationalInitTranslated[a]+=translation;
                    endFruitingDateSpeciesOperationalInitTranslated[a]+=translation;
                }
                chronologicalMemoryModel(
                        //Run ID
                        //int step=
                        s,

                        //Output variable
                        //string pathOfTheFileOutput=
                        pathOfTheFileOutputInitMdf,

                        //Time
                        //double timer=
                        timerInit,
                        //double previousTimer=
                        previousTimerInit,
                        //double cycleLength=
                        cycleLengthInit,

                        //Agent ability
                        //double visionRadius=
                        visionRadiusInit,
                        //double workingMemoryTime=
                        workingMemoryTimeInit,
                        //double timeLapseNotToReturn=
                        timeLapseNotToReturnInit,
                        //double coefficientDistanceToTime=
                        coefficientDistanceToTimeInit,
                        //string memoryType=
                        "Chronological_model_telescopic",
                        // bool extendibleArms,
                        true,
                        //double weightInFavourOfLongtermKnowledge,
                        weightInFavourOfLongtermKnowledgeInit,

                        //Agent location and success
                        //double xCoordinateAgent=
                        xCoordinateAgentInit,
                        //double yCoordinateAgent=
                        yCoordinateAgentInit,

                        //Initial environmental condition

                        //On Species
                        //int numberSpecies=
                        numberSpeciesInit,
                        //int numberIndividualsPerSpecies=
                        numberIndividualsPerSpeciesInit,
                        //vector<int>quantityFoodMaxTreeForSpecies,
                        quantityFoodMaxTreeForSpeciesInit,
                        //int numberTreesRareSpeciesInit,
                        numberTreesRareSpeciesInit,
                        //double fruitingLength=
                        fruitingLengthInit,
                        //int earliestStartFruitingDate=
                        earliestStartFruitingDateInit,
                        //int latestEndFruitingDate=
                        latestEndFruitingDateInit,
                        //double averageMinDistanceStart=
                        averageMinDistanceStartInit,
                        //double sdMinDistanceStart,
                        sdMinDistanceStartInit,
                        //double averageMinDistanceEnd=
                        averageMinDistanceEndInit,
                        //double sdMinDistanceEnd,
                        sdMinDistanceEndInit,
                        //std::vector<double> intraSpeciesSynchrony=
                        intraSpeciesSynchronyInit,
                        //double percentSpeciesWithLowSynchrony=
                        1.0,
                        //std::vector<int> startFruitingDateSpeciesOperational=
                        startFruitingDateSpeciesOperationalInitTranslated,
                        //std::vector<int> endFruitingDateSpeciesOperational=
                        endFruitingDateSpeciesOperationalInitTranslated,

                        //On Tree
                        //std::vector<double> xCoordinateTree=
                        xCoordinateTreeInit,
                        //std::vector<double> yCoordinateTree=
                        yCoordinateTreeInit,
                        //std::vector<double> foodQuantityTree=
                        foodQuantityTreeInit,
                        //std::vector<int> speciesTree=
                        speciesTreeInit,
                        //std::vector<int> endFruitingDateTree=
                        endFruitingDateTreeInit,
                        //std::vector<int> startFruitingDateTree=
                        startFruitingDateTreeInit,
                        //double probaTreeAlreadyDepleted,
                        0.0

                );

                chronologicalMemoryModelDelayRescued(
                        //Run ID
                        //int step=
                        s,

                        //Output variable
                        //string pathOfTheFileOutput=
                        pathOfTheFileOutputInitMdf,

                        //Time
                        //double timer=
                        timerInit,
                        //double previousTimer=
                        previousTimerInit,
                        //double cycleLength=
                        cycleLengthInit,
                        translation,

                        //Agent ability
                        //double visionRadius=
                        visionRadiusInit,
                        //double workingMemoryTime=
                        workingMemoryTimeInit,
                        //double timeLapseNotToReturn=
                        timeLapseNotToReturnInit,
                        //double coefficientDistanceToTime=
                        coefficientDistanceToTimeInit,
                        //string memoryType=
                        "Chronological_model_telescopic_rescued",
                        // bool extendibleArms,
                        true,
                        //double weightInFavourOfLongtermKnowledge,
                        weightInFavourOfLongtermKnowledgeInit,

                        //Agent location and success
                        //double xCoordinateAgent=
                        xCoordinateAgentInit,
                        //double yCoordinateAgent=
                        yCoordinateAgentInit,

                        //Initial environmental condition

                        //On Species
                        //int numberSpecies=
                        numberSpeciesInit,
                        //int numberIndividualsPerSpecies=
                        numberIndividualsPerSpeciesInit,
                        //vector<int>quantityFoodMaxTreeForSpecies,
                        quantityFoodMaxTreeForSpeciesInit,
                        //int numberTreesRareSpeciesInit,
                        numberTreesRareSpeciesInit,
                        //double fruitingLength=
                        fruitingLengthInit,
                        //int earliestStartFruitingDate=
                        earliestStartFruitingDateInit,
                        //int latestEndFruitingDate=
                        latestEndFruitingDateInit,
                        //double averageMinDistanceStart=
                        averageMinDistanceStartInit,
                        //double sdMinDistanceStart,
                        sdMinDistanceStartInit,
                        //double averageMinDistanceEnd=
                        averageMinDistanceEndInit,
                        //double sdMinDistanceEnd,
                        sdMinDistanceEndInit,
                        //std::vector<double> intraSpeciesSynchrony=
                        intraSpeciesSynchronyInit,
                        //double percentSpeciesWithLowSynchrony=
                        1.0,
                        //std::vector<int> startFruitingDateSpeciesOperational=
                        startFruitingDateSpeciesOperationalInitTranslated,
                        //std::vector<int> endFruitingDateSpeciesOperational=
                        endFruitingDateSpeciesOperationalInitTranslated,

                        //On Tree
                        //std::vector<double> xCoordinateTree=
                        xCoordinateTreeInit,
                        //std::vector<double> yCoordinateTree=
                        yCoordinateTreeInit,
                        //std::vector<double> foodQuantityTree=
                        foodQuantityTreeInit,
                        //std::vector<int> speciesTree=
                        speciesTreeInit,
                        //std::vector<int> endFruitingDateTree=
                        endFruitingDateTreeInit,
                        //std::vector<int> startFruitingDateTree=
                        startFruitingDateTreeInit,
                        //double probaTreeAlreadyDepleted,
                        0.0
                );
            }
        }
    }

cout << "Ended batch of simulations" << endl;

}

