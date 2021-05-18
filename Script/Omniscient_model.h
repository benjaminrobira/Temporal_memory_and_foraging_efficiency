#ifndef OMNISCIENT_MODEL_H_INCLUDED
#define OMNISCIENT_MODEL_H_INCLUDED


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
        std::vector<int>quantityFoodMaxTreeForSpecies,
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

);

#endif // OMNISCIENT_MODEL_H_INCLUDED
