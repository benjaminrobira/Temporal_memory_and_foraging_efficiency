# Temporal_memory_and_foraging :brain: 

*Note: English is not my native language. Hopefully, the following text will be understandable and without too many English errors.*

This repository contains codes and data for reproducing simulations of the manuscript: "*Foraging efficiency in temporally predictable environments: 
Is a long-term temporal memory really advantageous?*" by Benjamin Robira, Simon Benhamou, Shelly Masi, Violaine Llaurens and Louise Riotte-Lambert. 

:newspaper: The associated manuscript has been submitted and is currently under review. Hence, the code has gone through only a limited cleaning so as to allow readability. It might still include commented lines for explored - but not included in the manuscript (e.g. the competitive "Scenario 3") - parts, in case it is useful for the reviewing process.
The structure is not the best - it was done at the beginning of coding apprenticeship; my apology for that. Hopefully, once the manuscript is accepted, it will go through a substantial "layout" editing.

:file_folder: You will find two folders: scripts, with all the necessary codes to run the model, and output which corresponds to the data from our simulations. Normally, we used a fixed seed, so running the code should give the same output (yet potentially, results were individually re-run by scenario during revision, so that might not be systematic). Overall the low variability between runs should inevitably lead to similar results.

:computer: "Script": The code is implemented in C++ language. A codeblock project was created (`Memory_time_model.cbp`; .depend and .layout files are also associated), with Code::Blocks software (v 17.12) used as the programming interface. It stands as follows:

* `main.cpp`: the main function that, if compiled and run, will run the model for homogeneous/heterogeneous environment in time and space, as well as the "Delay experiment" for the chronological memory presented in Supplementary Material.
* `MEMORY_model`: the files with the code (.cpp) or the head (.h) that simulate the behaviour of a MEMORY type. The names are mostly self explicit: Associative and Chronological are named as in the paper (the correcting mechanism was added with the "delay" extension for the Chronological memory). The Omniscient memory corresponds to what is named "Prescient" in the memory (because we then estimated that omniscience was related to the spatial component, prescience to the temporal one). The "Null_model" is the "Basic" agent referred in the manuscript (with only a no-return memory avoiding back-tracking). Finally, the "Null-working-memory" corresponds to the "Sampler".
* `Food_calculation`: the file with the code (.cpp) or the head (.h) that implements the function to calculate currently available food at a given time and location given the phenology of the resource type, and previous history of the resource location. 
* `Confidence_estimation`: the file with the code (.cpp) or the head (.h) that implements the function for the agent endowed with a Chronological or Associative memory to estimate the fruiting probability at the species level based on private knowledge of the start/end of the productive period.

:chart_with_upwards_trend: "Output": The data are available for each tested scenario as zip. Each run has its own .txt output. Our interest was in the Foraging_success_tot variable (not the cumulative one, which summed the foraging success of the movement bouts)! Other variables are self explicit (and most of them are not of interest for reproduction of the work). Note that the synchrony is labelled between 0 and 1 while we speak about the desynchrony (and plotted with desynchrony) in the manuscript, where the relationships between synchrony (constrained between 0 and 1) with the desynchrony is explicited.
Hopefully, I did not forget any table or mislabelled file when merging everything!

:question::exclamation: You have a question? The code has got an error? I am disappointed it was not clear and/or incorrect. Please let me know at [benjamin.robira@normalesup.org](mailto:benjamin.robira@normalesup.org) and I hope we will be able to solve it together! 
