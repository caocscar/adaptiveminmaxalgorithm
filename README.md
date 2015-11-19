# adaptiveminmaxalgorithm
Algorithm for the adaptive min max methodology for background subtraction in Raman spectroscopy

## Basic Usage Instructions
1. Create a file called **filelist.txt** with names of all the files you want to process. 
2. Set graph flag to **1** (on) or **0** (off) in the m-file depending if you want to graph results or not.
3. Run matlab script.

## Options
As mentioned in the paper, you might want to change the **FSratio** threshold(s) for which polynomial order to
use based on your data and results. These thresholds are hardcoded in the *findorder* function. 
The defaults are [0.2, 0.75, 8.5, 55, 240].

## Reference
Cao, A, Pandya, AK, Serhatkulu, GK, Weber, RE, Dai, H, Thakur, JS, 
Naik, VM, Naik, R, Auner, GW, Rabah, R, Freeman, DC (2007).
"A robust method for automated background subtraction of tissue fluorescence." 
Journal of Raman Spectroscopy 38(9): 1199-1205. DOI: 10.1002/jrs.1753