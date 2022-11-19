#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <iostream>
class ProgressBar {
    public:
        void setIterationStartSize(int64_t _iterationStartSize){
            iterationStartSize = _iterationStartSize;
        }
        void setAlgorithmStartSize(int64_t _algorithmStartSize){
            algorithmStartSize = _algorithmStartSize;
        }
        #define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
        #define PBWIDTH 60

        void printIterationBar(int64_t batchProgress){
            float percentage = 1.0-(float(batchProgress)/float(iterationStartSize));
            int val = (int) (percentage * 100);
            int lpad = (int) (percentage * PBWIDTH);
            int rpad = PBWIDTH - lpad;
            printf("\rIteration %3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
            fflush(stdout);
        }
        void printAlgorithmBar(int64_t algorithmProgress){
            float percentage = 1.0-(float(algorithmProgress)/float(algorithmStartSize));
            int val = (int) (percentage * 100);
            int lpad = (int) (percentage * PBWIDTH);
            int rpad = PBWIDTH - lpad;
            printf("\rAlgorithm %3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
            fflush(stdout);
        }
    private:
        int64_t iterationStartSize;
        int64_t algorithmStartSize;

};
#endif /* PROGRESSBAR_H */