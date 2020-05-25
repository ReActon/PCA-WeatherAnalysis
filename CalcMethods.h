#include <string>
#include <vector>
#include <fstream>

#ifndef CALCMETHODS
#define CALCMETHODS

//global var declarations
float xMean;
float yMean;

float getMean(const std::vector<float> &vals);

float calculateCov(const std::vector<float> &firstIn, const std::vector<float> &secondIn);

float subtractX(float val);

float subtractY(float val);

Eigen::MatrixXd formMatrix(const std::vector<float> &xComps, const std::vector<float> &yComps);

void writetoF(Eigen::MatrixXd &matrix);


#endif
