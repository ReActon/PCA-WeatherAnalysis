/* Ryan Acton - Assignment 5*
 * PCA using Average Rainfall Data
 * Submitted 17/07/2020
 */
 
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "values.hpp"
#include "CalcMethods.h"

/* Method to calculate the mean from a given input vector */
float getMean(const std::vector<float> &vals) {
    float mean = 0;
    
     for (auto item: vals) {
         mean += item;
     }
    mean = mean/vals.size();
    return mean;
}

/* Functor-method to use with transform to subtract the x mean from each x val */
float subtractX(float val) {return (val - xMean);}

/* Functor-method to use with transform to subtract the y mean from each y val */
float subtractY(float val) {return (val - yMean);}


/* Method to calculate each covariance given x and y values from the separated vector inputs */
float calculateCov(const std::vector<float> &tfirstIn, const std::vector<float> &tsecondIn) {
    float total = 0;
    for (int i = 0; i < tfirstIn.size(); i++) {
        total += tfirstIn.at(i) * tsecondIn.at(i); 
    }
    total = total/(tfirstIn.size()-1);
    return total;
}

/* Method, using the Eigen library to calculate and return a covariance matrix, given x and y inputs */
Eigen::MatrixXd formMatrix(const std::vector<float> &xComps, const std::vector<float> &yComps) {
	
	Eigen::MatrixXd matrix(2,2);
	matrix(0, 0) = calculateCov(xComps, xComps);
	matrix(0, 1) = calculateCov(xComps, yComps);
	matrix(1, 0) = calculateCov(yComps, xComps);
	matrix(1, 1) = calculateCov(yComps, yComps);
	
	return matrix;
}

/* Method to write the values for eigenvalues, eigenvectors, covariance matrix and components to a file, given the matrix */
void writetoF(Eigen::MatrixXd &matrix) {
	
	std::string fileName = "answers.txt";
	std::ofstream outFile;
	outFile.open(fileName);
	
	Eigen::VectorXcd eivals = matrix.eigenvalues();
	
	Eigen::EigenSolver<Eigen::MatrixXd> cs(matrix);
	
	float totalVar = matrix(0, 0) + matrix(1, 1);
	
	if(outFile.is_open()) {
		
		//Formatting lines
		outFile << "Ryan Acton: Assignment 5 Answers" << std::endl;
		outFile << "--------------------------------------------" << std::endl << std::endl;
		
		//Eigenvalues for the PCs
		outFile << "Eigenvalues for Principal Component 1: ";
		outFile << eivals.row(0).real() << std::endl << std::endl;	//row 0 corresponds to first PC
		outFile << "Eigenvalues for Principal Component 2: ";
		outFile << eivals.row(1).real() << std::endl << std::endl;	//row 0 corresponds to second PC
		
		//Eigenvectors for the PCs
		outFile << "Eigenvectors for Principal Component 1: ";
		outFile << std::endl<< cs.eigenvectors().col(0).real() << std::endl << std::endl;
		outFile << "Eigenvectors for Principal Component 2: ";
		outFile << std::endl << cs.eigenvectors().col(1).real() << std::endl << std::endl;		
		
		//Covariance Matrix
		outFile << "Covariance Matrix is: " << std::endl;
		outFile << matrix << std::endl << std::endl;	
		
		//Total Variance
		outFile << "The Total Variance for the Covariance Matrix is: " << std::endl;
		outFile << totalVar << std::endl << std::endl;
		
		//Percentages
		outFile << "Principal Component 1 explains: " << std::endl;		//calculations divided the eigenvalue for the component by the variance and represent it as percentage
		outFile << ((eivals.row(0).real())/totalVar)*100 << "%" << std::endl << std::endl;  
		outFile << "Principal Component 2 explains: " << std::endl;
		outFile << ((eivals.row(1).real())/totalVar)*100 << "%" << std::endl;
		
		
	}
	
	outFile.close();

}

int main() {
	
	//empty vars to store the transformed vectors
    std::vector<float> xValResults;
    std::vector<float> yValResults;
	
	//setting the globale parameters to be used by the transform
	xMean = getMean(january);
	yMean = getMean(june);

    std::transform(january.begin(), january.end(), back_inserter(xValResults), subtractX);
    std::transform(june.begin(), june.end(), back_inserter(yValResults), subtractY);
	
	Eigen::MatrixXd matrix(2,2); 
	matrix = formMatrix(xValResults, yValResults);	//local copy of matrix to pass into the write method
	
	writetoF(matrix);

	
    return 0;
}
