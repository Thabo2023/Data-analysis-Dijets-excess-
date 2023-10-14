#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>

void CLTOPVAL() {
   // Read mass and CL values from the text file
   std::ifstream file("atlassig33.txt");
   std::vector<double> masses;
   std::vector<double> confidenceLevels;

   double mass, confidenceLevel;
   while (file >> mass >> confidenceLevel) {
       masses.push_back(mass);
       confidenceLevels.push_back(confidenceLevel);
   }
   file.close();

   // Loop through the data points and calculate the p-values
  
   int numPoints = masses.size();
   for (int i = 0; i < numPoints; ++i) {
       double mass = masses[i];
       double significance = confidenceLevels[i];

   
       double z = significance;
      

       cout<< mass <<" " << 2*(1 - ROOT::Math::normal_cdf(z))<< endl;
  
}

}

