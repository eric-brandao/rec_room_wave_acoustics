/*This is the C++ file to calculate the room FRF according to modal superposition
It is valid for a shoe-box room with hard walls - Calculations are done in frequency domain
%%% Author: Eric Brandão
%%% Last significant modification: 04/12/2018


%% Input variables
%%% general - all general stuff about the algorithm
%%% geometry - geometry of the room - vertexes and polygons
%%% sources - source properties
%%% receivers - point receivers
%%% results - will pass back to matlab
*/

// In the following, my interpretation of what is happening.
/* The C++ file starts with the definition of a class gateway named "class MexFunction" that has an "operator" void function.
This is the doorway though wich data comes from matlab workspace (inputs) and comes back (outputs)
This gateway class function will call other functions that acctually do stuff. For instance, the function "checkArguments"
perform a check on the type of input arguments. This is important because if the user provide inputs that the C++ file did
not foresee, there may be problems such as crash downs.
*/
#include "mex.hpp"
#include "mexAdapter.hpp"
#include<string>
#include<memory>
#include <iostream>
#include <vector>
#include <ctime>
#include <complex>      // std::complex, std::abs
//#include <cmath>
#include <D:\Work\UFSM\toolbox_ufsm\eigen-eigen-b3f3d4950030\Eigen/Dense>
#include "wavefrf.h"


// To simplify typing.
using matlab::mex::ArgumentList;
using namespace matlab::data;
using namespace matlab::engine;
using namespace std;

# define M_PI           3.14159265358979323846  /* pi */

class MexFunction : public matlab::mex::Function {
public:
	void operator()(ArgumentList outputs, ArgumentList inputs) {

		// check arguments and throw errors if necessary
		checkArguments(outputs, inputs);

		// get struct arrays to pass to the main function
		ArrayFactory factory;
		StructArray general(inputs[0]);
		StructArray geometry(inputs[1]);
		StructArray sources(inputs[2]);
		StructArray receivers(inputs[3]);
		StructArray results(inputs[4]); 

		// Call calculation function
		FRF_calculator(general, geometry, sources, receivers, results);

		// Pass data back to a struct in matlab
		outputs[0] = results;
	} // Operator ends here

	// Check arguments function - important to do in C++ and return errors
	void checkArguments(ArgumentList outputs, ArgumentList inputs) {
		std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
		ArrayFactory factory;
		StructArray general(inputs[0]);
		StructArray geometry(inputs[1]);
		StructArray sources(inputs[2]);
		StructArray receivers(inputs[3]);
		StructArray results(inputs[4]);

		// Check the number of inputs
		if (inputs.size() != 5) {
			matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
				0, std::vector<Array>({ factory.createScalar("Five inputs required in the following order: general, geometry, sources, receivers and results (that comes mostly empty)") }));
		}
		// Check the number of outputs
		if (outputs.size() > 1) {
			matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
				0,
				std::vector<Array>({ factory.createScalar("Too many outputs specified") }));
		}
		// Check if the inputs are all structures
		for (int i = 0; i < 4; i++) {
			if (inputs[i].getType() != ArrayType::STRUCT) {
				matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
					0, std::vector<Array>
					({ factory.createScalar("Input must be a structure") }));
			}
		}

		// Let's check inputs for Dimensions and Number of Fields: (general)
		auto numelem0 = general.getDimensions(); // get Dimensions
		size_t numFields0 = general.getNumberOfFields(); // get Number of Fields
		if (numelem0[0] != 1 || numelem0[1] != 1 || numFields0 != 14) {
			matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
				0, std::vector<Array>
				({ factory.createScalar("Variable (general) must be a 1x1 structure with 14 Fields: Check comments on original C++ code.") }));
		}

		// Let's check inputs for Dimensions and Number of Fields: (geometry)
		auto numelem1 = geometry.getDimensions(); // get Dimensions
		size_t numFields1 = geometry.getNumberOfFields(); // get Number of Fields
		if (numelem1[0] != 1 || numelem1[1] != 1 || numFields1 != 8) {
			matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
				0, std::vector<Array>
				({ factory.createScalar("Variable (geometry) must be a 1x1 structure with 8 Fields: Check comments on original C++ code.") }));
		}

		// Let's check inputs for Dimensions and Number of Fields: (sources)
		auto numelem2 = sources.getDimensions(); // get Dimensions
		size_t numFields2 = sources.getNumberOfFields(); // get Number of Fields
		if (numelem2[0] != 1 || numFields2 != 2) {
			matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
				0, std::vector<Array>
				({ factory.createScalar("Variable (sources) must be a 1xN structure array with 2 Fields: Check comments on original C++ code.") }));
		}

		// Let's check inputs for Dimensions and Number of Fields: (receivers)
		auto numelem3 = receivers.getDimensions(); // get Dimensions
		size_t numFields3 = receivers.getNumberOfFields(); // get Number of Fields
		if (numelem3[0] != 1 || numFields3 != 1) {
			matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
				0, std::vector<Array>
				({ factory.createScalar("Variable (receivers) must be a 1xNrec structure array with 1 Field: Check comments on original C++ code.") }));
		}

		// Let's check inputs for Dimensions and Number of Fields: (receivers)
		auto numelem4 = results.getDimensions(); // get Dimensions
		size_t numFields4 = results.getNumberOfFields(); // get Number of Fields
		if (numelem4[0] != 1 || numelem4[1] != 1 || numFields4 != 2) {
			matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
				0, std::vector<Array>
				({ factory.createScalar("Variable (results) must be a 1x1 structure  with 2 Fields: Check comments on original C++ code.") }));
		}

	} // Check arguments function ends here 

	// Here things really start - Here is all your calculations.
	void FRF_calculator(StructArray general, StructArray geometry, StructArray sources, StructArray receivers, StructArray& results) {
		// Initialization
		std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
		ArrayFactory factory;

		// General stuff
		auto generalfields = general.getFieldNames(); // Get the field names in the variable general
		std::vector<MATLABFieldIdentifier> generalfieldNames(generalfields.begin(), generalfields.end());
		// Gets arrays associated with "general" 
		const TypedArray<double> rho0 = general[0][generalfieldNames[3]];
		const TypedArray<double> c0 = general[0][generalfieldNames[4]];
		double c00 = c0[0]; // Helps complex
		const TypedArray<double> f_max = general[0][generalfieldNames[5]];
		
		const TypedArray<double> Fs = general[0][generalfieldNames[7]];
		const TypedArray<double> Dt = general[0][generalfieldNames[8]];
		const TypedArray<double> NFFT = general[0][generalfieldNames[9]];
		const TypedArray<double> freq = general[0][generalfieldNames[10]]; // Isto por exemplo é um array de frequências. Seria legal ter ele como Eigen.
		const TypedArray<double> idf = general[0][generalfieldNames[13]];
		int idfc = int(idf[0]);
		//Eigen::RowVectorXd Psi_xs(fnSize, 1);
		//TypedArray<Eigen::RowVectorXd> freq = general[0][generalfieldNames[10]];
		size_t freqSize = freq.getNumberOfElements();
		
		// Geometry stuff
		auto geometryfields = geometry.getFieldNames(); // Get the field names in the variable geometry
		std::vector<MATLABFieldIdentifier> geometryfieldNames(geometryfields.begin(), geometryfields.end());
		// Gets arrays associated with "geometry"
		const TypedArray<double> Lx = geometry[0][geometryfieldNames[0]];
		const TypedArray<double> Ly = geometry[0][geometryfieldNames[1]];
		const TypedArray<double> Lz = geometry[0][geometryfieldNames[2]];
		const TypedArray<double> T60 = geometry[0][geometryfieldNames[3]];
		const TypedArray<double> Volume = geometry[0][geometryfieldNames[4]];
		double Bn = 6.91 / T60[0];
		
		// Sources stuff
		auto sourcesfields = sources.getFieldNames(); // Get the field names in the variable sources
		std::vector<MATLABFieldIdentifier> sourcesfieldNames(sourcesfields.begin(), sourcesfields.end());
		// Gets arrays associated with "sources" - other relevant data is got inside for loops
		size_t Ns = sources.getNumberOfElements(); // Get the number of sources supplied by the user
		
		// Receivers's stuff
		auto receiversfields = receivers.getFieldNames(); // Get the field names in the variable receivers
		std::vector<MATLABFieldIdentifier> receiversfieldNames(receiversfields.begin(), receiversfields.end());
		// Gets arrays associated with "receivers" - other relevant data is got inside for loops
		size_t Nrec = receivers.getNumberOfElements(); // Get the number of receivers supplied by the user
		
	
		// Results stuff - This comes empty from matlab. Once the calculations are done the mex function fills the appropriate data and pass it back to matlab
		auto resultsfields = results.getFieldNames(); // Get the field names in the variable results
		std::vector<MATLABFieldIdentifier> resultsfieldNames(resultsfields.begin(), resultsfields.end());
		// Gets empty arrays and cells associated  with "results" - Only fn_table by now
		TypedArray<double> fn_table = results[0][resultsfieldNames[0]];
		size_t fnSize = fn_table.getNumberOfElements()/7;
		CellArray HwCell = std::move(results[0][resultsfieldNames[1]]);

		// Complex number definition
		const complex<double> i1(0.0, 1.0);
		
		// Copy the fn_table data to Eigen vectors (for future calculation)
		cout << "I am coppying data do Eigen" << endl;
		//typedef Eigen::Matrix<double, Eigen::Dynamic, 1> ArrayfnD;
		Eigen::MatrixXd fn(fnSize, 1);
		Eigen::MatrixXd nx(fnSize, 1);
		Eigen::MatrixXd ny(fnSize, 1);
		Eigen::MatrixXd nz(fnSize, 1);
		Eigen::MatrixXd Axyz(fnSize, 1);

		//cout << "size of nx is: (" << nx.rows() << ", " << nx.cols() << ");" << endl;

		for (int jfn = 0; jfn < fnSize; jfn++) {
			fn(jfn, 0) = fn_table[jfn][1];
			nx(jfn, 0) = fn_table[jfn][2];
			ny(jfn, 0) = fn_table[jfn][3];
			nz(jfn, 0) = fn_table[jfn][4];
			Axyz(jfn, 0) = fn_table[jfn][5];

		}
		Eigen::MatrixXcd kn2i(fnSize, 1); // collumn array
		kn2i = pow(2 * M_PI * fn.array() / c00, 2.0) + 2.0 * i1 * (2.0 * M_PI * fn.array() / c00) * Bn / c00;
		cout << "I finished coppying data do Eigen" << endl;
		//cout << "size of nx is: (" << nx.rows() << ", " << nx.cols() << ");" << endl;
		//cout << "size of Axyz is: (" << Axyz.rows() << ", " << Axyz.cols() << ");" << endl;
		//cout << "nx is: " << nx << endl;
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% /////////////////////
		// Now, we will calculate de spectrum of sound pressure.
		// Step 1: Calculate the modeal functions
		Eigen::MatrixXd Psi_Xr(fnSize, 1); // collumn array
		Eigen::MatrixXd Psi_Xs(fnSize, 1); // collumn array

		//Eigen::MatrixXcd Psi_XrM(freqSize,fnSize); // matrix


		for (size_t js = 0; js < Ns; js++) { // Loop over sources
			
			TypedArray<double> SourceCoord = std::move(sources[js][sourcesfieldNames[0]]); // Current source coord
			TypedArray<double> SourceQ = std::move(sources[js][sourcesfieldNames[1]]); // Current source coord
			double rho0QV = rho0[0] * SourceQ[0] / Volume[0];

			Psi_Xs.array() = Axyz.array() * cos(nx.array() * M_PI * SourceCoord[0] / Lx[0])  * cos(ny.array() * M_PI * SourceCoord[1] / Ly[0])
				* cos(nz.array() * M_PI * SourceCoord[2] / Lz[0]);

			for (size_t jrec = 0; jrec < Nrec; jrec++) { // Loop over receivers
				cout << "Performing calculations for source: " << js + 1 << "; and receiver: " << jrec+1 << endl;
				TypedArray<double> RecCoord = std::move(receivers[jrec][receiversfieldNames[0]]); // Current receiver coord

				Psi_Xr.array() = Axyz.array() * cos(nx.array() * M_PI * RecCoord[0] / Lx[0])  * cos(ny.array() * M_PI * RecCoord[1] / Ly[0])
					* cos(nz.array() * M_PI * RecCoord[2] / Lz[0]);

				// Loop over frequencies
				int posplot = 0;
				Array HwVec;
				HwVec = factory.createArray<complex<double>>({ freqSize, 1 }); // Create a Matlab array with the size of TimeCrossVec to pass to Matlab
																			   //Eigen::MatrixXcd HwE(freqSize, 1); // collumn array

				for (size_t jf = 0; jf < freqSize; jf++) {
					//DoProgress(jf, freqSize, posplot);
					DoProgress(jf, idfc, posplot);
					if (freq[jf] <= 1.5*f_max[0]) {
						
						Eigen::MatrixXcd Psi_XrC(fnSize, 1); // collumn array
						Psi_XrC.array() = Psi_Xr.array() / (pow(2 * M_PI * freq[jf] / c00, 2.0) - kn2i.array());

						Eigen::MatrixXcd Psi_ProdSum(1, 1);
						Psi_ProdSum = Psi_Xs.transpose() * Psi_XrC; // Product and sum of modal functions

						HwVec[jf] = i1 * rho0QV * (2.0 * M_PI * freq[jf]) * Psi_ProdSum(0);
					}
					//else
					//HwVec[jf] = 0.0 + i1 * 0.0;
				}

				// Matrix multiplication
				//HwE = Psi_XrM * Psi_Xr;
				//cout << HwE << endl;
				// Translate back to matlab data
				//for (size_t jf = 0; jf < freqSize; jf++) {
				//	HwVec[jf] = i1 * rho0QV * (2.0 * M_PI * freq[jf]) * HwE(jf);
				//}

				// pass spectrum to matlab
				HwCell[js][jrec] = HwVec;

			} // end of receivers loop
		} // end of sources loop 

		results[0][resultsfieldNames[1]] = std::move(HwCell);

		cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << endl;
		cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << endl;
		cout << "I'm done! Hope you like it!" << endl;

	} // End of the main function
}; // Mex Class ends here

/*
// Copy the fn_table data to Eigen vectors (for future calculation)
cout << "I am coppying data do Eigen" << endl;
//typedef Eigen::Matrix<double, Eigen::Dynamic, 1> ArrayfnD;
Eigen::MatrixXd fn(1, fnSize);
Eigen::MatrixXd nx(1, fnSize);
Eigen::MatrixXd ny(1, fnSize);
Eigen::MatrixXd nz(1, fnSize);
Eigen::MatrixXd Axyz(1, fnSize);

//cout << "size of nx is: (" << nx.rows() << ", " << nx.cols() << ");" << endl;

for (int jfn = 0; jfn < fnSize; jfn++) {
fn(0, jfn) = fn_table[jfn][1];
nx(0, jfn) = fn_table[jfn][2];
ny(0, jfn) = fn_table[jfn][3];
nz(0, jfn) = fn_table[jfn][4];
Axyz(0, jfn) = fn_table[jfn][5];

}
Eigen::MatrixXcd kn2i(1, fnSize); // collumn array
kn2i = pow(2 * M_PI * fn.array() / c00, 2.0) + 2.0 * i1 * (2.0 * M_PI * fn.array() / c00) * Bn / c00;
cout << "I finished coppying data do Eigen" << endl;
//cout << "size of nx is: (" << nx.rows() << ", " << nx.cols() << ");" << endl;
//cout << "size of Axyz is: (" << Axyz.rows() << ", " << Axyz.cols() << ");" << endl;
//cout << "nx is: " << nx << endl;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% /////////////////////
// Now, we will calculate de spectrum of sound pressure.
// Step 1: Calculate the modeal functions
Eigen::MatrixXd Psi_Xr(1, fnSize); // row array
Eigen::MatrixXd Psi_Xs(1, fnSize); // row array

Eigen::MatrixXcd Psi_XsM(freqSize, fnSize); // matrix


for (size_t js = 0; js < Ns; js++) { // Loop over sources
cout << "Performing calculations for source: " << js + 1 << ";" << endl;
TypedArray<double> SourceCoord = std::move(sources[js][sourcesfieldNames[0]]); // Current source coord
TypedArray<double> SourceQ = std::move(sources[js][sourcesfieldNames[1]]); // Current source coord
double rho0QV = rho0[0] * SourceQ[0] / Volume[0];

Psi_Xs.array() = Axyz.array() * cos(nx.array() * M_PI * SourceCoord[0] / Lx[0])  * cos(ny.array() * M_PI * SourceCoord[1] / Ly[0])
* cos(nz.array() * M_PI * SourceCoord[2] / Lz[0]);

for (size_t jrec = 0; jrec < Nrec; jrec++) { // Loop over receivers
TypedArray<double> RecCoord = std::move(receivers[jrec][receiversfieldNames[0]]); // Current receiver coord

Psi_Xr.array() = Axyz.array() * cos(nx.array() * M_PI * RecCoord[0] / Lx[0])  * cos(ny.array() * M_PI * RecCoord[1] / Ly[0])
* cos(nz.array() * M_PI * RecCoord[2] / Lz[0]);

// Loop over frequencies
int posplot = 0;
Eigen::MatrixXcd HwE(freqSize, 1); // collumn array

for (size_t jf = 0; jf < freqSize; jf++) {
DoProgress(jf, freqSize, posplot);

Psi_XsM.row(jf) = Psi_Xs.array() / (pow(2 * M_PI * freq[jf] / c00, 2.0) - kn2i.array());


//Eigen::MatrixXcd Psi_XrC(fnSize, 1); // collumn array
//Psi_XrC.array() = Psi_Xr.array() / (pow(2 * M_PI * freq[jf] / c00, 2.0) - kn2i.array());

//Eigen::MatrixXcd Psi_ProdSum(1, 1);
//Psi_ProdSum = Psi_Xs.transpose() * Psi_XrC; // Product and sum of modal functions

//HwVec[jf] = i1 * rho0QV * (2.0 * M_PI * freq[jf]) * Psi_ProdSum(0);
}

// Matrix multiplication
cout << "using Eigen - multiplying matrix..." << endl;
HwE = Psi_XsM * Psi_Xr.transpose();
//cout << HwE << endl;
// Translate back to matlab data
cout << "Translating back to matlab..." << endl;
Array HwVec;
HwVec = factory.createArray<complex<double>>({ freqSize, 1 }); // Create a Matlab array with the size of TimeCrossVec to pass to Matlab
for (size_t jf = 0; jf < freqSize; jf++) {
HwVec[jf] = i1 * rho0QV * (2.0 * M_PI * freq[jf]) * HwE(jf);
}

// pass spectrum to matlab
HwCell[js][jrec] = HwVec;

} // end of receivers loop
} // end of sources loop

results[0][resultsfieldNames[1]] = std::move(HwCell);

cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << endl;
cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << endl;
cout << "I'm done! Hope you like it!" << endl;
*/

