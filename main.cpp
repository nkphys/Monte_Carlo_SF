#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
using namespace std;
#include "Matrix.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "MCEngine.h"
#include "random"


int main(int argc, char *argv[]) {
    if (argc<2) { throw std::invalid_argument("USE:: executable inputfile"); }
    string inputfile = argv[1];

    bool check_Non_Int=false;

    Parameters Parameters_;
    Parameters_.Initialize(inputfile);

    Coordinates Coordinates_(Parameters_.lx, Parameters_.ly);

    mt19937_64 Generator_(Parameters_.RandomSeed);
    MFParams MFParams_(Parameters_,Coordinates_,Generator_);

    Hamiltonian Hamiltonian_(Parameters_,Coordinates_,MFParams_);

    Observables Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);




    if(check_Non_Int==true){
    //Parameters_.J_HUND=0.0;
   //Hamiltonian_.InteractionsCreate();
  //  Hamiltonian_.Ham_.print();
  // Hamiltonian_.Check_up_down_symmetry();
   //Hamiltonian_.Check_Hermiticity();
   // Hamiltonian_.Diagonalize('V');
   // Observables_.Get_Non_Interacting_dispersion();
    //Hamiltonian_.Ham_.print();
    // Observables_.Calculate_Akw();

    }


        MCEngine MCEngine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);

        Observables_.Initialize();     // Set All Observables to zero

        MCEngine_.RUN_MC();      // Monte-Carlo Engine


   

    cout << "--------THE END--------" << endl;
} // main
