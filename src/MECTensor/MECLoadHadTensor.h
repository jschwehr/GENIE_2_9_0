//____________________________________________________________________________
/*!

\class    genie::INukeHadroData

\brief    Singleton class to load & save mec hadron tensor tables - heavily
          based on MECLoadXSecFiles (which is based on INukeHadroData.cxx)

\author   Jackie Schwehr

\created  September 12, 2014

\notes    January 14, 2015 - Generalized for all neutrino flavors, Oxygen and Carbon


*/
//_____________________________________________________________________

#ifndef _MEC_LOAD_HAD_TENSOR_H_
#define _MEC_LOAD_HAD_TENSOR_H_

#include "Numerical/BLI2D.h"
#include "GHEP/GHepParticle.h"

namespace genie {

class MECLoadHadTensor
{
 public:
  static MECLoadHadTensor * Instance (int, int); // targetpdg, nupdg
  
  // public functions

  double XSecFullAll(int, int, double, double, double);
  double XSecFullpn(int, int, double, double, double) ;
  double XSecDeltaAll(int, int, double, double, double);
  double XSecDeltapn(int, int, double, double, double); 
  double MaxXSecAll(int, int, double);// targetpdg, nupdg, Enu
  double MaxXSecDelta(int, int, double); // targetpdg, nupdg, Enu

  double TotalXsecAtE(int targetpdg, int nupdg, double Enu);
  double GetTmuCostFromq0q3(double dq0, double dq3, double Enu, double lmass, double &tmu, double &cost, double &area);

  // Public Variables

  // hadron tensors
  // vector<grid> = five tensors (tensor[0-4].  vector< vector<grid>> four sets of five tables (tensor[0-3][0-4])
  const vector <genie::BLI2DNonUnifGrid *> HadTensorFullAllC12 (void) const {return HadTensorFullAll_C12_2DGrids;} // carbon
  const vector <genie::BLI2DNonUnifGrid *> HadTensorFullpnC12 (void) const {return HadTensorFullpn_C12_2DGrids;} // carbon
  const vector <genie::BLI2DNonUnifGrid *> HadTensorDeltaAllC12 (void) const {return HadTensorDeltaAll_C12_2DGrids;} // carbon
  const vector <genie::BLI2DNonUnifGrid *> HadTensorDeltapnC12 (void) const {return HadTensorDeltapn_C12_2DGrids;} // carbon

  const vector <genie::BLI2DNonUnifGrid *> HadTensorFullAllO16 (void) const {return HadTensorFullAll_O16_2DGrids;} // oxygen
  const vector <genie::BLI2DNonUnifGrid *> HadTensorFullpnO16 (void) const {return HadTensorFullpn_O16_2DGrids;} // oxygen
  const vector <genie::BLI2DNonUnifGrid *> HadTensorDeltaAllO16 (void) const {return HadTensorDeltaAll_O16_2DGrids;} // oxygen
  const vector <genie::BLI2DNonUnifGrid *> HadTensorDeltapnO16 (void) const {return HadTensorDeltapn_O16_2DGrids;} // oxygen

  
  const vector <genie::BLI2DNonUnifGrid *> HadTensorFullAllCa40 (void) const {return HadTensorFullAll_Ca40_2DGrids;} // calcium
  const vector <genie::BLI2DNonUnifGrid *> HadTensorFullpnCa40 (void) const {return HadTensorFullpn_Ca40_2DGrids;} // calcium
  const vector <genie::BLI2DNonUnifGrid *> HadTensorDeltaAllCa40 (void) const {return HadTensorDeltaAll_Ca40_2DGrids;} // calcium
  const vector <genie::BLI2DNonUnifGrid *> HadTensorDeltapnCa40 (void) const {return HadTensorDeltapn_Ca40_2DGrids;} // calcium
  
  // max xsec vectors
  const vector <double> EnuValues (void) const {return Enuvect;}

  const vector <double> MaxXSecAllC12v12 (void) const {return MaxXSecAllvectC12v12;}
  const vector <double> MaxXSecAllC12v14 (void) const {return MaxXSecAllvectC12v14;}
  const vector <double> MaxXSecAllC12v16 (void) const {return MaxXSecAllvectC12v16;}
  const vector <double> MaxXSecAllC12av12 (void) const {return MaxXSecAllvectC12av12;}
  const vector <double> MaxXSecAllC12av14 (void) const {return MaxXSecAllvectC12av14;}
  const vector <double> MaxXSecAllC12av16 (void) const {return MaxXSecAllvectC12av16;}

  const vector <double> MaxXSecDeltaC12v12 (void) const {return MaxXSecDeltavectC12v12;}
  const vector <double> MaxXSecDeltaC12v14 (void) const {return MaxXSecDeltavectC12v14;}
  const vector <double> MaxXSecDeltaC12v16 (void) const {return MaxXSecDeltavectC12v16;}
  const vector <double> MaxXSecDeltaC12av12 (void) const {return MaxXSecDeltavectC12av12;}
  const vector <double> MaxXSecDeltaC12av14 (void) const {return MaxXSecDeltavectC12av14;}
  const vector <double> MaxXSecDeltaC12av16 (void) const {return MaxXSecDeltavectC12av16;}

  const vector <double> MaxXSecAllO16v12 (void) const {return MaxXSecAllvectO16v12;}
  const vector <double> MaxXSecAllO16v14 (void) const {return MaxXSecAllvectO16v14;}
  const vector <double> MaxXSecAllO16v16 (void) const {return MaxXSecAllvectO16v16;}
  const vector <double> MaxXSecAllO16av12 (void) const {return MaxXSecAllvectO16av12;}
  const vector <double> MaxXSecAllO16av14 (void) const {return MaxXSecAllvectO16av14;}
  const vector <double> MaxXSecAllO16av16 (void) const {return MaxXSecAllvectO16av16;}

  const vector <double> MaxXSecDeltaO16v12 (void) const {return MaxXSecDeltavectO16v12;}
  const vector <double> MaxXSecDeltaO16v14 (void) const {return MaxXSecDeltavectO16v14;}
  const vector <double> MaxXSecDeltaO16v16 (void) const {return MaxXSecDeltavectO16v16;}
  const vector <double> MaxXSecDeltaO16av12 (void) const {return MaxXSecDeltavectO16av12;}
  const vector <double> MaxXSecDeltaO16av14 (void) const {return MaxXSecDeltavectO16av14;}
  const vector <double> MaxXSecDeltaO16av16 (void) const {return MaxXSecDeltavectO16av16;}

  
  const vector <double> MaxXSecAllCa40v12 (void) const {return MaxXSecAllvectCa40v12;}
  const vector <double> MaxXSecAllCa40v14 (void) const {return MaxXSecAllvectCa40v14;}
  const vector <double> MaxXSecAllCa40v16 (void) const {return MaxXSecAllvectCa40v16;}
  const vector <double> MaxXSecAllCa40av12 (void) const {return MaxXSecAllvectCa40av12;}
  const vector <double> MaxXSecAllCa40av14 (void) const {return MaxXSecAllvectCa40av14;}
  const vector <double> MaxXSecAllCa40av16 (void) const {return MaxXSecAllvectCa40av16;}

  const vector <double> MaxXSecDeltaCa40v12 (void) const {return MaxXSecDeltavectCa40v12;}
  const vector <double> MaxXSecDeltaCa40v14 (void) const {return MaxXSecDeltavectCa40v14;}
  const vector <double> MaxXSecDeltaCa40v16 (void) const {return MaxXSecDeltavectCa40v16;}
  const vector <double> MaxXSecDeltaCa40av12 (void) const {return MaxXSecDeltavectCa40av12;}
  const vector <double> MaxXSecDeltaCa40av14 (void) const {return MaxXSecDeltavectCa40av14;}
  const vector <double> MaxXSecDeltaCa40av16 (void) const {return MaxXSecDeltavectCa40av16;}


 private:

  // private offical stuff
  MECLoadHadTensor(int, int); // targetpdg, nupdg
  MECLoadHadTensor(const MECLoadHadTensor & shx);
  ~MECLoadHadTensor();
  static MECLoadHadTensor * fInstance;

  // private functions
  void LoadTensorTables(int targetpdg);
  void ReadHadTensorqzq0File(string filename, int nwpoints, int nqzpoints, int nq0points, double hadtensor_w_array[][57600]);
  //void MakeMaxXSecTables(int, int); //targetpdg, nupdg
  void WriteMaxXSecTables(int, int); // targetpdg, nupdg
  void ReadMaxXSecTables(int, int); // targetpdg, nupdg
  double MaxXSec(int, int, double, vector<double>); // targetpdg, nupdg, Enu, max xsec
  double XSec(int, int, double, double, double, vector <BLI2DNonUnifGrid *>); // targetpdg, nupdg, E, T, Costheta, HadTensor;


  // private varables

  // Hadron Tensors
  vector <genie::BLI2DNonUnifGrid *> HadTensorFullAll_C12_2DGrids; // carbon
  vector <genie::BLI2DNonUnifGrid *> HadTensorFullpn_C12_2DGrids; // carbon
  vector <genie::BLI2DNonUnifGrid *> HadTensorDeltaAll_C12_2DGrids; // carbon
  vector <genie::BLI2DNonUnifGrid *> HadTensorDeltapn_C12_2DGrids; // carbon

  vector <genie::BLI2DNonUnifGrid *> HadTensorFullAll_O16_2DGrids; // oxygen
  vector <genie::BLI2DNonUnifGrid *> HadTensorFullpn_O16_2DGrids; // oxygen
  vector <genie::BLI2DNonUnifGrid *> HadTensorDeltaAll_O16_2DGrids; // oxygen
  vector <genie::BLI2DNonUnifGrid *> HadTensorDeltapn_O16_2DGrids; // oxygen

  vector <genie::BLI2DNonUnifGrid *> HadTensorFullAll_Ca40_2DGrids; // oxygen
  vector <genie::BLI2DNonUnifGrid *> HadTensorFullpn_Ca40_2DGrids; // oxygen
  vector <genie::BLI2DNonUnifGrid *> HadTensorDeltaAll_Ca40_2DGrids; // oxygen
  vector <genie::BLI2DNonUnifGrid *> HadTensorDeltapn_Ca40_2DGrids; // oxygen


  // Vectors
  vector <double> Enuvect;

  vector <double> MaxXSecAllvectC12v12;
  vector <double> MaxXSecAllvectC12v14;
  vector <double> MaxXSecAllvectC12v16;
  vector <double> MaxXSecAllvectC12av12;
  vector <double> MaxXSecAllvectC12av14;
  vector <double> MaxXSecAllvectC12av16;

  vector <double> MaxXSecDeltavectC12v12;
  vector <double> MaxXSecDeltavectC12v14;
  vector <double> MaxXSecDeltavectC12v16;
  vector <double> MaxXSecDeltavectC12av12;
  vector <double> MaxXSecDeltavectC12av14;
  vector <double> MaxXSecDeltavectC12av16;

  vector <double> MaxXSecAllvectO16v12;
  vector <double> MaxXSecAllvectO16v14;
  vector <double> MaxXSecAllvectO16v16;
  vector <double> MaxXSecAllvectO16av12;
  vector <double> MaxXSecAllvectO16av14;
  vector <double> MaxXSecAllvectO16av16;

  vector <double> MaxXSecDeltavectO16v12;
  vector <double> MaxXSecDeltavectO16v14;
  vector <double> MaxXSecDeltavectO16v16;
  vector <double> MaxXSecDeltavectO16av12;
  vector <double> MaxXSecDeltavectO16av14;
  vector <double> MaxXSecDeltavectO16av16;

  vector <double> MaxXSecAllvectCa40v12;
  vector <double> MaxXSecAllvectCa40v14;
  vector <double> MaxXSecAllvectCa40v16;
  vector <double> MaxXSecAllvectCa40av12;
  vector <double> MaxXSecAllvectCa40av14;
  vector <double> MaxXSecAllvectCa40av16;

  vector <double> MaxXSecDeltavectCa40v12;
  vector <double> MaxXSecDeltavectCa40v14;
  vector <double> MaxXSecDeltavectCa40v16;
  vector <double> MaxXSecDeltavectCa40av12;
  vector <double> MaxXSecDeltavectCa40av14;
  vector <double> MaxXSecDeltavectCa40av16;



  // singleton cleaner
  struct Cleaner {
    void DummyMethodAndSilentCompiler(){}
    ~Cleaner(){
      if (MECLoadHadTensor::fInstance !=0){
	delete MECLoadHadTensor::fInstance;
	MECLoadHadTensor::fInstance = 0;
      }
    }
  };
  friend struct Cleaner;
};

} // genie namespace

#endif // _MEC_LOAD_HAD_TENSOR_H_
