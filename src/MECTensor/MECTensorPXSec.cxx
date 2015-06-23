/*

coppying heavily from MECPXSec.cxx



 */

// includes
#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/GBuild.h"
#include "Conventions/Units.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "MECTensor/MECTensorPXSec.h"
#include "MECTensor/MECLoadHadTensor.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"


using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________
MECTensorPXSec::MECTensorPXSec() :
  XSecAlgorithmI("genie::MECTensorPXSec")
{

}

//_____________________________________________________________
MECTensorPXSec::MECTensorPXSec(string config) :
  XSecAlgorithmI("genie::MECPSXec",config)
{

}

//_____________________________________________________________
MECTensorPXSec::~MECTensorPXSec()
{

}

//_____________________________________________________________
double MECTensorPXSec::XSec(
		     const Interaction * interaction, KinePhaseSpace_t kps) const
{
  // return double differential xsec

  // initial state
  int    tgtpdg = interaction->InitState().Tgt().Pdg();
  int    nupdg  = interaction->InitState().ProbePdg();
  double Enu    = interaction->InitState().ProbeE(kRfLab);
  TLorentzVector * v4Nu = interaction->InitState().GetProbeP4(kRfLab);

  // final state
  const Kinematics kinematics = interaction -> Kine();
  TLorentzVector v4lep = kinematics.FSLeptonP4();
  double Mlep = interaction->FSPrimLepton()->Mass();

  // kinematics 
  double Tmu = v4Nu->E()*v4Nu->E() - Mlep;
  double CosTheta = cos(v4lep.Theta() - v4Nu->Theta());
  
  // get hadron tensor
   MECLoadHadTensor * hadtensor = MECLoadHadTensor::Instance(tgtpdg, nupdg);

  double xsec = hadtensor->XSecFullAll(tgtpdg, nupdg, Enu, Tmu, CosTheta);
  //double xsec = XSecAllTargets(tgtpdg, nupdg, Enu, Tmu, CosTheta);

  return xsec;
}

//_____________________________________________________________
double MECTensorPXSec::Integral(const Interaction * interaction) const {
  // return integrated xsec at the given neutrino energy

  // initial state
  int    tgtpdg = interaction->InitState().Tgt().Pdg();
  int    nupdg  = interaction->InitState().ProbePdg();
  double Enu = interaction->InitState().ProbeE(kRfHitNucRest);

  // get hadron tensor
   MECLoadHadTensor * hadtensor = MECLoadHadTensor::Instance(tgtpdg, nupdg);

  double xSecAtE = hadtensor->TotalXsecAtE(tgtpdg, nupdg, Enu);
  //xSecAtE = TotalXSecAtEAllTargets(tgtpdg, nupdg, Enu);

  if (!xSecAtE) return 0;
  else return xSecAtE;  

}


//_____________________________________________________________
bool MECTensorPXSec::ValidProcess(const Interaction * interaction) const {
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info = interaction->ProcInfo();
  if(!proc_info.IsMECTensor()&&!proc_info.IsMECTensorPDD()) return false;

  else return true;
}

//_____________________________________________________________
void MECTensorPXSec::Configure(const Registry & config){
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_____________________________________________________________
void MECTensorPXSec::LoadConfig(void){

}

/*
//____________________________________________________________________
// extension to all targets:
double MECTensorPXSec::XSecAllTargets(int tgtpdg,int nupdg,double Enu,double Tmu,double CosTheta){
    
  double xsec;
  
  int pdghi, pdglo;
  int findPDG = findPDGHiLo(tgtpdg, &pdghi, &pdglo);

  // Have Target
  if (findPDG == 0){
    MECLoadHadTensor * hadtensor = MECLoadHadTensor::Instance(tgtpdg, nupdg);    
    xsec = hadtensor->XSecFullAll(tgtpdg, nupdg, Enu, Tmu, CosTheta);
  }

  // If don't have target, figure out if interpolate, extrapolate, or other
  else {
    
    // function had error:
    if (findPDG == -1){
      return -1;
    }
    // extrapolate above
    else if (pdghi == -1 && pdglo != -1 ){
      return -1;
    }
    // extrapolate below
    else if (pdghi != -1 && pdglo == -1 ){
      return -1;
    }
    
    // interpolate
    else if (pdghi != -1 && pdglo != -1 ){
      
      // load larger and smaller target tensors
      MECLoadHadTensor * hadtensorLo = MECLoadHadTensor::Instance(pdglo, nupdg);
      MECLoadHadTensor * hadtensorHi = MECLoadHadTensor::Instance(pdghi, nupdg);

      // pull xsecs 
      double xsecLo = hadtensorLo->XSecFullAll(pdglo, nupdg, Enu, Tmu, CosTheta);
      double xsecHi = hadtensorHi->XSecFullAll(pdghi, nupdg, Enu, Tmu, CosTheta);

      // linear interpolate
      xsec = interpolateAllTargets(tgtpdg, pdglo, pdghi, xsecLo, xsecHi);

    }// end else if interpolate
  }// end else if don't have target

  return xsec;
}

//____________________________________________________________________
// extension to all targets:
double MECTensorPXSec::TotalXSecAtEAllTargets(int tgtpdg, int nupdg, double Enu){

  double xSecAtE;
  int pdghi, pdglo;

  int findPDG = findPDGHiLo(tgtpdg, &pdghi, &pdglo);

  // Have Target
  if (findPDG == 0){
    MECLoadHadTensor * hadtensor = MECLoadHadTensor::Instance(tgtpdg, nupdg);
    xSecAtE = hadtensor->TotalXsecAtE(tgtpdg, nupdg, Enu);
  }

  // If don't have target, figure out if interpolate, extrapolate, or other
  else {

    // function had error:
    if (findPDG == -1){
      return -1;
    }

    // extrapolate above
    else if (pdghi == -1 && pdglo != -1 ){
      return -1;
    }

    // extrapolate below
    else if (pdghi != -1 && pdglo == -1 ){
      return -1;
    }
  
    // interpolate
    else if (pdghi != -1 && pdglo != -1 ){

      // load larger and smaller target tensors
      MECLoadHadTensor * hadtensorLo = MECLoadHadTensor::Instance(pdglo, nupdg);
      MECLoadHadTensor * hadtensorHi = MECLoadHadTensor::Instance(pdghi, nupdg);
      
      // pull xsecs 
      double xsecLo = hadtensorLo->TotalXsecAtE(pdglo, nupdg, Enu);
      double xsecHi = hadtensorHi->TotalXsecAtE(pdghi, nupdg, Enu);

      // linear interpolate
      xSecAtE = interpolateAllTargets(tgtpdg, pdglo, pdghi, xsecLo, xsecHi);
      
    }// end else if interpolate
  
  }// end else if don't have target
  return xSecAtE;
}
*/

//____________________________________________________________________
// interpolate xsec function
// for now, linear
double MECTensorPXSec::interpolateAllTargets(int tgt, int tgtlo, int tgthi, double xseclo, double xsechi){
  
  int tgta   = GetAFromPDG(tgt);
  int tgtloa = GetAFromPDG(tgtlo);
  int tgthia = GetAFromPDG(tgthi);

  double WeightA = (tgta - tgtloa) / (tgthia - tgtloa);
  double xsec = xseclo + (xsechi - xseclo)*WeightA;

  return xsec; 
}


//_____________________________________________________________________
// loop through return indice of first smallest target
// assumes array is sorted.
int MECTensorPXSec::findPDGHiLo(int tgt, int& tgthi, int& tgtlo){

  // Targets Available:
  int targets[3] = {1000060120, 1000080160, 1000200400};
  int size = 3;
  
  // check if lower than first target:
  if (tgt < targets[0]){
    tgtlo = -1;
    tgthi = targets[0];
    return 1;
  }
  // check if higher than last target:
  if (tgt > targets[size-1]){
    tgtlo = targets[size-1];
    tgthi = -1;
    return 1;
  }
  // check if within target list;
  int index;
  for (int i = 0; i < size; i++){
    // if target is in the list
    if (tgt == targets[i]){
      tgtlo = -1;
      tgthi = -1;
      return 0;
    }
    if (tgt > targets[i]){
      index=i;
    }
    else{
      tgtlo = targets[index];
      tgthi = targets[index+1];
      return 1;
    }
  }
  
  // failed?
  return -1;
}

//_____________________________________________________________
// get A (and Z) from genie 10 digit pdg:
// 10LZZZAAAI where AAA is the total baryon number, ZZZ is the total charge
int MECTensorPXSec::GetAFromPDG(int pdg){

  std::stringstream ss;
  string str;
  string stra;
  int a;

  ss << pdg;
  ss >> str;
  
  for (int i = 3; i <= 5; i++)  stra+= str[i];

  a = atoi(stra.c_str());

  return a;
}

int MECTensorPXSec::GetZFromPDG(int pdg){

  std::stringstream ss;
  string str;
  string strz;
  int z;

  ss << pdg;
  ss >> str;
  
  for (int i = 6; i <= 8; i++)  strz+= str[i];

  z = atoi(strz.c_str());

  return z;
}

void MECTensorPXSec::GetAZFromPDG(int pdg, int& a, int& z){
  std::stringstream ss;
  string str;
  string strz;
  string stra;

  ss << pdg;
  ss >> str;

  for (int i = 3; i <= 5; i++)  stra+= str[i];
  for (int i = 6; i <= 8; i++)  strz+= str[i];

  a = atoi(stra.c_str()); 
  z = atoi(strz.c_str());

  return ;
}
