#include "TRuntimeObjects.h"

#include <iostream>
#include <map>

#include <cstdio>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <sstream>
#include "TRandom.h"

#include "TPreserveGDirectory.h"
#include "TObject.h"
#include "TCaesar.h"
#include "TS800.h"
#include "TS800Hit.h"

#include "TChannel.h"
#include "GValue.h"
#include "TFile.h"
#include "GCutG.h"


void HandleUngated(TRuntimeObjects& obj);
void HandleS800(TRuntimeObjects &obj);
void HandleCaesar(TRuntimeObjects &obj, GCutG *in_cut, GCutG *out_cut);
void HandleCaesarAddback(TRuntimeObjects &obj, std::string dirname, TCaesar *caesar_ab);
void HandleMultiplicityAndCoincidenceMatrices(TRuntimeObjects& obj, std::string dirname, std::vector<double> energy_dc_singles);


void HandleCaesarAddback(TRuntimeObjects &obj, std::string dirname, TCaesar *caesar_ab){
  //Now loop over addback hits
  unsigned int num_addback_hits = caesar_ab->AddbackSize();
//  const int AB_ENERGY_THRESHOLD = 0;
  dirname += "_addback";
  std::string histname;
  std::vector<double> energies_addback_n0n1n2;
  energies_addback_n0n1n2.clear();

  for (unsigned int y=0; y < num_addback_hits; y++){
      TCaesarHit &hit = caesar_ab->GetAddbackHit(y);

      double energy_ab = hit.GetEnergy();
                  
      histname = "energy_addback";
      obj.FillHistogram(dirname,histname,
                        1024,0,8192,energy_ab);
      histname = "detnum_vs_energy_addback";
      obj.FillHistogram(dirname,histname,200,0,200,hit.GetAbsoluteDetectorNumber(),
                                         1125,-100,8192,energy_ab);
              
      if (hit.GetNumHitsContained() == 1 && !hit.is_garbage_addback){
        histname = "energy_addback_n0";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,energy_ab);
        energies_addback_n0n1n2.push_back(energy_ab);
      }
      else if (hit.GetNumHitsContained() == 2 && !hit.is_garbage_addback){
        histname = "energy_addback_n1";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,energy_ab);
        energies_addback_n0n1n2.push_back(energy_ab);

      }
      else if (hit.GetNumHitsContained() == 3 && !hit.is_garbage_addback){
        histname = "energy_addback_n2";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,energy_ab);
        energies_addback_n0n1n2.push_back(energy_ab);
      }
      else if(hit.is_garbage_addback){
        histname = "energy_addback_ng";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,energy_ab);
      }
      else {
        std::cout << "Weird event not meeting any criteria for addback" << std::endl;
        std::cout << "hit.is_garbage_addback    = " << hit.is_garbage_addback << std::endl;
        std::cout << "hit.GetNumHitsContained() = " << hit.GetNumHitsContained() << std::endl;
      }
  }//end for loop over addback hits


  HandleMultiplicityAndCoincidenceMatrices(obj, dirname, energies_addback_n0n1n2);
}

void HandleUngated(TRuntimeObjects& obj) {

  TCaesar  *caesar  = obj.GetDetector<TCaesar>();
  TS800    *s800    = obj.GetDetector<TS800>();

  std::string histname;
  std::string dirname;

  dirname = "Ungated";

  if(s800){
    histname = "S800_E_Up";
    obj.FillHistogram(dirname,histname,4101,-100,4000,s800->GetScint().GetEUp());
    if(s800->GetME1Size()){
      histname = "S800_Time_Up";
      obj.FillHistogram(dirname,histname,2000,-10,10,s800->GetMTof().fE1Up[0]);
    }
    if(s800->GetHodoscope().Size()){
      histname = "Hodoscope";
      obj.FillHistogram(dirname,histname,50,0,50,s800->GetHodoscope().GetHodoHit(0).GetChannel(),
                                       2048,0,8196,s800->GetHodoscope().GetHodoHit(0).GetCharge());
    }
  }

  if(!caesar)
    return;

  histname = "Fera_error_code";
  obj.FillHistogram(dirname,histname,5,0,5,caesar->GetFeraError());

  for(unsigned int y=0;y<caesar->Size();y++) {
      TCaesarHit &hit = caesar->GetCaesarHit(y);

      double charge = hit.GetCharge(); 
      histname = "charge_novalidcond";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,charge);
      double time   = hit.GetTime();
      histname = "time_novalidcond";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,time);

      if(hit.IsOverflow())
        continue;

      if(!hit.IsValid())
        continue;

      double energy = hit.GetEnergy();
      
      if(s800){
        double corr_time = time - ((s800->GetTrigger().GetS800Source())*(0.1/0.25));
        histname = "caesar_energy_vs_time_relative_to_s800";
        obj.FillHistogram(dirname,histname,2001,-2000,2000,corr_time,
                                           2149,-100,8192,energy);
      }
      histname = "energy";
      obj.FillHistogram(dirname,histname,
          1125,-100,8192,energy);
      histname = "charge";
      obj.FillHistogram(dirname,histname,
          8293,-100,8192,charge);
      histname = "time";
      obj.FillHistogram(dirname,histname,
          1001,-2000,2000,time);

      histname = "detnum_vs_energy";
      obj.FillHistogram(dirname,histname,200,0,200,hit.GetAbsoluteDetectorNumber(),
          1125,-100,8192,energy);
      histname = "detnum_vs_charge";
      obj.FillHistogram(dirname,histname,200,0,200,hit.GetAbsoluteDetectorNumber(),
          8293,-100,8192,charge);
      histname = "detnum_vs_time";
      obj.FillHistogram(dirname,histname,200,0,200,hit.GetAbsoluteDetectorNumber(),
          1001,-2000,2000,time);
  }//for loop over singles hits
  return;
}//end HandleUngated

int gates_loaded = 0;
std::vector<GCutG*> incoming_cuts = {};
std::vector<GCutG*> outgoing_cuts = {};
GCutG *tcut = 0;
void LoadGates(TRuntimeObjects &obj){
  TList *gates = &(obj.GetGates());
  TIter iter(gates);
  while(TObject *obj = iter.Next()) {
    GCutG *gate = (GCutG*)obj;
    std::string tag = gate->GetTag();
    if(!tag.compare("incoming")) {
      incoming_cuts.push_back(gate);
    } else if(!tag.compare("outgoing")) {
      outgoing_cuts.push_back(gate);
    } else if(!tag.compare("timeenergy")){
      tcut  = new GCutG(*gate);
    }
    gates_loaded++;
  }
}

void CheckGates(TRuntimeObjects &obj,std::vector<unsigned int> &incoming_passed, std::vector<unsigned int> &outgoing_passed){
  const TS800    *s800    = obj.GetDetector<TS800>();

  double obj_tof = s800->GetMTofXfpE1(); 
  double xfp = s800->GetMTofObjE1(); 
  for (unsigned int i = 0; i < incoming_cuts.size(); i++){
    if (incoming_cuts.at(i)->IsInside(obj_tof, xfp)){
      incoming_passed.push_back(i);
    }
  }

  double obj_corr = s800->GetMTofObjE1();
  double ic_sum = s800->GetIonChamber().GetAve();
  for (unsigned int i = 0; i < outgoing_cuts.size(); i++){
    if (outgoing_cuts.at(i)->IsInside(obj_corr,ic_sum)){
      outgoing_passed.push_back(i);
    }
  }

}

void HandleMultiplicityAndCoincidenceMatrices(TRuntimeObjects& obj, std::string dirname, std::vector<double> energy_dc_singles){

  for (unsigned int i = 0; i < energy_dc_singles.size(); i++){
    for (unsigned int j = i+1; j < energy_dc_singles.size(); j++){
        obj.FillHistogram(dirname,"gamma_gamma_matrix", 2048,0,8192, energy_dc_singles.at(i),
                                                        2048, 0, 8192,energy_dc_singles.at(j));
        obj.FillHistogram(dirname,"gamma_gamma_matrix", 2048,0,8192, energy_dc_singles.at(j),
                                                        2048, 0, 8192,energy_dc_singles.at(i));
    }
  }

  if (energy_dc_singles.size() ==1){
    obj.FillHistogram(dirname,"energy_dc_mult_1", 2048,0,8192, energy_dc_singles.at(0));
  }
  if (energy_dc_singles.size() ==2){
    obj.FillHistogram(dirname,"energy_dc_mult_2", 2048,0,8192, energy_dc_singles.at(0));
    obj.FillHistogram(dirname,"energy_dc_mult_2", 2048,0,8192, energy_dc_singles.at(1));
  }
  if (energy_dc_singles.size() == 3){
    obj.FillHistogram(dirname,"energy_dc_mult_3", 2048,0,8192, energy_dc_singles.at(0));
    obj.FillHistogram(dirname,"energy_dc_mult_3", 2048,0,8192, energy_dc_singles.at(1));
    obj.FillHistogram(dirname,"energy_dc_mult_3", 2048,0,8192, energy_dc_singles.at(2));
  }
  if (energy_dc_singles.size() > 3){
    for (auto energy_dc : energy_dc_singles){
      obj.FillHistogram(dirname, "energy_dc_greater_than_mult_3", 2048,0,8192, energy_dc);
    }
  }

}

void HandleCaesar(TRuntimeObjects &obj, GCutG *in_cut, GCutG *out_cut){
  std::string dirname(Form("caesar_gated_%s_%s", in_cut->GetName(), out_cut->GetName()));

  const TS800    *s800    = obj.GetDetector<TS800>();
  TCaesar    *caesar    = obj.GetDetector<TCaesar>();
  TCaesar  *caesar_ab = new TCaesar();//will be used to do addback correction 
                                      //only for the hits inside the time cut
  if(!caesar || !s800){
    return;
  }
  std::vector<double> energy_dc_singles = {};
  for (unsigned int i = 0; i < caesar->Size(); i++){
    TCaesarHit &hit = caesar->GetCaesarHit(i);
    double energy = hit.GetEnergy();
    double energy_dc = hit.GetDoppler();
    double time = hit.GetTime();
    int det_num = hit.GetAbsoluteDetectorNumber();
    double corr_time = time - ((s800->GetTrigger().GetS800Source())*(0.1/0.25));
    obj.FillHistogram(dirname,"caesar_energy_vs_time_relative_to_s800",
        2001,-2000,2000,corr_time,
        2149,-100,8192,energy);


    if (tcut){
      if (tcut->IsInside(corr_time, energy)){
        caesar_ab->InsertHit(hit);
        obj.FillHistogram(dirname,"energy_summary", 192,0,192, det_num,
            2149,-100,8192,energy);
        obj.FillHistogram(dirname,"time_summary", 192,0,192, det_num,
            2149,-2000,2000,time);
        obj.FillHistogram(dirname,"energy_dc_summary",192,0,192, det_num,
            2149,-100,8192,energy_dc);
        energy_dc_singles.push_back(energy_dc);
      }
    }
  }

  HandleMultiplicityAndCoincidenceMatrices(obj, dirname,energy_dc_singles);
  HandleCaesarAddback(obj, dirname,caesar_ab);
  return;
}

void HandleS800(TRuntimeObjects &obj){
  TS800    *s800    = obj.GetDetector<TS800>();
  
  if (s800){
    double obj_tof = s800->GetMTof().GetCorrelatedObjE1(); 
    double xfp_tof = s800->GetMTof().GetCorrelatedXfpE1(); 
    obj.FillHistogram("s800", "incoming_pid",  
        1000, -3000, 1000, obj_tof,
        1000, 1000, 3000, xfp_tof);

//    double obj_corr = s800->GetCorrTOF_OBJ_MESY();
    double ic_sum = s800->GetIonChamber().GetAve();
    obj.FillHistogram("s800", "uncorrected_outgoing_pid", 
        1000, -3000, 1000, obj_tof,
        4096, 0, 4096, ic_sum);

    double obj_corr = s800->GetMTofObjE1();
    obj.FillHistogram("s800", "corrected_outgoing_pid", 
        1000, -3000, 1000, obj_corr,
        4096, 0, 4096, ic_sum);

    //Now XFP/AFP corrections
//  std::vector<int> afp_corrs =  {1900, 1950, 2000, 2050, 2100, 2200};
//  std::vector<double> xfp_corrs =  {0.17, 0.175, 0.18, 0.2};
//  double afp = s800->GetAFP();//angle in focal plane
//  double xfp = s800->GetXFP();//x position in focal plane
//  for (auto afp_corr : afp_corrs){
//    for (auto xfp_corr : xfp_corrs){
//      double test_obj_corr = s800->GetTofE1_MTDC(afp_corr, xfp_corr);
//      obj.FillHistogram("s800_afpcorrs", Form("afp_vs_obj_afp_%d_xfp_%lf", afp_corr, xfp_corr), 
//          2000, -4000, 0, test_obj_corr,
//          1000, -0.1, 0.1, afp);
//      obj.FillHistogram("s800_xfpcorrs", Form("xfp_vs_obj_afp_%d_xfp_%lf", afp_corr, xfp_corr), 
//          2000, -4000, 0, test_obj_corr,
//          600, -300, 300, xfp);
//    }
//  }
  }
}
// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {

  TList *list  = &(obj.GetObjects());
  TList *gates = &(obj.GetGates());
  int numobj = list->GetSize();

  if (gates_loaded != gates->GetSize()){
    LoadGates(obj);
  }

  std::vector<unsigned int> incoming_passed_cuts;
  std::vector<unsigned int> outgoing_passed_cuts;
  if (incoming_cuts.size() || outgoing_cuts.size()){
    CheckGates(obj, incoming_passed_cuts, outgoing_passed_cuts);
  }

  HandleS800(obj);
//HandleUngated(obj);
//if(incoming_passed_cuts.size() && outgoing_passed_cuts.size()) {
//  for (auto incoming_passed : incoming_passed_cuts){
//    for (auto outgoing_passed : outgoing_passed_cuts){
//      HandleCaesar(obj,incoming_cuts.at(incoming_passed),outgoing_cuts.at(outgoing_passed));
//    }
//  }
//}
  

  if(numobj!=list->GetSize())
    list->Sort();
}
