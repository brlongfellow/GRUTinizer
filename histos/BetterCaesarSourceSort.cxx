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

int HandleUngated(TRuntimeObjects& obj) {

  TCaesar  *caesar  = obj.GetDetector<TCaesar>();
  TCaesar  *caesar_ab = new TCaesar();//will be used to do addback correction 
                                      //only for the hits inside the time cut
  TS800    *s800    = obj.GetDetector<TS800>();

  std::string histname;
  std::string dirname;

  dirname = "Ungated";

  if(s800){
    histname = "S800_E_Up";
    obj.FillHistogram(dirname,histname,4101,-100,4000,s800->GetScint().GetEUp());
    histname = "S800_E_Down";
    obj.FillHistogram(dirname,histname,4101,-100,4000,s800->GetScint().GetEDown());
    if(s800->GetME1Size()){
      histname = "S800_Time_Up";
      obj.FillHistogram(dirname,histname,2000,-10,10,s800->GetMTof().fE1Up[0]);
 //   histname = "S800_Time_Down";
 //   obj.FillHistogram(dirname,histname,2000,-10,10,s800->GetMTof().fE1Down[0]);
    }
    //if(s800->GetHodoscope().Size()){
      //histname = "Hodoscope";
      //obj.FillHistogram(dirname,histname,50,0,50,s800->GetHodoscope().GetHodoHit(0).GetChannel(),
      //                                 2048,0,8196,s800->GetHodoscope().GetHodoHit(0).GetCharge());
    //}
    for(unsigned int i=0; i<s800->GetHodoscope().Size(); i++){
      histname = "Hodoscope_Charge";
      obj.FillHistogram(dirname,histname,50,0,50,s800->GetHodoscope().GetHodoHit(i).GetChannel(),
                                       2048,0,8192,s800->GetHodoscope().GetHodoHit(i).GetCharge());
      histname = "Hodoscope_Energy";
      obj.FillHistogram(dirname,histname,50,0,50,s800->GetHodoscope().GetHodoHit(i).GetChannel(),
                                       2048,0,8192,s800->GetHodoscope().GetHodoHit(i).GetEnergy());
      histname = "Hodoscope_Time";
    obj.FillHistogram(dirname,histname,50,0,50,s800->GetHodoscope().GetHodoHit(i).GetChannel(),
                                       1000,-100000,100000,s800->GetHodoscope().GetHodoHit(i).GetTime());

    histname = "Hodoscope_TimeEnergy";
    obj.FillHistogram(dirname,histname,1024,0,8192,s800->GetHodoscope().GetHodoHit(0).GetEnergy(),
                                       1000,-100000,100000,s800->GetHodoscope().GetHodoHit(0).GetTime());
    }
  }

  if(!caesar)
    return false;

  histname = "Fera_error_code";
  obj.FillHistogram(dirname,histname,5,0,5,caesar->GetFeraError());

  //const int SINGLES_ENERGY_THRESHOLD = 150;
  const int AB_ENERGY_THRESHOLD = 0;

  //std::vector<double> energies_singles;
  std::vector<double> energies_addback;
  std::vector<double> energies_addback_n0;
  std::vector<double> energies_addback_n1;
  std::vector<double> energies_addback_n2;
  std::vector<double> energies_addback_ng;
  //energies_singles.clear();
  energies_addback.clear();

  for(unsigned int y=0;y<caesar->Size();y++) {
      TCaesarHit &hit = caesar->GetCaesarHit(y);

      double charge = hit.GetCharge(); 

      histname = "charge_novalidcondapplied";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,charge);

      if(hit.IsOverflow())
        continue;

      if(!hit.IsValid())
        continue;

      double energy = hit.GetEnergy();
      
      double time   = hit.GetTime();
     
      if(s800){
        double corr_time = time - ((s800->GetTrigger().GetS800Source())*(0.1/0.25));
        dirname = "Ungated";
        histname = "Caesar Energy vs Caesar Time Relative to LaBr";
        obj.FillHistogram(dirname,histname,2001,-2000,2000,corr_time,
                                           2149,-100,8192,energy);
        histname = "S800 Energy vs Caesar Energy";
        obj.FillHistogram(dirname,histname,2149,-100,8192,s800->GetScint().GetEUp(),
                                           2149,-100,8192,energy);
      }

      caesar_ab->InsertHit(hit);
        //if (energy > SINGLES_ENERGY_THRESHOLD){              
        //       energies_singles.push_back(energy);
        //}

        dirname = "Ungated";
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

  //Now loop over addback hits
  int num_addback_hits = caesar_ab->AddbackSize();
  for (int y=0; y < num_addback_hits; y++){
      TCaesarHit &hit = caesar_ab->GetAddbackHit(y);

      double energy_ab = hit.GetEnergy();
                  
      if (energy_ab > AB_ENERGY_THRESHOLD){    
        energies_addback.push_back(energy_ab);       
      }//For multiplicity purposes

      dirname = "Ungated";
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
        energies_addback_n0.push_back(energy_ab);

      }
      else if (hit.GetNumHitsContained() == 2 && !hit.is_garbage_addback){
        histname = "energy_addback_n1";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,energy_ab);
        energies_addback_n1.push_back(energy_ab);

      }
      else if (hit.GetNumHitsContained() == 3 && !hit.is_garbage_addback){
        histname = "energy_addback_n2";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,energy_ab);
        energies_addback_n2.push_back(energy_ab);
      }
      else if(hit.is_garbage_addback){
        histname = "energy_addback_ng";
        obj.FillHistogram(dirname,histname,
                          1125,-100,8192,energy_ab);
        energies_addback_ng.push_back(energy_ab);

      }
      else {
        std::cout << "Weird event not meeting any criteria for addback" << std::endl;
        std::cout << "hit.is_garbage_addback    = " << hit.is_garbage_addback << std::endl;
        std::cout << "hit.GetNumHitsContained() = " << hit.GetNumHitsContained() << std::endl;
      }
      
  }//end for loop over addback hits

  return 0;

}//end HandleUngated


// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {

  TList *list  = &(obj.GetObjects());
//  TList *gates = &(obj.GetGates());
  int numobj = list->GetSize();

  HandleUngated(obj);

//  if(gates_loaded!=gates->GetSize()) {
//    TIter iter(gates);
//    while(TObject *obj = iter.Next()) {
//      GCutG *gate = (GCutG*)obj;
//      std::string tag = gate->GetTag();
//      if(!tag.compare("incoming")) {
//        incoming_cuts.push_back(gate);
//      } else if(!tag.compare("outgoing")) {
//        outgoing_cuts.push_back(gate);
//      } else if(!tag.compare("timeenergy")){
//        timeenergy_cuts.push_back(gate);
//      }
//      gates_loaded++;
//    }
//  }
//
//  //printf("incoming.size() == %i\n",incoming_cuts.size());
//  int incoming_passed=-1;
//  int outgoing_passed=-1;
//  for(unsigned int x=0;x<incoming_cuts.size();x++) {
//    bool passed = OutgoingBeam(obj,incoming_cuts.at(x)); 
//    if(x!=0 && passed) {
//      incoming_passed = x;
//      break;
//    }
//  }
//  for(unsigned int x=0;x<outgoing_cuts.size();x++) {
//    bool passed = IncomingBeam(obj,outgoing_cuts.at(x)); 
//    if(x!=0 && passed) {
//      outgoing_passed = x;
//      break;
//    }
//  } 
//
//  if(incoming_passed>0 && outgoing_passed>0) {
//    HandleCaesar(obj,incoming_cuts.at(incoming_passed),outgoing_cuts.at(outgoing_passed));
//  }
//
  if(numobj!=list->GetSize())
    list->Sort();
}
