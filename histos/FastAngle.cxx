#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
//
// Note: from TRuntimeObjects:
//    TList& GetDetectors();
//    TList& GetObjects();
//

#include <cstdio>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include "TRandom.h"

#include "TObject.h"
#include "TFastScint.h"

float gamma1 = 1332;
float gamma2 = 1173;
//float gamma1 = 1836;
//float gamma2 = 898;
float tol1 = 30;
float tol2 = 30;
unsigned int   eventnum = 0;
bool first_time = true;

unsigned int depth = 2;
std::vector<TFastScint> lastfasts;

extern "C"

void MakeHistograms(TRuntimeObjects& obj) {
  TFastScint *fast = obj.GetDetector<TFastScint>();

  if(!fast)
    return;

  TList *list = &(obj.GetObjects());
  int numobj = list->GetSize(); 

  std::string dirname;
  std::string histname;

  dirname = "coincidence_energy_gated";

  //if(first_time){
    //for(int a=0;a<16;a++){    
      //if(a==0 || a==6 || a==12) continue;
      //for(int b=0;b<16;b++){
        //if(b==0 || b==6 || b==12) continue;
        //if(a==b) continue;
        //obj.FillHistogram(dirname,"Detector_Pairs_Per_Angle",
                          //180,0,180,fast->GetPosition(a).Angle(fast->GetPosition(b))*TMath::RadToDeg());
      //}
    //}
    //first_time = false; 
  //}//end if first_time through code, make normalization for number of counts based on detector pair angles

  for(unsigned int i=0;i<fast->Size();i++) {
    TFastScintHit hit = fast->GetLaBrHit(i);

    histname = "ChannelCharge";
    obj.FillHistogram(histname,4000,0,4000,hit.Charge(),
                             20,0,20,hit.GetChannel());

    histname = "ChannelEnergy";
    obj.FillHistogram(histname,4000,0,4000,hit.GetEnergy(),
                             20,0,20,hit.GetChannel());

    if(abs(hit.GetEnergy() - gamma1) < tol1){

      if(eventnum>depth) {
        for(unsigned int x=0;x<lastfasts.size()-1;x++) {
          TFastScint &last = lastfasts.at(x);
          for(unsigned int y=0;y<last.Size();y++) {
            if(last.GetLaBrHit(y).GetChannel()==hit.GetChannel()) 
            continue;
            if(last.GetLaBrHit(y).Charge()>3500) continue;
            double angle = hit.GetPosition().Angle(last.GetLaBrHit(y).GetPosition());
            histname = Form("energy_v_angle_gated_on_%5.1f_isotropic",gamma1);
            obj.FillHistogram(dirname,histname,180,0,180,angle*TMath::RadToDeg(),
                                               4000,0,4000,last.GetLaBrHit(y).GetEnergy());
            histname = Form("energy_v_angle_gated_on_%5.1f_in_channel_%i_isotropic",gamma1,hit.GetChannel());
            obj.FillHistogram(dirname,histname,180,0,180,angle*TMath::RadToDeg(),
                                               4000,0,4000,last.GetLaBrHit(y).GetEnergy());

            if((hit.GetChannel() == 1) || (hit.GetChannel() == 3) || (hit.GetChannel() == 7) || (hit.GetChannel() == 8) ||        (hit.GetChannel() == 9)){
               histname = Form("energy_v_angle_gated_on_%5.1f_in_channel_1,3,7,8,9_isotropic",gamma1);
               obj.FillHistogram(dirname,histname,180,0,180,angle*TMath::RadToDeg(),
                                                  4000,0,4000,last.GetLaBrHit(y).GetEnergy());
            }

            if((hit.GetChannel() == 1) || (hit.GetChannel() == 3) || (hit.GetChannel() == 5) || (hit.GetChannel() == 7) || (hit.GetChannel() == 8) || (hit.GetChannel() == 9) || (hit.GetChannel() == 10) || (hit.GetChannel() == 11) || (hit.GetChannel() == 13) || (hit.GetChannel() == 14)){
               histname = Form("energy_v_angle_gated_on_%5.1f_in_channel_1,3,5,7,8,9,10,11,13,14_isotropic",gamma1);
               obj.FillHistogram(dirname,histname,180,0,180,angle*TMath::RadToDeg(),
                                                  4000,0,4000,last.GetLaBrHit(y).GetEnergy());
            }
              
            }//end for loop over y
          }//end for loop over x 
        }//end if eventnum>depth

      for(unsigned int j=0;j<fast->Size();j++) {
        if(i==j)
          continue;

        //if((abs(hit.GetEnergy() - gamma1) < tol1) && (abs(hit2.GetEnergy() - gamma1) < tol1) && i>j){
          //continue;
        //}//if 1332 and 1332 in coincidence, only fill once      

        TFastScintHit hit2 = fast->GetLaBrHit(j);

        if(hit2.Charge()>3500) continue;
        
        histname = Form("energy_v_angle_gated_on_%5.1f",gamma1);
        obj.FillHistogram(dirname,histname,180,0,180,hit.GetPosition().Angle(hit2.GetPosition())*TMath::RadToDeg(),
                                           4000,0,4000,hit2.GetEnergy());

        if(abs(hit.GetEnergy() - gamma2) < tol2){
          histname = Form("%5.1f_and_%5.1f_coincidence_gamma_count_in_dets",gamma1,gamma2);
          obj.FillHistogram(dirname,histname,20,0,20,hit.GetChannel());
          obj.FillHistogram(dirname,histname,20,0,20,hit2.GetChannel());
        }

        histname = Form("energy_v_angle_gated_on_%5.1f_in_channel_%i",gamma1,hit.GetChannel());
        obj.FillHistogram(dirname,histname,180,0,180,hit.GetPosition().Angle(hit2.GetPosition())*TMath::RadToDeg(),
                                           4000,0,4000,hit2.GetEnergy());

        if((hit.GetChannel() == 1) || (hit.GetChannel() == 3) || (hit.GetChannel() == 7) || (hit.GetChannel() == 8) ||        (hit.GetChannel() == 9)){
               histname = Form("energy_v_angle_gated_on_%5.1f_in_channel_1,3,7,8,9",gamma1);
               obj.FillHistogram(dirname,histname,180,0,180,hit.GetPosition().Angle(hit2.GetPosition())*TMath::RadToDeg(),
                                                  4000,0,4000,hit2.GetEnergy());
        }

        if((hit.GetChannel() == 1) || (hit.GetChannel() == 3) || (hit.GetChannel() == 5) || (hit.GetChannel() == 7) || (hit.GetChannel() == 8) || (hit.GetChannel() == 9) || (hit.GetChannel() == 10) || (hit.GetChannel() == 11) || (hit.GetChannel() == 13) || (hit.GetChannel() == 14)){
               histname = Form("energy_v_angle_gated_on_%5.1f_in_channel_1,3,5,7,8,9,10,11,13,14",gamma1);
               obj.FillHistogram(dirname,histname,180,0,180,hit.GetPosition().Angle(hit2.GetPosition())*TMath::RadToDeg(),
                                                  4000,0,4000,hit2.GetEnergy());
        }

      }//end for loop over j
    }//end for loop over i
  }//end if gamma within tolerance of gamma1

  if(lastfasts.size()>depth) {
    lastfasts.erase(lastfasts.begin());
  }
  lastfasts.push_back(*fast);

  eventnum++;

  if(numobj!=list->GetSize())
    list->Sort();

}//end MakeHistograms
