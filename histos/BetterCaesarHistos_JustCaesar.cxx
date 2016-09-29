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

#include "TChannel.h"
#include "GValue.h"
#include "TFile.h"
#include "TCutG.h"


const int total_det_in_prev_rings[N_RINGS] = {0,10,24,48,72,96,120,144,168,182};
const double START_ANGLE = 3.2;
const double FINAL_ANGLE = 3.2;
const double ANGLE_STEPS = 0.1;
const int TOTAL_ANGLES = (FINAL_ANGLE-START_ANGLE)/ANGLE_STEPS + 1;

std::vector<double> angles;


bool init_called = false;

void Init() {
  double temp_angle = START_ANGLE;
  for (int i = 0; i < TOTAL_ANGLES; i++){
    angles.push_back(temp_angle);
    temp_angle += ANGLE_STEPS;
  }
}


// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  TCaesar  *caesar  = obj.GetDetector<TCaesar>();
  TS800    *s800    = obj.GetDetector<TS800>();

  TList *list = &(obj.GetObjects());
  int numobj = list->GetSize();
  if(caesar) {
    if(!init_called)
      Init();

    int valid_counter=0;

    for(unsigned int i=0;i<caesar->Size();i++) {
      TCaesarHit hit = caesar->GetCaesarHit(i);
      if(!hit.IsValid())
        continue;
      if(hit.IsOverflow())
        continue;
  
      valid_counter++;

      int ring = hit.GetRingNumber();

      obj.FillHistogram("s800","TrigBit",30,0,30,s800->GetTrigger().GetRegistr());
      
      std::string histname = Form("CASESAR_raw_time_TrigBit_%i",s800->GetTrigger().GetRegistr());

      obj.FillHistogram("Caesar",histname,4096,0,4096,hit.Time());
      //obj.FillHistogram("Caesar","GetTime_vs_GetDoppler",2000,0,2000,hit.GetTime(),
                                                         //2048,0,8192,hit.GetDoppler(beta,z_shift));
      //if(hit.GetDoppler(beta,z_shift)>300){
        //obj.FillHistogram("Caesar","GetTime_vs_GetAbsoluteDetectorNumber",200,0,200,hit.GetAbsoluteDetectorNumber(),
                                                                          //2000,0,2000,hit.GetTime());
      //}//end if Doppler > 300

      histname = Form("Detector_Charge_Summary_TrigBit_%i",s800->GetTrigger().GetRegistr());

      obj.FillHistogram("Caesar",histname,200,0,200,hit.GetDetectorNumber()+total_det_in_prev_rings[ring],
                                                  4096,0,4096,hit.Charge());

      histname = Form("Detector_Energy_Summary_TrigBit_%i",s800->GetTrigger().GetRegistr());

      obj.FillHistogram("Caesar",histname,200,0,200,hit.GetDetectorNumber()+total_det_in_prev_rings[ring],
                                                  2048,0,8192,hit.GetEnergy());

      //obj.FillHistogram("Caesar","Detector_Doppler_Summary",300,0,300,hit.GetDetectorNumber()+total_det_in_prev_rings[ring],
                                                  //1024,0,8192,hit.GetDoppler(beta,z_shift));

      //for(unsigned int j=0;j<caesar->Size();j++) {
        //if(i==j)
          //continue;
        //TCaesarHit hit2 = caesar->GetCaesarHit(j);
      //if(!hit2.IsValid())
        //continue;
      //if(hit2.IsOverflow())
        //continue;
        //obj.FillHistogram("Caesar","Energy_Coincidence_Matrix",2048,0,8192,hit.GetEnergy(),
                                                             //2048,0,8192,hit2.GetEnergy());

      //}//end for loop over coincidence hits

     

      /**************************************\
      | Lets think about the location of this|
      | if(s800) that is inside our loop over|
      | CAESAR Hits!!!!                      |
      \**************************************/


    }// end of caesar for hit loop
    obj.FillHistogram("Caesar","Multiplicity",300,0,300,valid_counter);
  
  } // I have a valid caesar event

  if(numobj!=list->GetSize())
    list->Sort();

}//end MakeHistograms
