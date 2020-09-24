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
#include "GCutG.h"


std::vector<GCutG*> incoming_cuts   = {0};
std::vector<GCutG*> outgoing_cuts   = {0};
std::vector<GCutG*> timeenergy_cuts; // = {0};
int gates_loaded=0;
long timestamp_old=0;
std::vector<long> timestamps = {0};

bool OutgoingBeam(TRuntimeObjects& obj,GCutG *incoming) {
  TS800    *s800    = obj.GetDetector<TS800>();

  if(!s800)
    return false;

  std::string histname;
  std::string dirname;
  std::string name;

  if(incoming){
    dirname = Form("outgoing_%s",incoming->GetName());
    name = incoming->GetName();
  }
  else{
    dirname = "outgoing";
  }


  //double objtac = s800->GetTof().GetTacOBJ();
  //double objmesy = s800->GetMTof().GetCorrelatedObjE1();
  double objmesy = s800->GetMTofObjE1();
  //double xfptac = s800->GetTof().GetTacXFP();
  //double xfpmesy = s800->GetMTof().GetCorrelatedXfpE1();
  double xfpmesy = s800->GetMTofXfpE1();
  long timestamp = s800->GetTimestamp();
  long counter = s800->GetEventCounter();
  timestamps.push_back(timestamp);

  //if(incoming) {
    //if(!incoming->IsInside(objmesy,xfpmesy))
      //return false;
  //}

  //ignore incoming gate condition for s44 since xfp tripped off for half the runs
  //x.compare(y) returns 0 if strings are same, <0 if x<y, >0 if x>y
  if(incoming&&name.compare("s44_incoming")) {
    if(!incoming->IsInside(objmesy,xfpmesy))
      return false;
  }

  double ic_sum         = s800->GetIonChamber().GetAve();
  //double objtac_corr    = s800->GetCorrTOF_OBJTAC();
  double afp            = s800->GetAFP();
  double bfp            = s800->GetBFP();
  double xfp_focalplane = s800->GetXFP(0);


  //std::cout << "---------------------------------------------------------" << std::endl;
  //std::cout << "timestamp_old = " << timestamps.at(timestamps.size()-2) << std::endl;
  //std::cout << "timestamp = " << timestamps.at(timestamps.size()-1) << std::endl;
  //std::cout << "timestamp_actual = " << timestamp << std::endl;
  //std::cout << "xfpmesy = " << xfpmesy << std::endl;

  //if(timestamp == timestamp_old)
    //std::cout << "TS are exactly the same" << std::endl;

  //std::cout << "---------------------------------------------------------" << std::endl;


  //if(bfp<-0.025) {
    //return false;
  //}


  histname = "ts_vs_eventcounter";
  obj.FillHistogram(dirname,histname,
                    1000,0,1000000000,timestamp-230768802,
        	    1000,0,100000,counter);

  obj.FillHistogram(dirname,"ts_diff",1000,-1000,100000,timestamp-timestamp_old);

  //std::cout << "ts_diff = " << timestamp-timestamp_old << std::endl;
  
  timestamp_old = timestamp;

  obj.FillHistogram(dirname,"E1_up_size",10,0,10,s800->GetME1Size());

  for(int i=0; i<s800->GetME1Size(); i++){
    obj.FillHistogram(dirname,"E1_up",1000,-1000,100000,s800->GetME1Up(i));
    if(i>0){
      obj.FillHistogram(dirname,"E1_up_mult>1",1000,-1000,100000,s800->GetME1Up(i));
    }
  }


  histname = "AFP_vs_OBJTOF";
  obj.FillHistogram(dirname,histname,
  		    1000,-3000,-1000,objmesy,
        	    1000,-1,1,afp); // check units of AFP

  histname = "BFP_vs_OBJTOF";
  obj.FillHistogram(dirname,histname,
  		    1000,-3000,-1000,objmesy,
        	    1000,-1,1,bfp); // check units of BFP
 
  histname = "XFP_vs_OBJTOF";
  obj.FillHistogram(dirname,histname,
        	    1000,-3000,-1000,objmesy,
        	    1000,-300,300,xfp_focalplane);

  histname = "IC_vs_OBJTOF_PID";
  obj.FillHistogram(dirname,histname,
        	    1000,-3000,-1000,objmesy,
        	    1000,-100,4000,ic_sum);

  histname = "AFP_vs_XFPTOF";
  obj.FillHistogram(dirname,histname,
  		    1000,1000,4000,xfpmesy,
        	    1000,-1,1,afp); // check units of AFP

  histname = "BFP_vs_XFPTOF";
  obj.FillHistogram(dirname,histname,
  		    1000,1000,4000,xfpmesy,
        	    1000,-1,1,bfp); // check units of BFP
 
  histname = "XFP_vs_XFPTOF";
  obj.FillHistogram(dirname,histname,
        	    1000,1000,4000,xfpmesy,
        	    1000,-300,300,xfp_focalplane);

  histname = "IC_vs_XFPTOF_PID";
  obj.FillHistogram(dirname,histname,
        	    1000,1000,4000,xfpmesy,
        	    1000,-100,4000,ic_sum);

  histname = "IC_vs_CRDC1Y";
  obj.FillHistogram(dirname,histname,
        	    500,-500,500,s800->GetCrdc(0).GetNonDispersiveY(),
        	    1000,-100,4000,ic_sum);

  histname = "IC_vs_CRDC2Y";
  obj.FillHistogram(dirname,histname,
        	    500,-500,500,s800->GetCrdc(1).GetNonDispersiveY(),
        	    1000,-100,4000,ic_sum);



  obj.FillHistogram(dirname,"CRDC1Y",500,-500,500,s800->GetCrdc(0).GetNonDispersiveY());
  obj.FillHistogram(dirname,"CRDC2Y",500,-500,500,s800->GetCrdc(1).GetNonDispersiveY());
  obj.FillHistogram(dirname,"CRDC1X",400,-400,400,s800->GetCrdc(0).GetDispersiveX());
  obj.FillHistogram(dirname,"CRDC2X",400,-400,400,s800->GetCrdc(1).GetDispersiveX());
  obj.FillHistogram(dirname,"TrigBit",30,0,30,s800->GetTrigger().GetRegistr());
  obj.FillHistogram(dirname,"S800_YTA",1000,-50,50,s800->GetYta());
  obj.FillHistogram(dirname,"S800_DTA",1000,-0.2,0.2,s800->GetDta());
  obj.FillHistogram(dirname,"ATA_vs_BTA",1000,-0.2,0.2,s800->GetAta(),
			                 1000,-0.2,0.2,s800->GetBta());
  obj.FillHistogram(dirname,"CRDC1_MaxPad",300,0,300,s800->GetCrdc(0).GetMaxPad(),
                                           2000,0,8000, s800->GetCrdc(0).GetMaxPadSum());
  obj.FillHistogram(dirname,"CRDC2_MaxPad",300,0,300,s800->GetCrdc(1).GetMaxPad(),
                                           2000,0,8000, s800->GetCrdc(1).GetMaxPadSum());

  histname = "AFP_vs_OBJTOF-XFPTOF";
  obj.FillHistogram(dirname,histname,
  		    1000,-6000,0,objmesy-xfpmesy,
        	    1000,-1,1,afp); // check units of AFP

  histname = "IC_sum";
  obj.FillHistogram(dirname,histname,1000,-100,4000,ic_sum);


  for(unsigned int i=0; i<s800->GetHodoscope().Size(); i++){
    histname = "Hodoscope";
    obj.FillHistogram(dirname,histname,50,0,50,s800->GetHodoscope().GetHodoHit(i).GetChannel(),
                                       1024,0,8192,s800->GetHodoscope().GetHodoHit(i).GetEnergy());
  }
  

  if(s800->GetTrigger().GetRegistr()==1){
    obj.FillHistogram(dirname,"S800_YTA_tb1",1000,-50,50,s800->GetYta());
    obj.FillHistogram(dirname,"S800_DTA_tb1",1000,-0.2,0.2,s800->GetDta());
    obj.FillHistogram(dirname,"ATA_vs_BTA_tb1",1000,-0.2,0.2,s800->GetAta(),
			                       1000,-0.2,0.2,s800->GetBta());
    obj.FillHistogram(dirname,"E1_up_tb1",1000,-1000,100000,s800->GetME1Up(0));
  }

  if(s800->GetTrigger().GetRegistr()==2){
    obj.FillHistogram(dirname,"S800_YTA_tb2",1000,-50,50,s800->GetYta());
    obj.FillHistogram(dirname,"S800_DTA_tb2",1000,-0.2,0.2,s800->GetDta());
    obj.FillHistogram(dirname,"ATA_vs_BTA_tb2",1000,-0.2,0.2,s800->GetAta(),
			                       1000,-0.2,0.2,s800->GetBta());
    obj.FillHistogram(dirname,"E1_up_tb2",1000,-1000,100000,s800->GetME1Up(0));
  }

  if(s800->GetTrigger().GetRegistr()==3){
    obj.FillHistogram(dirname,"S800_YTA_tb3",1000,-50,50,s800->GetYta());
    obj.FillHistogram(dirname,"S800_DTA_tb3",1000,-0.2,0.2,s800->GetDta());
    obj.FillHistogram(dirname,"ATA_vs_BTA_tb3",1000,-0.2,0.2,s800->GetAta(),
			                       1000,-0.2,0.2,s800->GetBta());
    obj.FillHistogram(dirname,"E1_up_tb3",1000,-1000,100000,s800->GetME1Up(0));
  }

  return true;
}

bool IncomingBeam(TRuntimeObjects& obj,GCutG *outgoing) {
  TS800    *s800    = obj.GetDetector<TS800>();
  if(!s800)
    return false;

  std::string histname;
  std::string dirname;
  if(outgoing)
    dirname = Form("incoming_%s",outgoing->GetName());
  else
    dirname = "incoming";
  
  double ic_sum         = s800->GetIonChamber().GetAve();
  //double objtac_corr    = s800->GetCorrTOF_OBJTAC();
  double objmesy_noncorr = s800->GetMTof().GetCorrelatedObjE1();
  double objmesy = s800->GetMTofObjE1();
  if(outgoing) {
    if(!outgoing->IsInside(objmesy,ic_sum))
      return false;
  }

  //double objtac = s800->GetTof().GetTacOBJ();
  //double xfptac = s800->GetTof().GetTacXFP();
  double xfpmesy_noncorr = s800->GetMTof().GetCorrelatedXfpE1();
  double xfpmesy = s800->GetMTofXfpE1();
  histname = "IncomingPID";
  obj.FillHistogram(dirname,histname,
                     1000,-3000,-1000,objmesy,
                     1000,1000,4000,xfpmesy);
  histname = "IncomingPID_noncorr";
  obj.FillHistogram(dirname,histname,
                     1000,-3000,-1000,objmesy_noncorr,
                     1000,1000,4000,xfpmesy_noncorr);
  return true;
}

//int HandleUngated(TRuntimeObjects& obj) {
//
//  TCaesar  *caesar  = obj.GetDetector<TCaesar>();
//  TCaesar  *caesar_ab = new TCaesar();//will be used to do addback correction 
//                                      //only for the hits inside the time cut
//
//  if(!caesar)
//    return false;
//
//  std::string histname;
//  std::string dirname;
//
//  const int SINGLES_ENERGY_THRESHOLD = 150;
//  const int AB_ENERGY_THRESHOLD = 0;
//
//  std::vector<double> energies_singles;
//  std::vector<double> energies_addback;
//  std::vector<double> energies_addback_n0;
//  std::vector<double> energies_addback_n1;
//  std::vector<double> energies_addback_n2;
//  std::vector<double> energies_addback_ng;
//  energies_singles.clear();
//  energies_addback.clear();
//
//  for(unsigned int y=0;y<caesar->Size();y++) {
//      TCaesarHit &hit = caesar->GetCaesarHit(y);
//
//      if(hit.IsOverflow())
//        continue;
//
//      double energy = hit.GetEnergy(); 
//
//      caesar_ab->InsertHit(hit);
//        if (energy > SINGLES_ENERGY_THRESHOLD){              
//               energies_singles.push_back(energy);
//        }
//
//        dirname = "Ungated";
//        histname = "energy";
//        obj.FillHistogram(dirname,histname,
//                          1024,0,8192,energy);
//         
//        histname = "detnum_vs_energy";
//        obj.FillHistogram(dirname,histname,200,0,200,hit.GetAbsoluteDetectorNumber(),
//                                           1024,0,8192,energy);
//
//  }//for loop over singles hits
//
//  //Now loop over addback hits
//  int num_addback_hits = caesar_ab->AddbackSize();
//  for (int y=0; y < num_addback_hits; y++){
//      TCaesarHit &hit = caesar_ab->GetAddbackHit(y);
//
//      double energy = hit.GetEnergy();
//                  
//      if (energy > AB_ENERGY_THRESHOLD){    
//        energies_addback.push_back(energy);       
//      }//For multiplicity purposes
//
//      dirname = "Ungated";
//      histname = "energy_addback";
//      obj.FillHistogram(dirname,histname,
//                        1024,0,8192,energy);
//              
//      if (hit.GetNumHitsContained() == 1 && !hit.is_garbage_addback){
//        histname = "energy_addback_n0";
//        obj.FillHistogram(dirname,histname,
//                          1024,0,8192,energy);
//        energies_addback_n0.push_back(energy);
//
//      }
//      else if (hit.GetNumHitsContained() == 2 && !hit.is_garbage_addback){
//        histname = "energy_addback_n1";
//        obj.FillHistogram(dirname,histname,
//                          1024,0,8192,energy);
//        energies_addback_n1.push_back(energy);
//
//      }
//      else if (hit.GetNumHitsContained() == 3 && !hit.is_garbage_addback){
//        histname = "energy_addback_n2";
//        obj.FillHistogram(dirname,histname,
//                          1024,0,8192,energy);
//        energies_addback_n2.push_back(energy);
//      }
//      else if(hit.is_garbage_addback){
//        histname = "energy_addback_ng";
//        obj.FillHistogram(dirname,histname,
//                          1024,0,8192,energy);
//        energies_addback_ng.push_back(energy);
//
//      }
//      else {
//        std::cout << "Weird event not meeting any criteria for addback" << std::endl;
//        std::cout << "hit.is_garbage_addback    = " << hit.is_garbage_addback << std::endl;
//        std::cout << "hit.GetNumHitsContained() = " << hit.GetNumHitsContained() << std::endl;
//      }
//      
//  }//end for loop over addback hits
//
//  return 0;
//
//}//end HandleUngated

int HandleCaesar(TRuntimeObjects& obj,GCutG *incoming,GCutG *outgoing) {
   
  TS800    *s800    = obj.GetDetector<TS800>();
  TCaesar  *caesar  = obj.GetDetector<TCaesar>();
  TCaesar  *caesar_ab = new TCaesar();//will be used to do addback correction 
                                      //only for the hits inside the time cut
  TCaesar  *caesar_ab_te = new TCaesar();

  std::string dirname; 
  std::string histname;
  std::string dirname2="";
  int cut=-1;

  //this loads the time-energy cut with the same title as the outgoing PID blob
  for(unsigned int x=0;x<timeenergy_cuts.size();x++) {
    if(strcmp(timeenergy_cuts.at(x)->GetTitle(),outgoing->GetName()) == 0){
      dirname2 = Form("caesar_%s_timeenergygated",outgoing->GetName());
      cut = x;
    }
  }

  //this enforces that the outgoing PID blob (AX_from_BY_incoming) is in the incoming gate (BY_incoming)
  const char *outgoing_name = outgoing->GetName();
  const char *toRemove = outgoing->GetTitle(); //going to remove AX from name
  size_t length = strlen(toRemove); //length of AX (can be 23Al or 26P, for example)
  outgoing_name += length+6; //remove AX_from_  
  if(strcmp(outgoing_name,incoming->GetName()) != 0)
    return false;

  dirname = Form("caesar_%s",outgoing->GetName());
  if(s800){
    obj.FillHistogram(dirname,"DTA_for_all_S800_hits",512,-0.2,0.2,s800->GetDta());
    obj.FillHistogram(dirname,"S800_YTA_for_all_S800_hits",1000,-50,50,s800->GetYta());
    obj.FillHistogram(dirname,"TrigBit",30,0,30,s800->GetTrigger().GetRegistr());
    histname = "IC_vs_OBJTOF_PID";
    obj.FillHistogram(dirname,histname,
        	      1000,-3000,-1000,s800->GetMTofObjE1(),
        	      1000,-100,4000,s800->GetIonChamber().GetAve());
  }

  if(!s800 || !caesar)
    return false;

  double beta    = GValue::Value(Form("BETA_%s",outgoing->GetTitle()));
  //double gamma   = 1/(sqrt(1-pow(beta,2)));
  double z_shift = GValue::Value("TARGET_SHIFT_Z");
  TVector3 track = s800->Track(); 

  const int SINGLES_ENERGY_THRESHOLD = 150;
  const int AB_ENERGY_THRESHOLD = 0;

  std::vector<double> energies_singles;
  std::vector<double> energies_singles_te;
  std::vector<double> energies_singles_nondop;
  std::vector<double> energies_singles_nondop_te;
  std::vector<double> energies_addback;
  std::vector<double> energies_addback_te;
  std::vector<double> energies_addback_n0;
  std::vector<double> energies_addback_n1;
  std::vector<double> energies_addback_n2;
  std::vector<double> energies_addback_ng;
  std::vector<double> energies_addback_n0_te;
  std::vector<double> energies_addback_n1_te;
  std::vector<double> energies_addback_n2_te;
  std::vector<double> energies_addback_ng_te;
  std::vector<double> time_singles;
  std::vector<double> time_singles_te;
  std::vector<TVector3> pos_singles;
  std::vector<TVector3> pos_singles_te;
  std::vector<double> detnum_singles;
  std::vector<double> vsnnum_singles;
  energies_singles.clear();
  energies_singles_te.clear();
  energies_addback.clear();
  energies_addback_te.clear();
  double calorimeter = 0;
  int counter = 0;
  int counter2 = 0;

    
  for(unsigned int y=0;y<caesar->Size();y++) {
      TCaesarHit &hit = caesar->GetCaesarHit(y);

      if(hit.IsOverflow())
        continue;

      if (hit.IsValid()){//only accept hits with both times and energies
     
          double energy_dc = hit.GetDoppler(beta,z_shift,&track);        
          double corr_time = caesar->GetCorrTime(hit,s800);

          caesar_ab->InsertHit(hit);
          if (energy_dc > SINGLES_ENERGY_THRESHOLD){              
                 energies_singles.push_back(energy_dc);
                 energies_singles_nondop.push_back(hit.GetEnergy());
                 time_singles.push_back(hit.Time());
                 pos_singles.push_back(hit.GetPosition(z_shift));
                 detnum_singles.push_back(hit.GetAbsoluteDetectorNumber());
                 vsnnum_singles.push_back(hit.GetVSN());
          }//No time-energy gate

          dirname = Form("caesar_%s",outgoing->GetName());
          histname = "energy_dc_pid_in_tcut";
          obj.FillHistogram(dirname,histname,
                            1024,0,8192,energy_dc);
          histname = "energy_lab_pid_in_tcut";
          obj.FillHistogram(dirname,histname,
                            1024,0,8192,hit.GetEnergy());
          if(s800->Track().Theta()<0.04){
            obj.FillHistogram(dirname,"energy_dc_scatterangle<40mrad",
                              1024,0,8192,energy_dc);
          }
          obj.FillHistogram(dirname,"DTA_energy_dc_pid_in_tcut",
			     512,-0.2,0.2,s800->GetDta(),
                             1024,0,8192,energy_dc);
          
          obj.FillHistogram(dirname,"corrtime_vs_doppler",2000,-2000,2000,corr_time,
                                                         1024,0,8192,energy_dc);

          obj.FillHistogram(dirname,"scatterangle_vs_doppler",1000,0,1,s800->Track().Theta(),
                                                         1024,0,8192,energy_dc);

          histname = "detnum_vs_doppler";
          obj.FillHistogram(dirname,histname,200,0,200,hit.GetAbsoluteDetectorNumber(),
			                     1024,0,8192,energy_dc);
          histname = "detnum_vs_energy";
          obj.FillHistogram(dirname,histname,200,0,200,hit.GetAbsoluteDetectorNumber(),
			                     1024,0,8192,hit.GetEnergy());

          if(y == 0){
            calorimeter = 0;
          }
          calorimeter += energy_dc;
          counter++;
          if(y == (caesar->Size() - 1)){
            obj.FillHistogram(dirname,"calorimeter_multiplicity",192,0,192,counter, 
                                                                 1024,0,8192,calorimeter);  
          }
  
          for(int beta_i=0;beta_i<500;beta_i++){
            double beta_use = 0.001*double(beta_i);
            histname = "FindBeta_Doppler";
            obj.FillHistogram(dirname,histname,
		        500,0,0.499,beta_use,
			1024,0,8192,hit.GetDoppler(beta_use,z_shift,&track));  
          }

          for(int z_i=0;z_i<601;z_i++){
            double z_use = -3 + z_i*0.01;
            TVector3 vec(0,0,z_use);
            obj.FillHistogram(dirname,"FindZShift_Doppler",
                      601,-3,3,z_use,
                      1024,0,8192,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vec).Angle(track)))
);  
          }

          obj.FillHistogram(dirname,"CRDC1Y",500,-500,500,s800->GetCrdc(0).GetNonDispersiveY());
          obj.FillHistogram(dirname,"CRDC2Y",500,-500,500,s800->GetCrdc(1).GetNonDispersiveY());
          obj.FillHistogram(dirname,"CRDC1_MaxPad",300,0,300,s800->GetCrdc(0).GetMaxPad(),
                                           2000,0,8000, s800->GetCrdc(0).GetMaxPadSum());
          obj.FillHistogram(dirname,"CRDC2_MaxPad",300,0,300,s800->GetCrdc(1).GetMaxPad(),
                                           2000,0,8000, s800->GetCrdc(1).GetMaxPadSum());
          obj.FillHistogram(dirname,"ATA_vs_BTA",1000,-0.2,0.2,s800->GetAta(),
			                 1000,-0.2,0.2,s800->GetBta());

          obj.FillHistogram(dirname,"S800_YTA",1000,-50,50,s800->GetYta());

          for(unsigned int i=0; i<s800->GetHodoscope().Size(); i++){
            histname = "Hodoscope_CAESAR";
            obj.FillHistogram(dirname,histname,1024,0,8192,energy_dc,
                                               1024,0,8192,s800->GetHodoscope().GetHodoHit(i).GetEnergy());
          }

          if (dirname2.length() && timeenergy_cuts.at(cut)->IsInside(corr_time,hit.GetDoppler(beta,z_shift,&track))){
           
              caesar_ab_te->InsertHit(hit);
              if (energy_dc > SINGLES_ENERGY_THRESHOLD){
           
                 energies_singles_te.push_back(energy_dc);
                 energies_singles_nondop_te.push_back(hit.GetEnergy());
                 time_singles_te.push_back(hit.Time());
                 pos_singles_te.push_back(hit.GetPosition(z_shift));
                
              }//For multiplicity purposes

              dirname = Form("caesar_te_%s",outgoing->GetName());
              histname = "energy_dc_pid_in_tcut";
              obj.FillHistogram(dirname,histname,
                                1024,0,8192,energy_dc);
              histname = "energy_lab_pid_in_tcut";
              obj.FillHistogram(dirname,histname,
                                1024,0,8192,hit.GetEnergy());
              if(s800->Track().Theta()<0.04){
                obj.FillHistogram(dirname,"energy_dc_scatterangle<40mrad",
                                  1024,0,8192,energy_dc);
              }
              obj.FillHistogram(dirname,"DTA_energy_dc_pid_in_tcut",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energy_dc);

              obj.FillHistogram(dirname,"corrtime_vs_doppler",2000,-2000,2000,corr_time,
                                                              1024,0,8192,energy_dc);

              obj.FillHistogram(dirname,"scatterangle_vs_doppler",1000,0,1,s800->Track().Theta(),
                                                         1024,0,8192,energy_dc);
              
              histname = "detnum_vs_doppler";
              obj.FillHistogram(dirname,histname,200,0,200,hit.GetAbsoluteDetectorNumber(),
			                         1024,0,8192,energy_dc);
              histname = "detnum_vs_energy";
              obj.FillHistogram(dirname,histname,200,0,200,hit.GetAbsoluteDetectorNumber(),
			                         1024,0,8192,hit.GetEnergy());

              if(y == 0){
                calorimeter = 0;
              }
              calorimeter += energy_dc;
              counter2++;
              if(y == (caesar->Size() - 1)){
                obj.FillHistogram(dirname,"calorimeter_multiplicity",192,0,192,counter2,   1024,0,8192,calorimeter);           
              }

              for(int beta_i=0;beta_i<500;beta_i++){
                double beta_use = 0.001*double(beta_i);
                histname = "FindBeta_Doppler";
                obj.FillHistogram(dirname,histname,
		            500,0,0.499,beta_use,
			    1024,0,8192,hit.GetDoppler(beta_use,z_shift,&track));  
              }

              for(int z_i=0;z_i<601;z_i++){
                double z_use = -3 + z_i*0.01;
                TVector3 vec(0,0,z_use);
                obj.FillHistogram(dirname,"FindZShift_Doppler",
                          601,-3,3,z_use,
                          1024,0,8192,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vec).Angle(track)))
);  
              }

              obj.FillHistogram(dirname,"CRDC1Y",500,-500,500,s800->GetCrdc(0).GetNonDispersiveY());
              obj.FillHistogram(dirname,"CRDC2Y",500,-500,500,s800->GetCrdc(1).GetNonDispersiveY());
              obj.FillHistogram(dirname,"CRDC1_MaxPad",300,0,300,s800->GetCrdc(0).GetMaxPad(),
                                           2000,0,8000, s800->GetCrdc(0).GetMaxPadSum());
              obj.FillHistogram(dirname,"CRDC2_MaxPad",300,0,300,s800->GetCrdc(1).GetMaxPad(),
                                           2000,0,8000, s800->GetCrdc(1).GetMaxPadSum());
              obj.FillHistogram(dirname,"ATA_vs_BTA",1000,-0.2,0.2,s800->GetAta(),
			                 1000,-0.2,0.2,s800->GetBta());

              obj.FillHistogram(dirname,"S800_YTA",1000,-50,50,s800->GetYta());

              for(unsigned int i=0; i<s800->GetHodoscope().Size(); i++){
                histname = "Hodoscope_CAESAR";
                obj.FillHistogram(dirname,histname,1024,0,8192,energy_dc,
                                                   1024,0,8192,s800->GetHodoscope().GetHodoHit(i).GetEnergy());
          }

          }//is inside time-energy cut     
       
      }//hit has both energy and time
  }//for loop over singles hits


  //Now loop over addback hits
  int num_addback_hits = caesar_ab->AddbackSize();
  for (int y = 0; y < num_addback_hits; y++){
      TCaesarHit &hit = caesar_ab->GetAddbackHit(y);

      //if(hit.IsOverflow())
        //continue; //already took overflow out

      if (hit.IsValid()){//only accept hits with both times and energies

              double energy_dc = hit.GetDoppler(beta,z_shift,&track);
              //double corr_time = caesar_ab->GetCorrTime(hit,s800);
            
              if (energy_dc > AB_ENERGY_THRESHOLD){    
        
                 energies_addback.push_back(energy_dc);   
            
              }//For multiplicity purposes

              dirname = Form("caesar_addback_%s",outgoing->GetName());
              histname = "ab_energy_dc_pid_in_tcut";
              obj.FillHistogram(dirname,histname,
                                1024,0,8192,energy_dc);

              obj.FillHistogram(dirname,"ab_energy_dc_pid_in_tcut_vs_scatterangle",1000,0,1,s800->Track().Theta(),
                                                                                   1024,0,8192,energy_dc);
              
              obj.FillHistogram(dirname,"DTA_ab_energy_dc_pid_in_tcut",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energy_dc);

              if(y == 0){
                calorimeter = 0;
              }
              calorimeter += energy_dc;
              if(y == (num_addback_hits - 1)){
                obj.FillHistogram(dirname,"calorimeter_multiplicity",192,0,192,num_addback_hits,1024,0,8192,calorimeter);           
              }

              if (hit.GetNumHitsContained() == 1 && !hit.is_garbage_addback){
                dirname = Form("caesar_addback_%s",outgoing->GetName());
                histname = "ab_energy_dc_pid_in_tcut_n0";
                obj.FillHistogram(dirname,histname,
                                  1024,0,8192,energy_dc);
                energies_addback_n0.push_back(energy_dc);

              }
              else if (hit.GetNumHitsContained() == 2 && !hit.is_garbage_addback){
                dirname = Form("caesar_addback_%s",outgoing->GetName());
                histname = "ab_energy_dc_pid_in_tcut_n1";
                obj.FillHistogram(dirname,histname,
                                  1024,0,8192,energy_dc);
                energies_addback_n1.push_back(energy_dc);

              }
              else if (hit.GetNumHitsContained() == 3 && !hit.is_garbage_addback){
                dirname = Form("caesar_addback_%s",outgoing->GetName());
                histname = "ab_energy_dc_pid_in_tcut_n2";
                obj.FillHistogram(dirname,histname,
                                  1024,0,8192,energy_dc);
                energies_addback_n2.push_back(energy_dc);
              }
              else if(hit.is_garbage_addback){
                dirname = Form("caesar_addback_%s",outgoing->GetName());
                histname = "ab_energy_dc_pid_in_tcut_ng";
                obj.FillHistogram(dirname,histname,
                                  1024,0,8192,energy_dc);
                energies_addback_ng.push_back(energy_dc);

              }
              else {
                std::cout << "Weird event not meeting any criteria for addback" << std::endl;
                std::cout << "hit.is_garbage_addback    = " << hit.is_garbage_addback << std::endl;
                std::cout << "hit.GetNumHitsContained() = " << hit.GetNumHitsContained() << std::endl;
              }
         
      }//hit has both energy and time
  }//end for loop over addback hits   

  //Fill multiplicity and coincidence histograms for singles
  unsigned int num_hits_singles = energies_singles.size();
  if (num_hits_singles == 1){
      dirname = Form("caesar_%s",outgoing->GetName());
      //histname = "energy_dc_mult_one";
      //obj.FillHistogram(dirname,histname,
      //                  1024,0,8192,energies_singles.at(0));
      //obj.FillHistogram(dirname,"angle_vs_energy_mult_one",180,0,180,TMath::RadToDeg()*pos_singles.at(0).Angle(track),
                                                         //1024,0,8192,energies_singles_nondop.at(0));
      obj.FillHistogram(dirname,"DTA_energy_dc_mult_one",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energies_singles.at(0));
  }
  //else if (num_hits_singles == 2){
      //dirname = Form("caesar_%s",outgoing->GetName());
      //histname = "energy_dc_mult_two";
      //obj.FillHistogram(dirname,histname,
      //                  1024,0,8192,energies_singles.at(0));
      //obj.FillHistogram(dirname,histname,
      //                  1024,0,8192,energies_singles.at(1));
      //obj.FillHistogram(dirname,"DTA_energy_dc_mult_two",
      //			        512,-0.2,0.2,s800->GetDta(),
      //                          1024,0,8192,energies_singles.at(0));
      //obj.FillHistogram(dirname,"DTA_energy_dc_mult_two",
      //		        512,-0.2,0.2,s800->GetDta(),
      //                          1024,0,8192,energies_singles.at(1));
  //}
  for (unsigned int i = 0; i < num_hits_singles; i++){
      for (unsigned int j = i+1; j < num_hits_singles; j++){
        dirname = Form("caesar_%s",outgoing->GetName());
        histname = "energy_dc_coincidence_matrix";
        obj.FillHistogram(dirname, histname,
                          1024,0,8192, energies_singles.at(i),
                          1024,0,8192, energies_singles.at(j));
        obj.FillHistogram(dirname, histname,
                          1024,0,8192, energies_singles.at(j),
                          1024,0,8192, energies_singles.at(i));
        histname = "poss_diff";
        obj.FillHistogram(dirname, histname,
                          500,0,1000,  (pos_singles.at(i)-pos_singles.at(j)).Mag());

        histname = "time_diff";
        obj.FillHistogram(dirname, histname,
                          1000,-1000,1000, time_singles.at(i)-time_singles.at(j));
        if(num_hits_singles==2){
          histname = "energy_dc_coincidence_matrix_multtwo";
          obj.FillHistogram(dirname, histname,
                            1024,0,8192, energies_singles.at(i),
                            1024,0,8192, energies_singles.at(j));
          obj.FillHistogram(dirname, histname,
                            1024,0,8192, energies_singles.at(j),
                            1024,0,8192, energies_singles.at(i));
        }
      }

    dirname = Form("caesar_%s",outgoing->GetName());
    histname = "multiplicity_energydc";
    obj.FillHistogram(dirname, histname,
                      192, 0, 192, num_hits_singles,
                      1024,0,8192, energies_singles.at(i));
  }

  if(num_hits_singles>0){  
    obj.FillHistogram(dirname,"DTA_energy_dc_first_hit_only",
		      512,-0.2,0.2,s800->GetDta(),
                      1024,0,8192,energies_singles.at(0));
  }

  if(detnum_singles.size()==2){
    dirname = Form("caesar_%s",outgoing->GetName());
    histname = "mult2_detnum_detnum";
    obj.FillHistogram(dirname, histname,
                      192, 0, 192, detnum_singles.at(0),
                      192, 0, 192, detnum_singles.at(1));
    obj.FillHistogram(dirname, histname,
                      192, 0, 192, detnum_singles.at(1),
                      192, 0, 192, detnum_singles.at(0));
    histname = "mult2_vsnnum_vsnnum";
    obj.FillHistogram(dirname, histname,
                      20, 0, 20, vsnnum_singles.at(0),
                      20, 0, 20, vsnnum_singles.at(1));
    obj.FillHistogram(dirname, histname,
                      20, 0, 20, vsnnum_singles.at(1),
                      20, 0, 20, vsnnum_singles.at(0));
  }

  //addback_mult == multiplicity in gates
  int addback_mult = energies_addback.size();
  int n0_mult = energies_addback_n0.size();
  int n1_mult = energies_addback_n1.size();
  int n2_mult = energies_addback_n2.size();
  int ng_mult = energies_addback_ng.size();
  if (addback_mult == 1){
      dirname = Form("caesar_addback_%s",outgoing->GetName());
      histname = "ab_energy_dc_mult_one";
      obj.FillHistogram(dirname,histname,
                        1024,0,8192,energies_addback.at(0));
      obj.FillHistogram(dirname,"DTA_ab_energy_dc_mult_one",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energies_addback.at(0));
      if (n0_mult ==1){
        histname = "ab_energy_dc_mult_one_n0";
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback_n0.at(0));
      }
      else if (n1_mult ==1){
        histname = "ab_energy_dc_mult_one_n1";
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback_n1.at(0));
      }
      else if (n2_mult ==1){
        histname = "ab_energy_dc_mult_one_n2";
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback_n2.at(0));
      }
  }//end addback_mult == 1
  else if (addback_mult == 2){
      dirname = Form("caesar_addback_%s",outgoing->GetName());
      histname = "ab_energy_dc_mult_two";
      obj.FillHistogram(dirname,histname,
                        1024,0,8192,energies_addback.at(0));
      obj.FillHistogram(dirname,histname,
                        1024,0,8192,energies_addback.at(1));
      obj.FillHistogram(dirname,"DTA_ab_energy_dc_mult_two",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energies_addback.at(0));
      obj.FillHistogram(dirname,"DTA_ab_energy_dc_mult_two",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energies_addback.at(1));
      for (int i = 0; i < addback_mult; i++){
        for (int j = i+1; j < addback_mult; j++){
          dirname = Form("caesar_addback_%s",outgoing->GetName());
          histname = "ab_energy_dc_coincidence_matrix_multtwo";
          obj.FillHistogram(dirname, histname,
              1024,0,8192, energies_addback.at(i),
              1024,0,8192, energies_addback.at(j));
          obj.FillHistogram(dirname, histname,
              1024,0,8192, energies_addback.at(j),
              1024,0,8192, energies_addback.at(i));

          if (n0_mult + n1_mult + n2_mult == 2){
            histname = "ab_energy_dc_coincidence_matrix_multtwo_nogarbage";
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback.at(i),
                1024,0,8192, energies_addback.at(j));
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback.at(j),
                1024,0,8192, energies_addback.at(i));
          }

          if (n0_mult == 2){
            histname = "ab_energy_dc_coincidence_matrix_multtwo_n0";
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback.at(i),
                1024,0,8192, energies_addback.at(j));
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback.at(j),
                1024,0,8192, energies_addback.at(i));
          }
          if (n0_mult + n1_mult == 2){
            histname = "ab_energy_dc_coincidence_matrix_multtwo_n0n1";
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback.at(i),
                1024,0,8192, energies_addback.at(j));
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback.at(j),
                1024,0,8192, energies_addback.at(i));
            //if(s800->Track().Theta()<0.04){
              //histname = "ab_energy_dc_coincidence_matrix_multtwo_n0n1<40mrad";  
              //obj.FillHistogram(dirname,histname,
                //1024,0,8192, energies_addback_te.at(i),
                //1024,0,8192, energies_addback_te.at(j));
              //obj.FillHistogram(dirname,histname,
                //1024,0,8192, energies_addback_te.at(j),
                //1024,0,8192, energies_addback_te.at(i));
            //}
          }
        }//end for
      }//end for
  
      if (n0_mult + n1_mult + n2_mult == 2){
        histname = "ab_energy_dc_multtwo_nogarbage";
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback.at(0));
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback.at(1));
        if (n2_mult ==0){
          if (n1_mult == 0){
            histname = "ab_energy_dc_multtwo_n0";
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback.at(0));
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback.at(1));
          }//pure n0
          else{
            histname = "ab_energy_dc_multtwo_n1";
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback.at(0));
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback.at(1));
          }
        }//no n2 events
        else{
            histname = "ab_energy_dc_multtwo_n2";
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback.at(0));
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback.at(1));
        }
      } 
  }//addback_mult == 2

  for (int i = 0; i < addback_mult; i++){
      for (int j = i+1; j < addback_mult; j++){
        dirname = Form("caesar_addback_%s",outgoing->GetName());
        histname = "ab_energy_dc_coincidence_matrix";
        obj.FillHistogram(dirname, histname,
                          1024,0,8192, energies_addback.at(i),
                          1024,0,8192, energies_addback.at(j));
        obj.FillHistogram(dirname, histname,
                          1024,0,8192, energies_addback.at(j),
                          1024,0,8192, energies_addback.at(i));
        if (ng_mult == 0){
          histname = "ab_energy_dc_coincidence_matrix_nogarbage";
          obj.FillHistogram(dirname, histname,
              1024,0,8192, energies_addback.at(i),
              1024,0,8192, energies_addback.at(j));
          obj.FillHistogram(dirname, histname,
              1024,0,8192, energies_addback.at(j),
              1024,0,8192, energies_addback.at(i));
        }
      }

    dirname = Form("caesar_addback_%s",outgoing->GetName());
    histname = "ab_multiplicity_energydc";
    obj.FillHistogram(dirname, histname,
                      192, 0, 192, addback_mult,
                      1024,0,8192, energies_addback.at(i));
  }

  //dirname = Form("caesar_addback_%s",outgoing->GetName());
  //histname = "ab_multiplicity";
  //obj.FillHistogram(dirname, histname,
  //                  192, 0, 192, addback_mult);
  if(addback_mult>0){
    obj.FillHistogram(dirname,"DTA_energy_dc_first_hit_only",
                      512,-0.2,0.2,s800->GetDta(),
                      1024,0,8192,energies_addback.at(0));
  }

  //Now loop over addback_te hits
  int num_addback_hits_te = caesar_ab_te->AddbackSize();
  for (int y = 0; y < num_addback_hits_te; y++){
      TCaesarHit &hit = caesar_ab_te->GetAddbackHit(y);

      //if(hit.IsOverflow())
        //continue; //already took overflow out

      if (hit.IsValid()){//only accept hits with both times and energies

          double energy_dc = hit.GetDoppler(beta,z_shift,&track);
          double corr_time = caesar_ab_te->GetCorrTime(hit,s800);

   
          if (dirname2.length() && timeenergy_cuts.at(cut)->IsInside(corr_time,hit.GetDoppler(beta,z_shift,&track))){
            
              if (energy_dc > AB_ENERGY_THRESHOLD){    
        
                 energies_addback_te.push_back(energy_dc);   
            
              }//For multiplicity purposes

              dirname = Form("caesar_addback_te_%s",outgoing->GetName());
              histname = "ab_energy_dc_pid_in_tcut";
              obj.FillHistogram(dirname,histname,
                                1024,0,8192,energy_dc);
              
              obj.FillHistogram(dirname,"ab_energy_dc_pid_in_tcut_vs_scatterangle",1000,0,1,s800->Track().Theta(),
                                                                                   1024,0,8192,energy_dc);
              
              obj.FillHistogram(dirname,"DTA_ab_energy_dc_pid_in_tcut",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energy_dc);

              if(y == 0){
                calorimeter = 0;
              }
              calorimeter += energy_dc;
              if(y == (num_addback_hits - 1)){
                obj.FillHistogram(dirname,"ab_calorimeter_multiplicity",192,0,192,num_addback_hits_te,
                                       1024,0,8192,calorimeter);          
              }

              if (hit.GetNumHitsContained() == 1 && !hit.is_garbage_addback){
                dirname = Form("caesar_addback_te_%s",outgoing->GetName());
                histname = "ab_energy_dc_pid_in_tcut_n0";
                obj.FillHistogram(dirname,histname,
                                  1024,0,8192,energy_dc);
                energies_addback_n0_te.push_back(energy_dc);

              }
              else if (hit.GetNumHitsContained() == 2 && !hit.is_garbage_addback){
                dirname = Form("caesar_addback_te_%s",outgoing->GetName());
                histname = "ab_energy_dc_pid_in_tcut_n1";
                obj.FillHistogram(dirname,histname,
                                  1024,0,8192,energy_dc);
                energies_addback_n1_te.push_back(energy_dc);

              }
              else if (hit.GetNumHitsContained() == 3 && !hit.is_garbage_addback){
                dirname = Form("caesar_addback_te_%s",outgoing->GetName());
                histname = "ab_energy_dc_pid_in_tcut_n2";
                obj.FillHistogram(dirname,histname,
                                  1024,0,8192,energy_dc);
                energies_addback_n2_te.push_back(energy_dc);
              }
              else if(hit.is_garbage_addback){
                dirname = Form("caesar_addback_te_%s",outgoing->GetName());
                histname = "ab_energy_dc_pid_in_tcut_ng";
                obj.FillHistogram(dirname,histname,
                                  1024,0,8192,energy_dc);
                energies_addback_ng_te.push_back(energy_dc);

              }
              else {
                std::cout << "Weird event not meeting any criteria for addback" << std::endl;
                std::cout << "hit.is_garbage_addback    = " << hit.is_garbage_addback << std::endl;
                std::cout << "hit.GetNumHitsContained() = " << hit.GetNumHitsContained() << std::endl;
              }
          }//is inside time-energy cut
      }//hit has both energy and time
  }//end for loop over addback_te hits

  //Fill multiplicity and coincidence histograms for singles_te
  unsigned int num_hits_singles_te = energies_singles_te.size();
  if (num_hits_singles_te == 1){
      dirname = Form("caesar_te_%s",outgoing->GetName());
      //histname = "energy_dc_mult_one";
      //obj.FillHistogram(dirname,histname,
      //                  1024,0,8192,energies_singles_te.at(0));
      //obj.FillHistogram(dirname,"angle_vs_energy_mult_one",180,0,180,TMath::RadToDeg()*pos_singles_te.at(0).Angle(track),
                                                         //1024,0,8192,energies_singles_nondop_te.at(0));
      obj.FillHistogram(dirname,"DTA_energy_dc_mult_one",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energies_singles.at(0));
  }
  //else if (num_hits_singles_te == 2){
      //dirname = Form("caesar_te_%s",outgoing->GetName());
      //histname = "energy_dc_mult_two";
      //obj.FillHistogram(dirname,histname,
      //                  1024,0,8192,energies_singles_te.at(0));
      //obj.FillHistogram(dirname,histname,
      //                  1024,0,8192,energies_singles_te.at(1));
  //}
  for (unsigned int i = 0; i < num_hits_singles_te; i++){
      for (unsigned int j = i+1; j < num_hits_singles_te; j++){
        dirname = Form("caesar_te_%s",outgoing->GetName());
        histname = "energy_dc_coincidence_matrix";
        obj.FillHistogram(dirname, histname,
                          1024,0,8192, energies_singles_te.at(i),
                          1024,0,8192, energies_singles_te.at(j));
        obj.FillHistogram(dirname, histname,
                          1024,0,8192, energies_singles_te.at(j),
                          1024,0,8192, energies_singles_te.at(i));
        histname = "poss_diff";
        obj.FillHistogram(dirname, histname,
                          500,0,1000,  (pos_singles_te.at(i)-pos_singles_te.at(j)).Mag());

        histname = "time_diff";
        obj.FillHistogram(dirname, histname,
                          1000,-1000,1000, time_singles_te.at(i)-time_singles_te.at(j));
        if(num_hits_singles_te==2){
          histname = "energy_dc_coincidence_matrix_multtwo";
          obj.FillHistogram(dirname, histname,
                            1024,0,8192, energies_singles_te.at(i),
                            1024,0,8192, energies_singles_te.at(j));
          obj.FillHistogram(dirname, histname,
                            1024,0,8192, energies_singles_te.at(j),
                            1024,0,8192, energies_singles_te.at(i));
        }
      }

    dirname = Form("caesar_te_%s",outgoing->GetName());
    histname = "multiplicity_energydc";
    obj.FillHistogram(dirname, histname,
                      192, 0, 192, num_hits_singles_te,
                      1024,0,8192, energies_singles_te.at(i));
  }

  //dirname = Form("caesar_te_%s",outgoing->GetName());
  //histname = "multiplicity";
  //obj.FillHistogram(dirname, histname,
  //                  192, 0, 192, num_hits_singles_te);

  if(num_hits_singles_te>0){
    obj.FillHistogram(dirname,"DTA_energy_dc_first_hit_only",
		      512,-0.2,0.2,s800->GetDta(),
                      1024,0,8192,energies_singles_te.at(0));
  }

  //addback_mult_te == multiplicity in gates with time-energy cut
  int addback_mult_te = energies_addback_te.size();
  int n0_mult_te = energies_addback_n0_te.size();
  int n1_mult_te = energies_addback_n1_te.size();
  int n2_mult_te = energies_addback_n2_te.size();
  int ng_mult_te = energies_addback_ng_te.size();
  if (addback_mult_te == 1){
      dirname = Form("caesar_addback_te_%s",outgoing->GetName());
      histname = "ab_energy_dc_mult_one";
      obj.FillHistogram(dirname,histname,
                        1024,0,8192,energies_addback_te.at(0));
      obj.FillHistogram(dirname,"DTA_ab_energy_dc_mult_one",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energies_addback_te.at(0));
      if (n0_mult_te ==1){
        histname = "ab_energy_dc_mult_one_n0";
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback_n0_te.at(0));
      }
      else if (n1_mult_te ==1){
        histname = "ab_energy_dc_mult_one_n1";
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback_n1_te.at(0));
      }
      else if (n2_mult_te ==1){
        histname = "ab_energy_dc_mult_one_n2";
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback_n2_te.at(0));
      }
  }
  else if (addback_mult_te == 2){
      dirname = Form("caesar_addback_te_%s",outgoing->GetName());
      histname = "ab_energy_dc_mult_two";
      obj.FillHistogram(dirname,histname,
                        1024,0,8192,energies_addback_te.at(0));
      obj.FillHistogram(dirname,histname,
                        1024,0,8192,energies_addback_te.at(1));
      obj.FillHistogram(dirname,"DTA_ab_energy_dc_mult_two",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energies_addback_te.at(0));
      obj.FillHistogram(dirname,"DTA_ab_energy_dc_mult_two",
			        512,-0.2,0.2,s800->GetDta(),
                                1024,0,8192,energies_addback_te.at(1));
      for (int i = 0; i < addback_mult_te; i++){
        for (int j = i+1; j < addback_mult_te; j++){
          dirname = Form("caesar_addback_te_%s",outgoing->GetName());
          histname = "ab_energy_dc_coincidence_matrix_multtwo";
          obj.FillHistogram(dirname, histname,
              1024,0,8192, energies_addback_te.at(i),
              1024,0,8192, energies_addback_te.at(j));
          obj.FillHistogram(dirname, histname,
              1024,0,8192, energies_addback_te.at(j),
              1024,0,8192, energies_addback_te.at(i));

          if (n0_mult_te + n1_mult_te + n2_mult_te == 2){
            histname = "ab_energy_dc_coincidence_matrix_multtwo_nogarbage";
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback_te.at(i),
                1024,0,8192, energies_addback_te.at(j));
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback_te.at(j),
                1024,0,8192, energies_addback_te.at(i));
          }

          if (n0_mult_te == 2){
            histname = "ab_energy_dc_coincidence_matrix_multtwo_n0";
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback_te.at(i),
                1024,0,8192, energies_addback_te.at(j));
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback_te.at(j),
                1024,0,8192, energies_addback_te.at(i));
          }
          if (n0_mult_te + n1_mult_te == 2){
            histname = "ab_energy_dc_coincidence_matrix_multtwo_n0n1";
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback_te.at(i),
                1024,0,8192, energies_addback_te.at(j));
            obj.FillHistogram(dirname,histname,
                1024,0,8192, energies_addback_te.at(j),
                1024,0,8192, energies_addback_te.at(i));
            //if(s800->Track().Theta()<0.04){
              //histname = "ab_energy_dc_coincidence_matrix_multtwo_n0n1<40mrad";  
              //obj.FillHistogram(dirname,histname,
                //1024,0,8192, energies_addback_te.at(i),
                //1024,0,8192, energies_addback_te.at(j));
              //obj.FillHistogram(dirname,histname,
                //1024,0,8192, energies_addback_te.at(j),
                //1024,0,8192, energies_addback_te.at(i));
            //}
          }
        }
      }
      if (n0_mult_te + n1_mult_te + n2_mult_te == 2){
        histname = "ab_energy_dc_multtwo_nogarbage";
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback_te.at(0));
        obj.FillHistogram(dirname,histname,
                          1024,0,8192,energies_addback_te.at(1));
        if (n2_mult_te ==0){
          if (n1_mult_te == 0){
            histname = "ab_energy_dc_multtwo_n0";
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback_te.at(0));
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback_te.at(1));
          }//pure n0
          else{
            histname = "ab_energy_dc_multtwo_n1";
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback_te.at(0));
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback_te.at(1));
          }
        }//no n2 events
        else{
            histname = "ab_energy_dc_multtwo_n2";
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback_te.at(0));
            obj.FillHistogram(dirname,histname,
                              1024,0,8192,energies_addback_te.at(1));
        }
      }
  }//addback_mult == 2
  for (int i = 0; i < addback_mult_te; i++){
      for (int j = i+1; j < addback_mult_te; j++){
        dirname = Form("caesar_addback_te_%s",outgoing->GetName());
        histname = "ab_energy_dc_coincidence_matrix";
        obj.FillHistogram(dirname, histname,
                          1024,0,8192, energies_addback_te.at(i),
                          1024,0,8192, energies_addback_te.at(j));
        obj.FillHistogram(dirname, histname,
                          1024,0,8192, energies_addback_te.at(j),
                          1024,0,8192, energies_addback_te.at(i));
        if (ng_mult_te == 0){
          histname = "ab_energy_dc_coincidence_matrix_nogarbage";
          obj.FillHistogram(dirname, histname,
              1024,0,8192, energies_addback_te.at(i),
              1024,0,8192, energies_addback_te.at(j));
          obj.FillHistogram(dirname, histname,
              1024,0,8192, energies_addback_te.at(j),
              1024,0,8192, energies_addback_te.at(i));
        }
      }

    dirname = Form("caesar_addback_te_%s",outgoing->GetName());
    histname = "ab_multiplicity_energydc";
    obj.FillHistogram(dirname, histname,
                      192, 0, 192, addback_mult_te,
                      1024,0,8192, energies_addback_te.at(i));
  }

  //dirname = Form("caesar_addback_te_%s",outgoing->GetName());
  //histname = "ab_multiplicity";
  //obj.FillHistogram(dirname, histname,
  //                    192, 0, 192, addback_mult_te);

  if(addback_mult_te>0){
    obj.FillHistogram(dirname,"DTA_energy_dc_first_hit_only",
		      512,-0.2,0.2,s800->GetDta(),
                      1024,0,8192,energies_addback_te.at(0));
  }

  //if(caesar) { caesar->Clear(); }
  //if(caesar_ab) { caesar_ab->Clear(); }
  //if(caesar_ab_te) { caesar_ab_te->Clear(); }
  
  return true;

}//end HandleCaesar


// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {

  TList *list  = &(obj.GetObjects());
  TList *gates = &(obj.GetGates());
  int numobj = list->GetSize();

//  HandleUngated(obj);

  if(gates_loaded!=gates->GetSize()) {
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      if(!tag.compare("incoming")) {
        incoming_cuts.push_back(gate);
      } else if(!tag.compare("outgoing")) {
        outgoing_cuts.push_back(gate);
      } else if(!tag.compare("timeenergy")){
        timeenergy_cuts.push_back(gate);
      } else {
        std::cout << "Unknown gate: " << gate->GetName() << std::endl;
      }
      gates_loaded++;
    }
  }

  //printf("incoming.size() == %i\n",incoming_cuts.size());
  //int incoming_passed=-1;
  //int outgoing_passed=-1;
  std::vector<int> incoming_passed;
  std::vector<int> outgoing_passed;

  for(unsigned int x=0;x<incoming_cuts.size();x++) {
    bool passed = OutgoingBeam(obj,incoming_cuts.at(x)); 
    if(x!=0 && passed) {
      //incoming_passed = x;
      incoming_passed.push_back(x);
      break;
    }
  }
  //for(unsigned int x=0;x<outgoing_cuts.size();x++) {
    //bool passed = IncomingBeam(obj,outgoing_cuts.at(x)); 
    //if(x!=0 && passed) {
      //outgoing_passed = x;
      //break;
    //}
  //} 

  for(unsigned int x=0;x<outgoing_cuts.size();x++) {
    bool passed = true; //false;
    if(x==0 || (x!=0)) { //&& incoming_passed>0)) {
      passed = IncomingBeam(obj,outgoing_cuts.at(x));
    }
    if(x!=0 && passed) {
      outgoing_passed.push_back(x);
      //break;
    }
  }

  //if(incoming_passed>0 && outgoing_passed>0) {
    //HandleCaesar(obj,incoming_cuts.at(incoming_passed),outgoing_cuts.at(outgoing_passed));
  //}

  for(unsigned int i=0; i<incoming_passed.size(); i++){
    for(unsigned int j=0; j<outgoing_passed.size(); j++){
      HandleCaesar(obj,incoming_cuts.at(incoming_passed.at(i)),outgoing_cuts.at(outgoing_passed.at(j)));
      //std::cout << outgoing_cuts.at(outgoing_passed.at(j))->GetName() << std::endl;
    }
  }

  if(numobj!=list->GetSize())
    list->Sort();
}
