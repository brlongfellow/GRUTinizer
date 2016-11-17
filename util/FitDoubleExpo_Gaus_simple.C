#include <iostream>
#include "TF1.h"
#include "TMath.h"
#include "TH1.h"

//Double_t gBgConstant, gBgSlope, gContent, gMean, gContent_1, gMean_1, gContent_2, gMean_2, gSigma, gSigma_1, gSigma_2, gBinW, gChi2pNDF;

//int temp_glob = 0;




class tenHistHolder{
  public:
    tenHistHolder(TH1* hist1, TH1* hist2, TH1* hist3, TH1* hist4, TH1* hist5,
        TH1* hist6, TH1* hist7, TH1* hist8, TH1* hist9, TH1* hist10){
      myupdate=0;
      this->hist1 = hist1;
      this->hist2 = hist2;
      this->hist3 = hist3;
      this->hist4 = hist4;
      this->hist5 = hist5;
      this->hist6 = hist6;
      this->hist7 = hist7;
      this->hist8 = hist8;
      this->hist9 = hist9;
      this->hist10 = hist10;
    }
    
    int myupdate;

    TH1* hist1;
    TH1* hist2;
    TH1* hist3;
    TH1* hist4;
    TH1* hist5;
    TH1* hist6;
    TH1* hist7;
    TH1* hist8;
    TH1* hist9;
    TH1* hist10;
      

    double histValue(TH1* m_hist, double x){
      int binNum = m_hist->GetXaxis()->FindBin(x); //gHist->GetBin() does not respect rebinning.

      int nBins = m_hist->GetNbinsX();
      int kevPerBin = m_hist->GetXaxis()->GetXmax()/nBins;
      int curBinX = m_hist->GetBinCenter(binNum);
      int nextBinX = m_hist->GetBinCenter(binNum+1);
      int prevBinX = m_hist->GetBinCenter(binNum-1);

      if (x > prevBinX && x <= curBinX){
        double leftDiff = x - prevBinX;
        double rightDiff = curBinX - x;

        leftDiff = 1.0 - leftDiff/(double)kevPerBin;   //These numbers are now less than 1
        rightDiff = 1.0 - rightDiff/(double)kevPerBin; //and a measure of how close it is to that bin
        double binContentLeft = m_hist->GetBinContent(binNum-1);
        double binContentRight = m_hist->GetBinContent(binNum);
        return (leftDiff*binContentLeft+rightDiff*binContentRight);
      }

      else if (x > curBinX && x < nextBinX){
        double leftDiff = x - curBinX;
        double rightDiff = nextBinX - x;

        leftDiff = 1.0 - leftDiff/(double)kevPerBin;
        rightDiff = 1.0 - rightDiff/(double)kevPerBin;
        double binContentLeft = m_hist->GetBinContent(binNum);
        double binContentRight = m_hist->GetBinContent(binNum+1);
        return (leftDiff*binContentLeft+rightDiff*binContentRight);
      }
      std::cout << "FAILED IN HISTVALUE!" << std::endl;
      return m_hist->GetBinContent(binNum);
    }

    virtual double operator() (double *x, double *par){
      if((myupdate%500)==0) {
        printf("     myupdate = %i               \r");
        fflush(stdout);
      }
      myupdate++;
      return (par[0]*(histValue(hist1, x[0])+histValue(hist2,x[0])*par[1]) + 
          histValue(hist3,x[0])*par[2] + histValue(hist4,x[0])*par[3] + 
          histValue(hist5,x[0])*par[4] + histValue(hist6,x[0])*par[5] +
          histValue(hist7,x[0])*par[6] + histValue(hist8,x[0])*par[7] +
          histValue(hist9,x[0])*par[8] + histValue(hist10,x[0])*par[9]);
    }
};

class DoubleExpTenHist : tenHistHolder{
  /*
     par[0]  histogram 1 scaling factor
     par[1]  histogram 2 scaling factor
     par[2]  histogram 3 scaling factor
     par[3]  histogram 4 scaling factor
     par[4]  histogram 5 scaling factor
     par[5]  histogram 6 scaling factor
     par[6]  histogram 7 scaling factor
     par[7]  histogram 8 scaling factor
     par[8]  histogram 9 scaling factor
     par[9]  histogram 10 scaling factor
     par[10]  double exponential scaling factor
     par[11]  initial exponential value
     par[12]  exponential decay constant
     par[13]  2nd initial exponential value
     par[14]  2nd exponential decay constant
     */
  public:
    DoubleExpTenHist(TH1 *hist1, TH1 *hist2, TH1 *hist3,TH1 *hist4, TH1 *hist5, TH1 *hist6, TH1* hist7, TH1* hist8, TH1* hist9, TH1* hist10): tenHistHolder(hist1,hist2,hist3,hist4,hist5,hist6,hist7,hist8,hist9,hist10){}
    virtual double operator()(double *x, double *par){
      return (par[0]*histValue(hist1, x[0]) + par[1]*histValue(hist2,x[0]) +
          par[2]*histValue(hist3,x[0]) + par[3]*histValue(hist4,x[0]) + 
          par[4]*histValue(hist5,x[0]) + par[5]*histValue(hist6,x[0]) +
          par[6]*histValue(hist7,x[0]) + par[7]*histValue(hist8,x[0]) +
          par[8]*histValue(hist9,x[0]) + par[9]*histValue(hist10,x[0]) +  
          par[10]*(par[11]*TMath::Exp(par[12]*x[0])+ par[13]*TMath::Exp(par[14]*x[0])));
    }

};

TF1 *FitDoubleExpTenHist(TH1 *hist_to_fit, TH1 *geant_1, TH1 *geant_2,
    TH1 *geant_3, TH1 *geant_4, TH1 *geant_5, TH1* geant_6,
    TH1 *geant_7, TH1 *geant_8, TH1 *geant_9, TH1* geant_10,
    //                         Double_t gLowX, Double_t gUpX, Double_t *init=NULL){
    Double_t gLowX, Double_t gUpX){
      if(gLowX > gUpX){
        std::cout << "Your range is illogical" << std::endl;
        return NULL;
      }

      hist_to_fit->GetXaxis()->SetRangeUser(gLowX,gUpX);
      int gBinW = hist_to_fit->GetBinWidth(1);

      DoubleExpTenHist* deh10 = new DoubleExpTenHist(geant_1, geant_2, geant_3, geant_4, 
                                                     geant_5, geant_6, geant_7, 
                                                     geant_8, geant_9, geant_10);
      TF1* fitfunc = new TF1("double_exp_ten_hist",deh10, 0, 1, 14,"DoubleExpTenHist");
      TF1* fitfunc_to_draw = new TF1("double_exp_ten_hist_check",deh10, 0, 1, 14,"DoubleExpTenHist");

      //if (init==NULL){
      //std::cout << "Need initial parameters. Exiting." << std::endl;
      //return NULL;
      //}
      //if (init!=NULL){
      //fitfunc->SetParameters(init);
      //}

      //fitfunc->FixParameter(1,init[1]);
      //for (int i = 7; i < 11; i++){
      //fitfunc->FixParameter(i, init[i]);
      //}

      ////fitfunc->SetParLimits(0, 9e-04, 1.2e-03);
      //fitfunc->SetParLimits(2, 0, 1);
      //fitfunc->SetParLimits(3, 0, 1);
      ////  fitfunc->SetParLimits(4, 0, 1);
      //fitfunc->SetParLimits(5, 0, 1);
      //fitfunc->FixParameter(4,init[4]);
      ////  fitfunc->FixParameter(2,init[2]);
      ////  fitfunc->FixParameter(5,init[5]);

      ////fitfunc->SetRange(gLowX, gUpX);
      fitfunc->SetLineColor(4);
      fitfunc->SetLineWidth(3);
      hist_to_fit->SetLineColor(kBlue);
      geant_1->SetLineColor(kGreen+4);
      geant_2->SetLineColor(kRed);
      geant_3->SetLineColor(kMagenta);
      geant_4->SetLineColor(kOrange);
      geant_5->SetLineColor(kSpring);
      geant_6->SetLineColor(kViolet);
      geant_7->SetLineColor(kCyan);
      geant_8->SetLineColor(kYellow);
      geant_9->SetLineColor(kPink);
      geant_10->SetLineColor(kTeal);

      fitfunc->SetParameter(0, 0.2);
      fitfunc->SetParameter(1, 0.05);
      fitfunc->SetParameter(2, 0.2);
      fitfunc->SetParameter(3, 0.2);
      fitfunc->SetParameter(4, 1.5);
      fitfunc->SetParameter(5, 1);
      fitfunc->SetParameter(6, 1);
      fitfunc->SetParameter(7, 0.5);
      fitfunc->SetParameter(8, 0.1);
      fitfunc->SetParameter(9, 1);
      fitfunc->SetParameter(10, 1);
      fitfunc->SetParameter(11, 162808);
      fitfunc->SetParameter(12, -0.00210025);
      fitfunc->SetParameter(13, 486.425);
      fitfunc->SetParameter(14, 0.00161493);

      fitfunc->SetParName(0, "1082 Scaling    ");
      fitfunc->SetParName(1, "1320 Scaling    ");
      fitfunc->SetParName(2, "1848 Scaling    ");
      fitfunc->SetParName(3, "1906 Scaling    ");
      fitfunc->SetParName(4, "2360 Scaling    ");
      fitfunc->SetParName(5, "2560 Scaling    ");
      fitfunc->SetParName(6, "2930 Scaling    ");
      fitfunc->SetParName(7, "3754 Scaling    ");
      fitfunc->SetParName(8, "670  Scaling    ");
      fitfunc->SetParName(9, "824  Scaling    ");
      fitfunc->SetParName(10, "Double Exp. Scaling   ");
      fitfunc->SetParName(11, "Exp.Scaling         ");
      fitfunc->SetParName(12, "Exp.Decay.Const       ");
      fitfunc->SetParName(13, "2nd Exp. Scaling    ");
      fitfunc->SetParName(14, "2nd Exp. Decay Const  ");

      //hist_to_fit->Fit(fitfunc, "PRMEQ");
      hist_to_fit->Sumw2();
      hist_to_fit->Fit(fitfunc, "PMEOQ","",gLowX,gUpX);
      fitfunc_to_draw->SetRange(0,4096);
      fitfunc_to_draw->SetParameters(fitfunc->GetParameters());
      fitfunc_to_draw->SetChisquare(fitfunc->GetChisquare());
      fitfunc_to_draw->SetParError(0, fitfunc->GetParError(0));

      return fitfunc_to_draw;
    }
