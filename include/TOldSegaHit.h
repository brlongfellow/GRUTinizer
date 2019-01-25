#ifndef __TOLDSEGAHIT_H__
#define __TOLDSEGAHIT_H__


#include <TDetectorHit.h>

#include <algorithm>
#include <TMath.h>
#include <TOldSega.h>


class TOldSegaHit : public TDetectorHit {

  public: 
    TOldSegaHit() { }
    TOldSegaHit(const TOldSegaHit &rhs) { rhs.Copy(*this); }
  
    void Copy(TObject &rhs) const;      //  { TDetectorHit::Copy(rhs); }
    void Print(Option_t *opt="") const; //  { TDetectorHit::Print(opt); }
    void Clear(Option_t *opt=""); 

    TDetectorHit &GetSegment(int i)  { return segments.at(i); }
    int Size() const                { return segments.size(); }

    void AddSegment(int id, float charge);

    int   GetSegId(int i=0)  const { if(Size()>i) return (segments.at(i).Address())&0x000000ff; return -1; }
    int   GetSegChg(int i)   const { if(Size()>i) return (segments.at(i).Charge()); return -1; }
    float GetSegEng(int i)   const { if(Size()>i) return (segments.at(i).GetEnergy()); return -1.0; }


    int GetDetId()       const { return (Address()&0x0000ff00)>>8; } 
    int GetDetNum()      const { if(TChannel::GetChannel(Address())) return TChannel::GetChannel(Address())->GetNumber(); return -1; }

    void Sort() { std::sort(segments.begin(),segments.end()); }

    TVector3 GetPosition() const;
    //TVector3 GetPosition(double z_shift) const;

    int GetLayer() const;

    double GetDoppler(double beta, const TVector3 *track=0) const {
      if(Size()<1)
        return 0.0;
      if(track==0) {
        track = &BeamUnitVec;
      }
      TChannel *chan = TChannel::GetChannel(Address());
      double theta;
      if(chan && chan->UseThetaOffset()) {
        theta = chan->ThetaOffset();
        //printf("using theta = %f\n ",theta);
      } else {
        theta = 0.0;
      }
      double tmp = 0.0;
      double gamma = 1/(sqrt(1-pow(beta,2)));
      tmp = GetEnergy()*gamma *(1 - beta*TMath::Cos(TOldSegaHit::GetPosition().Angle(*track) + theta));
      return tmp;
    }
    
  private:
    std::vector<TDetectorHit> segments;

  ClassDef(TOldSegaHit,1)
};



#endif

