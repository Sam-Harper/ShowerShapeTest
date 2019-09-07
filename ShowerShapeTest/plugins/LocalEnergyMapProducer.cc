#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TFile.h"
#include "TH2D.h"
#include <string>


class LocalEnergyMapProducer : public edm::EDAnalyzer {

private:
  bool fillFromEles_;
  bool fillFromPhos_;
  bool fillFromSCs_;
  bool genMatch_;

  edm::EDGetTokenT<edm::View<reco::GsfElectron> > elesToken_;
  edm::EDGetTokenT<edm::View<reco::Photon> > phosToken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebSCsToken_;
  std::vector<edm::EDGetTokenT<reco::SuperClusterCollection>> scTokens_;
  edm::EDGetTokenT<EcalRecHitCollection> ebRecHitsToken_;
  edm::EDGetTokenT<EcalRecHitCollection> eeRecHitsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPartsToken_;
  int nrBarrel_;
  int nrEndcap_;
  TH2* energyMapBarrel_;
  TH2* energyMapEndcap_;

  LocalEnergyMapProducer(const LocalEnergyMapProducer& rhs)=delete;
  LocalEnergyMapProducer& operator=(const LocalEnergyMapProducer& rhs)=delete;

public:
  explicit LocalEnergyMapProducer(const edm::ParameterSet& iPara);
  virtual ~LocalEnergyMapProducer();
  
private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& iSetup);
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endRun(edm::Run const& iRun, edm::EventSetup const&){}
  virtual void endJob();
  template<typename T>
  void setToken(edm::EDGetTokenT<T>& token,const edm::ParameterSet& iPara,const std::string& tagName){
    token = consumes<T>(iPara.getParameter<edm::InputTag>(tagName));
  }
  template<typename T>
  void setToken(std::vector<edm::EDGetTokenT<T> > & tokens,const edm::ParameterSet& iPara,const std::string& tagName){
    for(auto& tag: iPara.getParameter<std::vector<edm::InputTag> >(tagName)){
      tokens.push_back(consumes<T>(tag));
    }
  }
  void fillLocalEnergyMap(DetId seedId,const EcalRecHitCollection& ebHits,const EcalRecHitCollection& eeHits);
  void fillLocalEnergyMapBarrel(EBDetId seedId,const EcalRecHitCollection& hits);
  void fillLocalEnergyMapEndcap(EEDetId seedId,const EcalRecHitCollection& hits);

  static float getEnergy(DetId hitId,const EcalRecHitCollection& hits);
  template<typename DetIdType>
  static float calE5x5(DetIdType seedId,const EcalRecHitCollection& hits);
  static std::vector<const reco::GenParticle*> filterGenParts(edm::Handle<reco::GenParticleCollection> genParts,const std::vector<int>& allowedPIDs);
  static bool isBestGenPartMatch(const reco::SuperCluster& sc,const std::vector<const reco::GenParticle*> genParts,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles);
};



LocalEnergyMapProducer::LocalEnergyMapProducer(const edm::ParameterSet& iPara):
  fillFromEles_(iPara.getParameter<bool>("fillFromEles")),
  fillFromPhos_(iPara.getParameter<bool>("fillFromPhos")),
  fillFromSCs_(iPara.getParameter<bool>("fillFromSCs")),
  genMatch_(iPara.getParameter<bool>("genMatch")),
  nrBarrel_(0),nrEndcap_(0)
{
  setToken(elesToken_,iPara,"elesTag");
  setToken(phosToken_,iPara,"phosTag");
  setToken(scTokens_,iPara,"scTags");
  setToken(ebRecHitsToken_,iPara,"ebRecHitsTag");
  setToken(eeRecHitsToken_,iPara,"eeRecHitsTag"); 
  setToken(genPartsToken_,iPara,"genPartsTag");
}

LocalEnergyMapProducer::~LocalEnergyMapProducer()
{

}


void LocalEnergyMapProducer::beginJob()
{
  edm::Service<TFileService> fs;
  fs->file().cd();
  energyMapBarrel_ = new TH2D("energyMapBarrel",";local #ieta;local #iphi",5,-2.5,2.5,5,-2.5,2.5);
  energyMapEndcap_ = new TH2D("energyMapEndcap",";local ix;local iy",5,-2.5,2.5,5,-2.5,2.5);
  energyMapBarrel_->SetDirectory(&fs->file());
  energyMapEndcap_->SetDirectory(&fs->file());
} 

void LocalEnergyMapProducer::beginRun(const edm::Run& run,const edm::EventSetup& iSetup)
{ 
 
}

namespace {
  template<typename T> 
  edm::Handle<T> getHandle(const edm::Event& iEvent,const edm::EDGetTokenT<T>& token)
  {
    edm::Handle<T> handle;
    iEvent.getByToken(token,handle);
    return handle;
  }
  template<typename T> 
  std::vector<edm::Handle<T> > getHandle(const edm::Event& iEvent,const std::vector<edm::EDGetTokenT<T> >& tokens)
  {
    std::vector<edm::Handle<T> > handles;
    for(auto& token : tokens){
      edm::Handle<T> handle;
      iEvent.getByToken(token,handle);
      handles.emplace_back(std::move(handle));
    }
    return handles;
  }
}


void LocalEnergyMapProducer::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup)
{
  auto elesHandle = getHandle(iEvent,elesToken_);
  auto phosHandle = getHandle(iEvent,phosToken_);
  auto ebRecHitsHandle = getHandle(iEvent,ebRecHitsToken_);
  auto eeRecHitsHandle = getHandle(iEvent,eeRecHitsToken_);
  auto scHandles = getHandle(iEvent,scTokens_);
  auto genParts = getHandle(iEvent,genPartsToken_);

  std::vector<const reco::GenParticle*> monoPoles = filterGenParts(genParts,{-4110000,4110000});

  if(fillFromEles_){
    for(const reco::GsfElectron& ele : *elesHandle){
      fillLocalEnergyMap(ele.superCluster()->seed()->seed(),
			 *ebRecHitsHandle,*eeRecHitsHandle);
    }
  }else if(fillFromPhos_){
    for(const reco::Photon& pho : *phosHandle){
      //put some selection on the photon
      fillLocalEnergyMap(pho.superCluster()->seed()->seed(),
			 *ebRecHitsHandle,*eeRecHitsHandle);
    }
  }else if(fillFromSCs_){
    for(auto& scHandle : scHandles){
      for(const reco::SuperCluster& sc : *scHandle){
	std::cout <<" sc filling "<<std::endl;
	if(!genMatch_ || isBestGenPartMatch(sc,monoPoles,scHandles)){
	  fillLocalEnergyMap(sc.seed()->seed(),*ebRecHitsHandle,*eeRecHitsHandle);
	}
      }
    }
  }
}
		

void LocalEnergyMapProducer::endJob()
{ 
  if(nrBarrel_!=0) energyMapBarrel_->Scale(1./nrBarrel_);
  if(nrEndcap_!=0) energyMapEndcap_->Scale(1./nrEndcap_);

}

void LocalEnergyMapProducer::fillLocalEnergyMap(DetId seedId,const EcalRecHitCollection& ebHits,const EcalRecHitCollection& eeHits)
{
  if(seedId.subdetId()==EcalBarrel) fillLocalEnergyMapBarrel(seedId,ebHits);
  else if(seedId.subdetId()==EcalEndcap) fillLocalEnergyMapEndcap(seedId,eeHits);
  else throw cms::Exception("LogicError") << "error sub det "<<seedId.subdetId()<<" is neither barrel nor endcap"<<std::endl;
}


float LocalEnergyMapProducer::getEnergy(DetId hitId,const EcalRecHitCollection& hits)
{
  if(hitId!=DetId(0)){
    EcalRecHitCollection::const_iterator it = hits.find(hitId);
    if(it!=hits.end()) return it->energy();
  }
  //didnt exist or not valid id
  return 0.;
}

template<typename DetIdType>
float LocalEnergyMapProducer::calE5x5(DetIdType seedId,const EcalRecHitCollection& hits)
{
  float e5x5=0.;
  for(int iEtaOrXNr=-2;iEtaOrXNr<=2;iEtaOrXNr++){
    for(int iPhiOrYNr=-2;iPhiOrYNr<=2;iPhiOrYNr++){
      DetId id = seedId.offsetBy(iEtaOrXNr,iPhiOrYNr);
      float energy = getEnergy(id,hits);
      e5x5+=energy;
    }
  }
  return e5x5;
}


bool LocalEnergyMapProducer::isBestGenPartMatch(const reco::SuperCluster& sc,const std::vector<const reco::GenParticle*> genParts,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles)
{
  float maxDR2 = 0.2*0.2;
  const reco::GenParticle* bestGenMatch = nullptr;
  for(auto& genPart : genParts){
    float dR2 = reco::deltaR2(sc.eta(),sc.phi(),genPart->eta(),genPart->phi());
    if(dR2<maxDR2){
      bestGenMatch = genPart;
      maxDR2 = dR2;
    }
  }
  
  if(bestGenMatch){
    for(auto& scHandle : scHandles){
      for(auto& otherSC : *scHandle){
	if(&otherSC != &sc){
	  float dR2 = reco::deltaR2(otherSC.eta(),otherSC.phi(),bestGenMatch->eta(),bestGenMatch->phi());
	  if(dR2<maxDR2) return false; //other SC matches better
	}
      }
    }
    return true; //match found and no better match found
  }else{
    return false; //no match found
  }
}


void LocalEnergyMapProducer::fillLocalEnergyMapBarrel(EBDetId seedId,const EcalRecHitCollection& hits)
{
  float e5x5 = calE5x5(seedId,hits); 
  if(e5x5==0) return;
  
  nrBarrel_++;

  for(int iEtaNr=-2;iEtaNr<=2;iEtaNr++){
    for(int iPhiNr=-2;iPhiNr<=2;iPhiNr++){
      DetId id = seedId.offsetBy(iEtaNr,iPhiNr);
      float energy = getEnergy(id,hits);
      energyMapBarrel_->Fill(iEtaNr,iPhiNr,energy/e5x5);
    }
  }
}
	
void LocalEnergyMapProducer::fillLocalEnergyMapEndcap(EEDetId seedId,const EcalRecHitCollection& hits)
{
  float e5x5 = calE5x5(seedId,hits); 
  if(e5x5==0) return;
  
  nrEndcap_++;

  for(int iXNr=-2;iXNr<=2;iXNr++){
    for(int iYNr=-2;iYNr<=2;iYNr++){
      DetId id = seedId.offsetBy(iXNr,iYNr);
      float energy = getEnergy(id,hits);
      energyMapEndcap_->Fill(iXNr,iYNr,energy/e5x5);
    }
  }
}
std::vector<const reco::GenParticle*> LocalEnergyMapProducer::filterGenParts(edm::Handle<reco::GenParticleCollection> genParts,const std::vector<int>& allowedPIDs)
{
  std::vector<const reco::GenParticle*> filteredParts;
  for(auto& genPart: * genParts){
    if(genPart.isPromptFinalState() && 
       std::find(allowedPIDs.begin(),allowedPIDs.end(),genPart.pdgId())!=allowedPIDs.end()){
      filteredParts.push_back(&genPart);
    }
  }
  return filteredParts;
}


//define this as a plug-in
DEFINE_FWK_MODULE(LocalEnergyMapProducer);
  

