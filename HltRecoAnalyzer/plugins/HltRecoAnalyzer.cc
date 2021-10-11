// -*- C++ -*-
//
// Package:    DeepTauStudy/HltRecoAnalyzer
// Class:      HltRecoAnalyzer
//
/**\class HltRecoAnalyzer HltRecoAnalyzer.cc DeepTauStudy/HltRecoAnalyzer/plugins/HltRecoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Anna Mascellani
//         Created:  Mon, 17 May 2021 14:32:53 GMT
//
//

// system include files
#include <memory>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TNamed.h"
#include <vector>
#include <map>
#include "TROOT.h"
#include "TChain.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "TauTriggerTools/Common/interface/AnalysisTypes.h"
#include "TauTriggerTools/Common/interface/GenTruthTools.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class HltRecoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HltRecoAnalyzer(const edm::ParameterSet&);

  // Analysis tree to be filled
    TTree *HltTree;

private:
  // void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  int run, lumi, evt;
  int trigflag;
  std::vector<float> trig_tau_pt; /* pt of the hlt tau */ 
  std::vector<float> trig_tau_eta; /* eta of the hlt tau */ 
  std::vector<float> trig_tau_phi; /* phi of the hlt tau */ 
  std::vector<float> trig_tau_e; /* energy of the hlt tau */
  std::vector<float> off_tau_pt; /* pt of the offline tau */ 
  std::vector<float> off_tau_eta; /* eta of the offline tau */ 
  std::vector<float> off_tau_phi; /* phi of the offline tau */ 
  std::vector<float> off_tau_e; /* energy of the offline tau */ 
  std::vector<float> off_deepTau_VSe; /* off deepTau vs e */
  std::vector<float> off_deepTau_VSmu; /* off deepTau vs mu */
  std::vector<float> off_deepTau_VSjet; /* off deepTau vs jet */
  std::vector<bool> off_decayModeFinding; /* new decayModeFinding flag */
  std::vector<int> off_decayMode; /* decay mode of the offline tau */
  std::vector<float> off_tau_dz; /* dz of the offline tau */
  std::vector<float> off_gen_tau_pt; /* pt of the gen tau */ 
  std::vector<float> off_gen_tau_eta; /* eta of the gen tau */ 
  std::vector<float> off_gen_tau_phi; /* phi of the gen tau */ 
  std::vector<float> off_gen_tau_e; /* energy of the gen tau */ 
  std::vector<int> off_lepton_gen_match; /* lepton gen match of the offline tau */


  edm::EDGetTokenT<edm::TriggerResults> hltresults_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > trig_objs_;
  edm::EDGetTokenT<std::vector<pat::Tau>> taus_offline_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> gen_particles_;
  const std::vector<reco::GenParticle>* genParticles;

  TFile* m_file; // pointer to Histogram file
  std::string _HistName; // Name of Histogram file
  std::string trigName_chosen;
  std::string filter_to_pass; // Name of last HLT filter to be passed by hlt objects
};

//
// constructors and destructor
//
HltRecoAnalyzer::HltRecoAnalyzer(const edm::ParameterSet& iConfig): 
  hltresults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltresults"))),
  trig_objs_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("trig_objects"))),
  taus_offline_(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus_offline"))),
  gen_particles_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen_particles"))),
  trigName_chosen(iConfig.getParameter<std::string>("hlt_path")),
  filter_to_pass(iConfig.getParameter<std::string>("filter_to_pass"))
  {

  // open the tree file
  edm::ParameterSet runParameters = iConfig.getParameter<edm::ParameterSet>("run_parameters");
  _HistName = runParameters.getUntrackedParameter<std::string>("HistogramFile", "test.root");
  m_file = new TFile(_HistName.c_str(), "RECREATE");
  if (m_file)
      m_file->cd();
  
  // Initialize the tree
  HltTree = new TTree("HltTree", "");

  TString trigName = trigName_chosen;

  HltTree->Branch("run", &run);
  HltTree->Branch("evt", &evt);
  HltTree->Branch("lumi", &lumi);
  HltTree->Branch(trigName, &trigflag);
  HltTree->Branch("tau_trig_pt", &trig_tau_pt); 
  HltTree->Branch("tau_trig_eta", &trig_tau_eta); 
  HltTree->Branch("tau_trig_phi", &trig_tau_phi);
  HltTree->Branch("tau_trig_e", &trig_tau_e); 
  HltTree->Branch("tau_off_pt", &off_tau_pt); 
  HltTree->Branch("tau_off_eta", &off_tau_eta); 
  HltTree->Branch("tau_off_phi", &off_tau_phi); 
  HltTree->Branch("tau_off_e", &off_tau_e); 
  HltTree->Branch("tau_off_deepTau_VSe", &off_deepTau_VSe);
  HltTree->Branch("tau_off_deepTau_VSmu", &off_deepTau_VSmu);
  HltTree->Branch("tau_off_deepTau_VSjet", &off_deepTau_VSjet);
  HltTree->Branch("tau_off_decayModeFinding", &off_decayModeFinding);
  HltTree->Branch("tau_off_decayMode", &off_decayMode);
  HltTree->Branch("tau_off_dz", &off_tau_dz);
  HltTree->Branch("gen_tau_off_pt", &off_gen_tau_pt);
  HltTree->Branch("gen_tau_off_eta", &off_gen_tau_eta);
  HltTree->Branch("gen_tau_off_phi", &off_gen_tau_phi);
  HltTree->Branch("gen_tau_off_e", &off_gen_tau_e);
  HltTree->Branch("tau_off_lepton_gen_match", &off_lepton_gen_match);
}

//
// member functions
//

// ------------ method called for each event  ------------
void HltRecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  evt = iEvent.id().event();
  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();
  std::cout << evt << "\t" << run << "\t" << lumi << std::endl;

  off_tau_pt.clear();
  off_tau_eta.clear();
  off_tau_phi.clear();
  off_tau_e.clear();
  off_deepTau_VSe.clear();
  off_deepTau_VSmu.clear();
  off_deepTau_VSjet.clear();
  off_decayModeFinding.clear();
  off_decayMode.clear();
  off_tau_dz.clear();
  off_gen_tau_pt.clear();
  off_gen_tau_eta.clear();
  off_gen_tau_phi.clear();
  off_gen_tau_e.clear();
  off_lepton_gen_match.clear();

  edm::Handle<edm::TriggerResults> hltresults;
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trig_objs;
  edm::Handle<std::vector<pat::Tau>> taus_offline;
  edm::Handle<std::vector<reco::GenParticle>> gen_particles;

  // iEvent.getByLabel(hltresults_, hltresults);
  iEvent.getByToken(hltresults_, hltresults);
  iEvent.getByToken(trig_objs_, trig_objs);
  // std::cout << "trig objs handle is valid: " << trig_objs.isValid() << std::endl;
  iEvent.getByToken(taus_offline_, taus_offline);
  iEvent.getByToken(gen_particles_, gen_particles);

  genParticles = gen_particles.isValid() ? gen_particles.product() : nullptr;
  float default_value = -999.;
  int default_int_value = -999;

  // std::cout << hltresults.isValid() << std::endl;

  if (hltresults.isValid()) {
    int ntrigs = hltresults->size();
    if (ntrigs==0){std::cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;}

    edm::TriggerNames const& triggerNames = iEvent.triggerNames(*hltresults);
    bool found = false;
    for (int itrig = 0; itrig != ntrigs; ++itrig){
      std::string trigName=triggerNames.triggerName(itrig);
      if (trigName==trigName_chosen){
        found = true;
        bool accept = hltresults->accept(itrig);
        if (accept){trigflag = 1;}
        else {trigflag = 0;}
      }
    }
    std::cout << "Trigger flag: " << trigflag << std::endl;
    if (!found){std::cout << "%HLTInfo -- No trigger result found for " << trigName_chosen << std::endl;}

  }

  // now we look whether the hlt filter is passed
  // before we do this we need to make a new collection of trig objects
  // with their filters unpacked
  // we have to make a new copy as the unpacking modifies them and CMSSW
  // forbids (for very good reasons) modification of products in the event
  std::vector<pat::TriggerObjectStandAlone> trig_objs_unpacked;
  for(auto& trigObj : *trig_objs){
    trig_objs_unpacked.push_back(trigObj);
    trig_objs_unpacked.back().unpackFilterLabels(iEvent,*hltresults);
  }

  //now we actually check whether the hlt filter is passed
  for(const auto trigObj : trig_objs_unpacked){
    bool passed_filter = false;
    const std::vector<std::string>& objFilters = trigObj.filterLabels();
    int n_trig_filters = objFilters.size();
    for(int i=0; i<n_trig_filters; i++){
      if(objFilters.at(i)==filter_to_pass) {passed_filter = true;}
    }
    // and store hlt object p4
    if(passed_filter && trigObj.type(trigger::TriggerObjectType::TriggerTau)) {
      std::cout <<std::endl <<"Tau passing filter! " <<std::endl;
      trig_tau_pt.push_back(trigObj.polarP4().pt());
      std::cout << trig_tau_pt.back() << std::endl;
      trig_tau_eta.push_back(trigObj.polarP4().eta());
      std::cout << trig_tau_eta.back() << std::endl;
      trig_tau_phi.push_back(trigObj.polarP4().phi());
      std::cout << trig_tau_phi.back() << std::endl;
      trig_tau_e.push_back(trigObj.polarP4().e());
      std::cout << trig_tau_e.back() << std::endl;
      }
  }//end loop over trigger objects

  // Store offline taus info
  if(taus_offline.isValid()) {
    for(size_t off_tau_index = 0; off_tau_index < taus_offline->size(); ++off_tau_index) {
      const pat::Tau& offline_tau = taus_offline->at(off_tau_index);

      if(genParticles) {
        const auto off_gen_match = analysis::gen_truth::LeptonGenMatch(offline_tau.polarP4(), *genParticles);
        off_lepton_gen_match.push_back(static_cast<int>(off_gen_match.match));
        off_gen_tau_pt.push_back(static_cast<float>(off_gen_match.visible_p4.pt()));
        off_gen_tau_eta.push_back(static_cast<float>(off_gen_match.visible_p4.eta()));
        off_gen_tau_phi.push_back(static_cast<float>(off_gen_match.visible_p4.phi()));
        off_gen_tau_e.push_back(static_cast<float>(off_gen_match.visible_p4.e()));
      } else {
        off_lepton_gen_match.push_back(default_int_value);
        off_gen_tau_pt.push_back(default_value);
        off_gen_tau_eta.push_back(default_value);
        off_gen_tau_phi.push_back(default_value);
        off_gen_tau_e.push_back(default_value);
      }

      off_tau_pt.push_back(offline_tau.polarP4().pt());
      // std::cout << offline_tau.polarP4().pt() << std::endl;
      off_tau_eta.push_back(offline_tau.polarP4().eta());
      off_tau_phi.push_back(offline_tau.polarP4().phi());
      off_tau_e.push_back(offline_tau.polarP4().e());
      off_deepTau_VSe.push_back(offline_tau.tauID("byDeepTau2017v2p1VSeraw"));
      off_deepTau_VSmu.push_back(offline_tau.tauID("byDeepTau2017v2p1VSmuraw"));
      off_deepTau_VSjet.push_back(offline_tau.tauID("byDeepTau2017v2p1VSjetraw"));
      off_decayModeFinding.push_back(static_cast<bool>(offline_tau.tauID("decayModeFindingNewDMs")));
      off_decayMode.push_back(static_cast<int>(offline_tau.decayMode()));
      const pat::PackedCandidate* leadChargedHadrCand = dynamic_cast<const pat::PackedCandidate*>(offline_tau.leadChargedHadrCand().get());
      off_tau_dz.push_back(static_cast<float>(leadChargedHadrCand->dz()));
    }
  }
  
  // std::cout << " Ending Event Analysis" << std::endl;
  // After analysis, fill the variables tree
  if (m_file)
    m_file->cd();
  HltTree->Fill();
}

// // ------------ method called once each job just before starting event loop  ------------
// void HltRecoAnalyzer::beginJob() {
//   // please remove this method if not needed
// }

// ------------ method called once each job just after ending the event loop  ------------
void HltRecoAnalyzer::endJob() {
  if (m_file)
      m_file->cd();
  HltTree->Write();
  delete HltTree;
  HltTree = 0;
  
  if (m_file) {         // if there was a tree file...
    m_file->Write();    // write out the branches
    delete m_file;      // close and delete the file
    m_file = 0;         // set to zero to clean up
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HltRecoAnalyzer);