// -*- C++ -*-
//
// Package: RootupleChib
// Class: RootupleChib
//
/**\class RootupleChib RootupleChib.cc RootupleChib/RootupleChib/src/RootupleChib.cc
 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author: Alessandro Degano,32 1-C13,+41227678098
// Created: Tue Feb 19 14:32:37 CET 2013
// $Id: RootupleChib.cc,v 1.4 2013/06/21 12:49:49 gdujany Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


//
// class declaration
//

class RootupleChib : public edm::EDAnalyzer {
public:
    explicit RootupleChib(const edm::ParameterSet&);
    ~RootupleChib();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    
    // ----------member data ---------------------------
    std::string file_name;
    edm::InputTag chi_cand_Label;
    edm::InputTag pi0_comb_Label;
    edm::InputTag ups_Label;
    edm::InputTag refit1S_Label;
    edm::InputTag refit2S_Label;
    edm::InputTag refit3S_Label;
    edm::InputTag primaryVertices_Label;
    edm::InputTag triggerResults_Label;
    bool isMC_;
    
    UInt_t run;
    UInt_t event;

	TLorentzVector chi_p4;
    TLorentzVector dimuon_p4;
    TLorentzVector muonP_p4;
    TLorentzVector muonN_p4;
    TLorentzVector photon_p4;
    
    TLorentzVector rf1S_chi_p4;
    TLorentzVector rf1S_dimuon_p4;
    TLorentzVector rf1S_muonP_p4;
    TLorentzVector rf1S_muonN_p4;
    TLorentzVector rf1S_photon_p4;
    
    TLorentzVector rf2S_chi_p4;
    TLorentzVector rf2S_dimuon_p4;
    TLorentzVector rf2S_muonP_p4;
    TLorentzVector rf2S_muonN_p4;
    TLorentzVector rf2S_photon_p4;

    /*
    Double_t invm1S;
    Double_t invm2S;
    Double_t invm3S;
    Double_t chib_mass;
    Double_t chib_pt;
    Double_t chib_eta;
    Double_t chib_phi;
    Double_t dimuon_mass;
    Double_t dimuon_rapidity;
    Double_t dimuon_pt;
    Double_t photon_eta;
    Double_t photon_pt;
    */
    
    Double_t ele_lowerPt_pt;
    Double_t ele_higherPt_pt;
    Double_t ctpv;
    Double_t ctpv_error;
    Double_t pi0_abs_mass;
    Double_t Y1S_nsigma;
    Double_t Y2S_nsigma;
    Double_t Y3S_nsigma;
    Double_t conv_vertex;
    Double_t dz;
    Double_t numPrimaryVertices;
    Int_t trigger;

    //Double_t rf1S_chib_mass;
    Double_t probFit1S;
    Int_t rf1S_rank;
    /*
    Double_t rf1S_chib_pt;
    Double_t rf1S_chib_eta;
    Double_t rf1S_dimuon_mass;
    Double_t rf1S_dimuon_rapidity;
    Double_t rf1S_dimuon_pt;
    Double_t rf1S_photon_eta;
    Double_t rf1S_photon_pt;
    */
    
    //Double_t rf2S_chib_mass;
    Double_t probFit2S;
    //Double_t rf3S_chib_mass;
    Double_t probFit3S;
    
    /*
    Double_t ups_mass;
    Double_t ups_rapidity;
    Double_t ups_pt;
    */

    TTree* chib_tree;
    //TTree* upsilon_tree;
    //TLorentzVector lorVect;
    
    /*
    TLorentzVector chib_p4;
    Int_t chib_pdgId;
    TLorentzVector Upsilon_p4;
    Int_t Upsilon_pdgId;
    TLorentzVector photon_p4;
    TLorentzVector muP_p4;
    TLorentzVector muM_p4;
    */
    
    TLorentzVector gen_chi_p4;
    Int_t gen_chi_pdgId;
    TLorentzVector gen_dimuon_p4;
    Int_t gen_dimuon_pdgId;
    TLorentzVector gen_photon_p4;
    TLorentzVector gen_muP_p4;
    TLorentzVector gen_muM_p4;
    

    
};




//
// constants, enums and typedefs
//

//
// static data member definitions
//
static const double pi0_mass = 0.1349766;
static const Double_t Y1SMass = 9.46030;
static const Double_t Y2SMass = 10.02326;
static const Double_t Y3SMass = 10.3552;

// 2011 par
//static const double Y_sig_par_A = 0.058;
//static const double Y_sig_par_B = 0.047;
//static const double Y_sig_par_C = 0.22;

// 2012 par
static const double Y_sig_par_A = 62.62;
static const double Y_sig_par_B = 56.3;
static const double Y_sig_par_C = -20.77;

//
// constructors and destructor
//
RootupleChib::RootupleChib(const edm::ParameterSet& iConfig):
chi_cand_Label(iConfig.getParameter<edm::InputTag>("chi_cand")),
pi0_comb_Label(iConfig.getParameter<edm::InputTag>("pi0_comb")),
ups_Label(iConfig.getParameter<edm::InputTag>("ups_cand")),
refit1S_Label(iConfig.getParameter<edm::InputTag>("refit1S")),
refit2S_Label(iConfig.getParameter<edm::InputTag>("refit2S")),
refit3S_Label(iConfig.getParameter<edm::InputTag>("refit3S")),
primaryVertices_Label(iConfig.getParameter<edm::InputTag>("primaryVertices")),
triggerResults_Label(iConfig.getParameter<edm::InputTag>("TriggerResults")),
isMC_(iConfig.getParameter<bool>("isMC"))
{
    edm::Service<TFileService> fs;
    chib_tree = fs->make<TTree>("chibTree","Tree of chib");
    
    
    chib_tree->Branch("run" , &run, "run/I");
    chib_tree->Branch("event", &event, "event/I");
    




