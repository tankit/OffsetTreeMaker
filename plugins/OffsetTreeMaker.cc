// -*- C++ -*-
//
// Package:    treemaker/OffsetTreeMaker
// Class:      OffsetTreeMaker
//
/**\class OffsetTreeMaker OffsetTreeMaker.cc treemaker/OffsetTreeMaker/plugins/OffsetTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  charles harrington
//         Created:  Mon, 09 Nov 2015 17:09:43 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Common/interface/Ref.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "parsePileUpJSON2.h"
#include <vector>
#include "TMath.h"
//root files
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>

using namespace std;

const int ETA_BINS = 82;
const int ETA_BINS_GME = 18;
const int PHI_BINS_GME = 11;
const double GRID_AREA = (10./ETA_BINS_GME)*(2*M_PI/PHI_BINS_GME);
const int MAXNPV = 150;
const int MAXJETS = 4;


float muWeight_Run2017D[80] = { 0.000000,0.000000,0.000002,0.000019,0.000066,0.000127,0.000154,0.000164,0.000159,0.000164,0.000201,0.000251,0.000383,0.000707,0.001715,0.005040,0.012966,0.023951,0.033498,0.040964,0.046665,0.049625,0.051241,0.054287,0.059305,0.064711,0.068642,0.070111,0.068788,0.064589,0.057980,0.049998,0.041697,0.033764,0.026589,0.020402,0.015300,0.011246,0.008106,0.005715,0.003926,0.002617,0.001687,0.001050,0.000631,0.000366,0.000206,0.000112,0.000058,0.000030,0.000014,0.000007,0.000003,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000 };
float muWeight_Run2017E[80] = { 0.000000,0.000001,0.000010,0.000019,0.000041,0.000084,0.000103,0.000134,0.000142,0.000149,0.000179,0.000347,0.000791,0.001573,0.002727,0.004371,0.006719,0.009951,0.014011,0.018429,0.022376,0.025295,0.027196,0.028305,0.028876,0.029197,0.029491,0.029891,0.030495,0.031348,0.032376,0.033378,0.034091,0.034331,0.034082,0.033479,0.032694,0.031851,0.031002,0.030153,0.029294,0.028422,0.027538,0.026641,0.025709,0.024696,0.023538,0.022170,0.020548,0.018665,0.016562,0.014324,0.012059,0.009878,0.007876,0.006118,0.004637,0.003436,0.002494,0.001777,0.001247,0.000863,0.000591,0.000401,0.000271,0.000182,0.000123,0.000082,0.000055,0.000037,0.000025,0.000017,0.000012,0.000008,0.000005,0.000004,0.000002,0.000002,0.000001,0.000001 };

float etabins[ETA_BINS+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

float phibins[PHI_BINS_GME+1] = 
  {-3.142, -2.57, -1.999, -1.428, -0.8568, -0.2856, 0.2856, 0.8568, 1.428, 1.999, 2.57, 3.142};
class OffsetTreeMaker : public edm::EDAnalyzer {
  public:
    explicit OffsetTreeMaker(const edm::ParameterSet&);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    int getEtaIndex(float eta);
    enum Flavor{
      chm = 0, chu, nh, ne, hfh, hfe, lep, untrk, numFlavors, X //undefined
    };
    Flavor getFlavor(int id);

    int counter;
    TFile* root_file;
    TTree* tree;

    TH1F* h;
    TH2F* h2_GME; // 2d ET_eta_phi histo 
    TH2F* h2_finnereta;
//    TH2F *h2_GME_gen, *h2_finnereta_gen; // for GenParticle
    TH2F *h2_twopi, *h2_twopi_gen, *h2_twopi_det, *h2_twopi_scaled; // for phi from 0 to 2*Pi
    TRandom3* rand;

    int nEta;
    float energy[ETA_BINS], eRMS[ETA_BINS], et[ETA_BINS], etMED[ETA_BINS], etMEAN[ETA_BINS];
    float et_gme[ETA_BINS_GME][PHI_BINS_GME];
//    float et_gme_gen[ETA_BINS_GME][PHI_BINS_GME];
    UChar_t f[numFlavors][ETA_BINS];  //energy fraction by flavor
//    float et_gen[ETA_BINS], etMED_gen[ETA_BINS], etMEAN_gen[ETA_BINS];
    float et_twopi[ETA_BINS_GME][PHI_BINS_GME], et_twopi_gen[ETA_BINS_GME][PHI_BINS_GME], et_twopi_det[ETA_BINS_GME][PHI_BINS_GME], et_twopi_scaled[ETA_BINS_GME][PHI_BINS_GME];
    float rho_gme, rho_gme_gen, rho_gme_det, rho_gme_scaled;

    ULong64_t event;
    int run, lumi, bx;
    float mu, muWeight;
    float rho, rhoC0, rhoCC;
    float rhoCentral, rhoCentralCalo;
    float weight;
    int nGenParticles;
    int idHardScat1, idHardScat2, idxHardScat;

    int nPVall, nPV;
    float pv_ndof[MAXNPV], pv_z[MAXNPV], pv_rho[MAXNPV];

    float ht;
    int nJets;
    float jet_eta[MAXJETS], jet_phi[MAXJETS], jet_pt[MAXJETS], jet_area[MAXJETS];
    float jet_ch[MAXJETS], jet_nh[MAXJETS], jet_ne[MAXJETS], jet_hfh[MAXJETS], jet_hfe[MAXJETS], jet_lep[MAXJETS];

    vector<int> pf_type;
    vector<float> pf_pt, pf_eta, pf_phi, pf_et;
    vector<int> particle_id, particle_charge;
    vector<float> particle_et, particle_pt, particle_eta, particle_phi;
    vector<float> particle_et_scaled;
    vector<bool> isPromptFinalState;

    TString RootFileName_;
    string puFileName_;
    int numSkip_;
    bool isMC_, writeCands_, writeParticles_;

    edm::EDGetTokenT< vector<reco::Vertex> > pvTag_;
//    edm::EDGetTokenT< vector<reco::Track> > trackTag_;
    edm::EDGetTokenT< vector<PileupSummaryInfo> > muTag_;
    edm::EDGetTokenT< vector<pat::PackedCandidate> > pfTag_;
    edm::EDGetTokenT< vector<pat::PackedGenParticle> > genTag_;
    edm::EDGetTokenT< vector<reco::GenParticle> > genparticlesTag_;
    edm::EDGetTokenT< GenEventInfoProduct > generatorTag_;
    edm::EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT<double> rhoC0Tag_;
    edm::EDGetTokenT<double> rhoCCTag_;
    edm::EDGetTokenT<double> rhoCentralTag_;
    edm::EDGetTokenT<double> rhoCentralCaloTag_;
    edm::EDGetTokenT< vector<pat::Jet> > pfJetTag_;
};

OffsetTreeMaker::OffsetTreeMaker(const edm::ParameterSet& iConfig)
{
  numSkip_ = iConfig.getParameter<int> ("numSkip");
  RootFileName_ = iConfig.getParameter<string>("RootFileName");
  puFileName_ = iConfig.getParameter<string>("puFileName");
  isMC_ = iConfig.getParameter<bool>("isMC");
  writeCands_ = iConfig.getParameter<bool>("writeCands");
  writeParticles_ = iConfig.getParameter<bool>("writeParticles");
  pvTag_ = consumes< vector<reco::Vertex> >( iConfig.getParameter<edm::InputTag>("pvTag") );
//  trackTag_ = consumes< vector<reco::Track> >( iConfig.getParameter<edm::InputTag>("trackTag") );
  muTag_ = consumes< vector<PileupSummaryInfo> >( iConfig.getParameter<edm::InputTag>("muTag") );
  pfTag_ = consumes< vector<pat::PackedCandidate> >( iConfig.getParameter<edm::InputTag>("pfTag") );
  genTag_ = consumes< vector<pat::PackedGenParticle> >( iConfig.getParameter<edm::InputTag>("genTag") );
  genparticlesTag_ = consumes< vector<reco::GenParticle> >( iConfig.getParameter<edm::InputTag>("GenParticles") );
  generatorTag_ = consumes<GenEventInfoProduct>( iConfig.getParameter<edm::InputTag>("Generator") );
  rhoTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoTag") );
  rhoC0Tag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoC0Tag") );
  rhoCCTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoCCTag") );
  rhoCentralTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoCentralTag") );
  rhoCentralCaloTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoCentralCaloTag") );
  pfJetTag_ = consumes< vector<pat::Jet> >( iConfig.getParameter<edm::InputTag>("pfJetTag") );
}

// ------------ method called once each job just before starting event loop  ------------
void  OffsetTreeMaker::beginJob() {

  root_file = new TFile(RootFileName_,"RECREATE");
  tree = new TTree("T","Offset Tree");

  counter = -1;
  h = new TH1F("mu", "mu", 100, 0, 50);
  rand = new TRandom3;
  h2_GME = new TH2F("GME_2D", "GME_2D_yPhi_xEta", ETA_BINS_GME, -5.0, 5.0, PHI_BINS_GME , -1*M_PI , M_PI);
  h2_finnereta = new TH2F("finnereta_ET", "yPhi_xEta", ETA_BINS, etabins, PHI_BINS_GME, phibins);
//  h2_GME_gen = new TH2F("GME_2D_gen", "GME_2D_yPhi_xEta", ETA_BINS_GME, -5.0, 5.0, PHI_BINS_GME , -1*M_PI , M_PI);
//  h2_finnereta_gen = new TH2F("finnereta_ET_gen", "yPhi_xEta", ETA_BINS, etabins, PHI_BINS_GME, phibins);
  h2_twopi = new TH2F("GME_twopi", "GME_2D_yPhi_xEta", ETA_BINS_GME, -5.0, 5.0, PHI_BINS_GME , 0 , 2*M_PI);
  h2_twopi_gen = new TH2F("GME_twopi_gen", "GME_2D_yPhi_xEta", ETA_BINS_GME, -5.0, 5.0, PHI_BINS_GME , 0 , 2*M_PI);
  h2_twopi_det = new TH2F("GME_twopi_det", "GME_2D_yPhi_xEta", ETA_BINS_GME, -5.0, 5.0, PHI_BINS_GME , 0 , 2*M_PI);
  h2_twopi_scaled = new TH2F("GME_twopi_scaled", "GME_2D_yPhi_xEta", ETA_BINS_GME, -5.0, 5.0, PHI_BINS_GME , 0 , 2*M_PI);

  if (!isMC_){
    parsePileUpJSON2( puFileName_ );

    tree->Branch("run", &run, "run/I");
    tree->Branch("lumi", &lumi, "lumi/I");
    tree->Branch("bx", &bx, "bx/I");
    tree->Branch("event", &event, "event/l");
  }

  if (writeCands_) {
    tree->Branch("pf_type", "std::vector<int>",   &pf_type);
    tree->Branch("pf_pt",   "std::vector<float>", &pf_pt);
    tree->Branch("pf_eta",  "std::vector<float>", &pf_eta);
    tree->Branch("pf_phi",  "std::vector<float>", &pf_phi);
    tree->Branch("pf_et",  "std::vector<float>", &pf_et);
  }

  tree->Branch("mu", &mu, "mu/F");
  tree->Branch("muWeight", &muWeight, "muWeight/F");

  tree->Branch("rho",   &rho,   "rho/F");
  tree->Branch("rhoC0", &rhoC0, "rhoC0/F");
  tree->Branch("rhoCC", &rhoCC, "rhoCC/F");
  tree->Branch("rhoCentral", &rhoCentral, "rhoCentral/F");
  tree->Branch("rhoCentralCalo", &rhoCentralCalo, "rhoCentralCalo/F");

  tree->Branch("nPVall",  &nPVall, "nPVall/I");
  tree->Branch("nPV",     &nPV,    "nPV/I");
  tree->Branch("pv_ndof", pv_ndof, "pv_ndof[nPVall]/F");
  tree->Branch("pv_z",    pv_z,    "pv_z[nPVall]/F");
  tree->Branch("pv_rho",  pv_rho,  "pv_rho[nPVall]/F");

  tree->Branch("nEta",   &nEta,  "nEta/I");
  tree->Branch("energy", energy, "energy[nEta]/F");
  tree->Branch("et",     et,     "et[nEta]/F");
  tree->Branch("eRMS",   eRMS,   "eRMS[nEta]/F");

  tree->Branch("etMED",     etMED,     "etMED[nEta]/F");
  tree->Branch("etMEAN",    etMEAN,    "etMEAN[nEta]/F");

  tree->Branch("et_gme",       et_gme,     "et_gme[18][11]/F");
  tree->Branch("et_twopi",     et_twopi,     "et_twopi[18][11]/F");
  tree->Branch("rho_gme",   &rho_gme,  "rho_gme/F");
  if (isMC_) {
//    tree->Branch("et_gme_gen",   et_gme_gen, "et_gme_gen[18][11]/F");
    tree->Branch("et_twopi_gen", et_twopi_gen, "et_twopi_gen[18][11]/F");
    tree->Branch("rho_gme_gen", &rho_gme_gen, "rho_gme_gen/F");
    tree->Branch("et_twopi_det", et_twopi_det, "et_twopi_det[18][11]/F");
    tree->Branch("rho_gme_det", &rho_gme_det, "rho_gme_det/F");
    tree->Branch("et_twopi_scaled", et_twopi_scaled, "et_twopi_scaled[18][11]/F");
    tree->Branch("rho_gme_scaled", &rho_gme_scaled, "rho_gme_scaled/F");
//    tree->Branch("et_gen",     et_gen,     "et_gen[nEta]/F");
//    tree->Branch("etMED_gen",  etMED_gen,  "etMED_gen[nEta]/F");
//    tree->Branch("etMEAN_gen", etMEAN_gen, "etMEAN_gen[nEta]/F");
    tree->Branch("weight",        &weight,        "weight/F");
    tree->Branch("nGenParticles", &nGenParticles, "nGenParticles/I");
    tree->Branch("idHardScat1", &idHardScat1, "idHardScat1/I");
    tree->Branch("idHardScat2", &idHardScat2, "idHardScat2/I");
    tree->Branch("idxHardScat", &idxHardScat, "idxHardScat/I");
    if (writeParticles_) {
      tree->Branch("particle_id",   "std::vector<int>",   &particle_id);
      tree->Branch("particle_pt",   "std::vector<float>", &particle_pt);
      tree->Branch("particle_eta",  "std::vector<float>", &particle_eta);
      tree->Branch("particle_phi",  "std::vector<float>", &particle_phi);
      tree->Branch("particle_et",   "std::vector<float>", &particle_et);
      tree->Branch("particle_charge", "std::vector<int>", &particle_charge);
      tree->Branch("isPromptFinalState", "std::vector<bool>", &isPromptFinalState);
      tree->Branch("particle_et_scaled",   "std::vector<float>", &particle_et_scaled);
    }
  }

  tree->Branch("fchm",   f[chm],   "fchm[nEta]/b");
  tree->Branch("fchu",   f[chu],   "fchu[nEta]/b");
  tree->Branch("fnh",    f[nh],    "fnh[nEta]/b");
  tree->Branch("fne",    f[ne],    "fne[nEta]/b");
  tree->Branch("fhfh",   f[hfh],   "fhfh[nEta]/b");
  tree->Branch("fhfe",   f[hfe],   "fhfe[nEta]/b");
  tree->Branch("flep",   f[lep],   "flep[nEta]/b");
  tree->Branch("funtrk", f[untrk], "funtrk[nEta]/b");


  tree->Branch("ht", &ht, "ht/F");
  tree->Branch("nJets",    &nJets,   "nJets/I");
  tree->Branch("jet_eta",  jet_eta,  "jet_eta[nJets]/F");
  tree->Branch("jet_phi",  jet_phi,  "jet_phi[nJets]/F");
  tree->Branch("jet_pt",   jet_pt,   "jet_pt[nJets]/F");
  tree->Branch("jet_area", jet_area, "jet_area[nJets]/F");

  tree->Branch("jet_ch",  jet_ch,  "jet_ch[nJets]/F");
  tree->Branch("jet_nh",  jet_nh,  "jet_nh[nJets]/F");
  tree->Branch("jet_ne",  jet_ne,  "jet_ne[nJets]/F");
  tree->Branch("jet_hfh", jet_hfh, "jet_hfh[nJets]/F");
  tree->Branch("jet_hfe", jet_hfe, "jet_hfe[nJets]/F");
  tree->Branch("jet_lep", jet_lep, "jet_lep[nJets]/F");
}

// ------------ method called for each event  ------------
void OffsetTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  ++counter;
  if (counter%numSkip_ != 0) return;

//------------ Pileup ------------//

  if (isMC_){
    edm::Handle< vector<PileupSummaryInfo> > pileups;
    iEvent.getByToken(muTag_, pileups);

    mu = pileups->at(1).getTrueNumInteractions();
  }
  else{
    run = int(iEvent.id().run());
    lumi = int(iEvent.getLuminosityBlock().luminosityBlock());
    bx = iEvent.bunchCrossing();
    event = iEvent.id().event();

    mu = getAvgPU( run, lumi );
    int muI = mu;
    if (run>=302031 && run<=302663) muWeight = muWeight_Run2017D[muI];
    else if (run>=303824 && run<=304797) muWeight = muWeight_Run2017E[muI];
    else muWeight = 1;
    if (mu==0) return;
  }

//------------ Primary Vertices ------------//

  edm::Handle< vector<reco::Vertex> > primaryVertices;
  iEvent.getByToken(pvTag_, primaryVertices);

  nPVall = primaryVertices->size();
  nPV = 0;

  for (int i = 0; i != nPVall; ++i){
    reco::Vertex pv = primaryVertices->at(i);

    if( !pv.isFake() && pv.ndof() > 4 && pv.z() <= 24 && pv.position().rho() <= 2 ) ++nPV;

    if(i >= MAXNPV) continue;
    pv_ndof[i] = pv.ndof();
    pv_z[i] = pv.z();
    pv_rho[i] = pv.position().rho();
  }
  if(nPVall > MAXNPV) nPVall = MAXNPV;

//------------ Rho ------------//

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoTag_, rhoHandle);
  rho = *rhoHandle;

  edm::Handle<double> rhoC0Handle;
  iEvent.getByToken(rhoC0Tag_, rhoC0Handle);
  rhoC0 = *rhoC0Handle;

  edm::Handle<double> rhoCCHandle;
  iEvent.getByToken(rhoCCTag_, rhoCCHandle);
  rhoCC = *rhoCCHandle;

  edm::Handle<double> rhoCentralHandle;
  iEvent.getByToken(rhoCentralTag_, rhoCentralHandle);
  rhoCentral = *rhoCentralHandle;

  edm::Handle<double> rhoCentralCaloHandle;
  iEvent.getByToken(rhoCentralCaloTag_, rhoCentralCaloHandle);
  rhoCentralCalo = *rhoCentralCaloHandle;

//------------ Gen Particles -----------//

  if (isMC_) {

    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorTag_, genInfo);
    weight = genInfo->weight();
    //if( genInfo->weight() > 0 ) weight = 1.0; else weight = -1.0;

    edm::Handle<std::vector<reco::GenParticle>> pruned;
    iEvent.getByToken(genparticlesTag_, pruned);
    //const std::vector<reco::GenParticle>>* genps_coll = pruned.product();

    edm::Handle<std::vector<pat::PackedGenParticle>> particles;
    iEvent.getByToken(genTag_, particles);

    nGenParticles = particles->size();
    if( particles->size() == 0 ) return;

//    memset(et_gen, 0, sizeof(et_gen));
//    memset(etMED_gen, 0, sizeof(etMED_gen));
//    memset(etMEAN_gen, 0, sizeof(etMEAN_gen));
    memset(et_twopi_gen,  0, sizeof(et_twopi_gen));
    memset(et_twopi_det,  0, sizeof(et_twopi_det));
    memset(et_twopi_scaled,  0, sizeof(et_twopi_scaled));

    //h2_GME_gen->Reset(); h2_finnereta_gen->Reset();
    h2_twopi_gen->Reset(); h2_twopi_det->Reset(); h2_twopi_scaled->Reset();

    particle_id.clear(); particle_et.clear(); particle_pt.clear(); particle_eta.clear(); particle_phi.clear(); particle_charge.clear();
    isPromptFinalState.clear(); particle_et_scaled.clear();

    int id_hard1 = 0, id_hard2 = 0, n_hard = 0;
    for (size_t i = 0; i < pruned->size(); ++i ) {
      const reco::GenParticle & p = (*pruned)[i];
      if (p.status() == 21) { //incoming particles of the hardest subprocess
        if (n_hard==0) id_hard1 = p.pdgId();
        else if (n_hard==1) id_hard2 = p.pdgId();
        n_hard++;
      }
    }
    int idx_hard = -1;
    if (id_hard1 == 21 && id_hard2 == 21) idx_hard = 2; //gg
    else if (id_hard1 == 21 || id_hard2 == 21) idx_hard = 1; //gq
    else if (id_hard1 != 0 && id_hard2 != 0 && fabs(id_hard1) < 6 && fabs(id_hard2) < 6) idx_hard = 0; //qq

    int nGenPar = 0; //float etEM = 0, etCH = 0, etNH = 0;
    //vector<reco::GenParticle>::const_iterator i_particle, endparticle = particles->end();
    vector<pat::PackedGenParticle>::const_iterator i_particle, endparticle = particles->end();
    for (i_particle = particles->begin(); i_particle != endparticle; ++i_particle) {
//      int etaIndex = getEtaIndex( i_particle->eta() );
      int pdgIdAbs = abs(i_particle->pdgId());
      float newPhi = (i_particle->phi()>=0) ? i_particle->phi() : i_particle->phi() + 2*M_PI;
      if (pdgIdAbs != 12 && pdgIdAbs != 14 && pdgIdAbs != 16) {
//        et_gen[etaIndex] += i_particle->et();
//        h2_GME_gen->Fill(i_particle->eta(),i_particle->phi(),i_particle->et());
//        h2_finnereta_gen->Fill(i_particle->eta(),i_particle->phi(),i_particle->et());
        h2_twopi_gen->Fill(i_particle->eta(),newPhi,i_particle->et());
      }
      int charge = i_particle->charge();
      float et = i_particle->et(); //updated on Nov 18, 2019
      //if (pdgIdAbs == 11 || pdgIdAbs == 22) etEM += et; //electron and photon
      //if (pdgIdAbs == 211 || pdgIdAbs == 321) etCH += et; //pion and kaon
      //if (pdgIdAbs == 111 || pdgIdAbs == 311) etNH += et; //pi0 and K0
      float et_scaled = et; //added on Dec 9th, 2019
      if (charge == 0 && pdgIdAbs != 22 && pdgIdAbs != 12 && pdgIdAbs != 14 && pdgIdAbs != 16) {
        float et_nhHHe = 0.45*et, et_nhHCAL = 0.55*0.5*et, et_nhECAL = 0.55*0.5*et;
        float et_nhECAL_scaled = et_nhECAL * 1.076 * (1-1.403*pow(et_nhECAL,0.646-1));
        et_scaled = et_nhHHe + et_nhHCAL + et_nhECAL_scaled;
      }

      if (pdgIdAbs == 11 || pdgIdAbs == 13 || pdgIdAbs == 15 || ((pdgIdAbs < 11 || pdgIdAbs > 16) && ((pdgIdAbs == 22 && et > 0.3) || (charge!=0 && et > 0.4) || (charge==0 && pdgIdAbs != 22 && et_scaled > 3)))) h2_twopi_scaled->Fill(i_particle->eta(),newPhi,et_scaled);

      if ((pdgIdAbs < 11 || pdgIdAbs > 16) && ((pdgIdAbs == 22 && et < 0.3) || (charge!=0 && et < 0.4) || (charge==0 && pdgIdAbs != 22 && et < 3))) continue;
      ++nGenPar;

      if (writeParticles_) {
        particle_id.push_back( i_particle->pdgId() );
        particle_et.push_back( i_particle->et() );
        particle_pt.push_back( i_particle->pt() );
        particle_eta.push_back( i_particle->eta() );
        particle_phi.push_back( i_particle->phi() );
        particle_charge.push_back( i_particle->charge() );
        isPromptFinalState.push_back( i_particle->isPromptFinalState() );
        particle_et_scaled.push_back( et );
      }

      if (pdgIdAbs != 12 && pdgIdAbs != 14 && pdgIdAbs != 16) h2_twopi_det->Fill(i_particle->eta(),newPhi,et);
    }
    nGenParticles = nGenPar;
    idHardScat1 = id_hard1; idHardScat2 = id_hard2; idxHardScat = idx_hard;

//    for (int ieta = 1; ieta != (ETA_BINS+1); ++ieta){
//      vector<double> x_gen; double et_sum_gen = 0;
//      for (int iphi = 1; iphi != PHI_BINS_GME+1; ++iphi){
//        x_gen.push_back(h2_finnereta_gen->GetBinContent(ieta, iphi));
//        et_sum_gen += h2_finnereta_gen->GetBinContent(ieta, iphi);
//      }
//      sort(x_gen.begin(),x_gen.end());
//      float median_gen = x_gen[5]; // for 11 phi bins
//      etMED_gen[ieta-1] = median_gen;
//      etMEAN_gen[ieta-1] = float(et_sum_gen)/(11);
//    }
  }

//------------ PF Particles ------------//

  edm::Handle< vector<pat::PackedCandidate> > pfCandidates;
  iEvent.getByToken(pfTag_, pfCandidates);

  memset(energy, 0, sizeof(energy));    //reset arrays to zero
  memset(eRMS,   0, sizeof(eRMS));
  memset(et, 0, sizeof(et));
  memset(etMED, 0, sizeof(etMED));
  memset(etMEAN, 0, sizeof(etMEAN));

  nEta = ETA_BINS;

  float eFlavor[numFlavors][ETA_BINS] = {};
  float e2[ETA_BINS] = {};  //energy squared
  int nPart[ETA_BINS] = {}; //number of particles per eta bin
 
  memset(et_gme,        0, sizeof(et_gme));
  memset(et_twopi,      0, sizeof(et_twopi));

  pf_type.clear(); pf_pt.clear(); pf_eta.clear(); pf_phi.clear(); pf_et.clear();
  h2_GME->Reset();
  h2_finnereta->Reset();
  h2_twopi->Reset();

  vector<pat::PackedCandidate>::const_iterator i_pf, endpf = pfCandidates->end();
  for (i_pf = pfCandidates->begin(); i_pf != endpf; ++i_pf) {

    int etaIndex = getEtaIndex( i_pf->eta() );
    Flavor flavor = getFlavor( fabs(i_pf->pdgId()) );

    if (etaIndex == -1 || flavor == X) continue;
/*
    bool attached = false;
    reco::TrackRef pftrack( i_pf->trackRef() );

    if (flavor == chm && !pftrack.isNull() ) { //check charged hadrons ONLY
      
      vector<reco::Vertex>::const_iterator i_pv, endpv = primaryVertices->end();
      for (i_pv = primaryVertices->begin(); i_pv != endpv && !attached; ++i_pv) {
        
        if ( !i_pv->isFake() && i_pv->ndof() >= 4 && fabs(i_pv->z()) < 24 ) {

          reco::Vertex::trackRef_iterator i_vtxTrk, endvtxTrk = i_pv->tracks_end();
          for(i_vtxTrk = i_pv->tracks_begin(); i_vtxTrk != endvtxTrk && !attached; ++i_vtxTrk) {
              
            reco::TrackRef vtxTrk(i_vtxTrk->castTo<reco::TrackRef>());
            if (vtxTrk == pftrack)
              attached = true;
          } 
        }
      }
      if (!attached) flavor = chu; //unmatched charged hadron
    }*/
    float e = i_pf->energy();

    energy[etaIndex] += e;
    et[etaIndex] += i_pf->et();
    eFlavor[flavor][etaIndex] += e;

    h2_GME->Fill(i_pf->eta(),i_pf->phi(),i_pf->et());
    float newPhi = (i_pf->phi()>=0) ? i_pf->phi() : i_pf->phi() + 2*M_PI;
    h2_twopi->Fill(i_pf->eta(),newPhi,i_pf->et());
    h2_finnereta->Fill(i_pf->eta(),i_pf->phi(),i_pf->et());
    e2[etaIndex] += (e*e);
    nPart[etaIndex] ++;

    if (writeCands_) {
      pf_type.push_back( static_cast<int>(flavor) );
      pf_pt.push_back( i_pf->pt() );
      pf_eta.push_back( i_pf->eta() );
      pf_phi.push_back( i_pf->phi() );
      pf_et.push_back( i_pf->et() );
    }
  }

  vector<double> xg(ETA_BINS_GME*PHI_BINS_GME); vector<double> xg_gen(ETA_BINS_GME*PHI_BINS_GME);
  vector<double> xg_det(ETA_BINS_GME*PHI_BINS_GME); vector<double> xg_scaled(ETA_BINS_GME*PHI_BINS_GME);
  for (int ieta = 1; ieta != (ETA_BINS_GME+1); ++ieta){
    for (int iphi = 1; iphi != (PHI_BINS_GME+1); ++iphi){
      int igrid = PHI_BINS_GME*(ieta-1)+(iphi-1);
      et_gme[ieta-1][iphi-1] =  h2_GME->GetBinContent(ieta, iphi);
      et_twopi[ieta-1][iphi-1] =  h2_twopi->GetBinContent(ieta, iphi);
      xg[igrid] = et_twopi[ieta-1][iphi-1]/GRID_AREA;
      if (isMC_) {
//        et_gme_gen[ieta-1][iphi-1] =  h2_GME_gen->GetBinContent(ieta, iphi);
        et_twopi_gen[ieta-1][iphi-1] =  h2_twopi_gen->GetBinContent(ieta, iphi);
        xg_gen[igrid] = et_twopi_gen[ieta-1][iphi-1]/GRID_AREA;
        et_twopi_det[ieta-1][iphi-1] =  h2_twopi_det->GetBinContent(ieta, iphi);
        xg_det[igrid] = et_twopi_det[ieta-1][iphi-1]/GRID_AREA;
        et_twopi_scaled[ieta-1][iphi-1] =  h2_twopi_scaled->GetBinContent(ieta, iphi);
        xg_scaled[igrid] = et_twopi_scaled[ieta-1][iphi-1]/GRID_AREA;
      }
    }
  }
  sort(xg.begin(),xg.end()); rho_gme = 0.5*(xg[98]+xg[99]); // out of 11*18=198 entries
  if (isMC_) {
    sort(xg_gen.begin(),xg_gen.end()); rho_gme_gen = 0.5*(xg_gen[98]+xg_gen[99]);
    sort(xg_det.begin(),xg_det.end()); rho_gme_det = 0.5*(xg_det[98]+xg_det[99]);
    sort(xg_scaled.begin(),xg_scaled.end()); rho_gme_scaled = 0.5*(xg_scaled[98]+xg_scaled[99]);
  }

  for (int ieta = 1; ieta != (ETA_BINS+1); ++ieta){
    vector<double> x;
    double et_sum = 0;
    for (int iphi = 1; iphi != PHI_BINS_GME+1; ++iphi){
      x.push_back(h2_finnereta->GetBinContent(ieta, iphi));
      et_sum += h2_finnereta->GetBinContent(ieta, iphi);
    }
    sort(x.begin(),x.end()); 
    float median = x[5]; // for 11 phi bins
    etMED[ieta-1] = median;
    etMEAN[ieta-1] = float(et_sum)/(11);
  }


//------------ Tracks ------------//
/*
  edm::Handle< vector<reco::Track> > tracks;
  iEvent.getByToken(trackTag_, tracks);

  vector<reco::Track>::const_iterator i_trk, endtrk = tracks->end();
  for (i_trk = tracks->begin(); i_trk != endtrk; ++i_trk) {

    if ( !i_trk->quality(reco::Track::tight) ) continue;
    bool matched = false;

    vector<pat::PackedCandidate>::const_iterator i_pf, endpf = pfCandidates->end();
    for (i_pf = pfCandidates->begin();  i_pf != endpf && !matched; ++i_pf) {

      if ( &(*i_trk) == i_pf->trackRef().get() )
        matched = true;      
    }
    if (matched) continue;

    int etaIndex = getEtaIndex( i_trk->eta() );
    if (etaIndex == -1) continue;

    float e = i_trk->p();

    energy[etaIndex] += e;
    et[etaIndex] += i_trk->pt();
    eFlavor[untrk][etaIndex] += e;
    e2[etaIndex] += (e*e);
    nPart[etaIndex] ++;
  }
*/
  for (int i=0; i != nEta; ++i){

    for (int flav = 0; flav != numFlavors; ++flav){
      UChar_t f_value; float eFlav = eFlavor[flav][i]; float E = energy[i];

      if (eFlav == 0)      f_value = 0;
      else if (eFlav == E) f_value = 255;
      else                 f_value = int( eFlav * 256 / E );

      f[flav][i] = f_value;
    }

    nPart[i] == 0 ? eRMS[i] = 0 : eRMS[i] = sqrt( e2[i]/nPart[i] );
  }


//------------ PF Jets ------------//

  edm::Handle< vector<pat::Jet> > pfJets;
  iEvent.getByToken(pfJetTag_, pfJets);

  ht = 0;
  vector<pat::Jet>::const_iterator i_jet, endjet = pfJets->end();
  for (i_jet = pfJets->begin(); i_jet != endjet; ++i_jet) {

    float pt = i_jet->pt();
    if (pt > 10) ht += pt;
  }

  pfJets->size()<MAXJETS ? nJets = pfJets->size() : nJets = MAXJETS;
  for (int i=0; i != nJets; ++i){
    pat::Jet jet = pfJets->at(i);

    jet_eta[i] = jet.eta();
    jet_phi[i] = jet.phi();
    jet_pt[i] = jet.pt();
    jet_area[i] = jet.jetArea();

    jet_ch[i] = jet.chargedHadronEnergyFraction();
    jet_nh[i] = jet.neutralHadronEnergyFraction();
    jet_ne[i] = jet.photonEnergyFraction();
    jet_hfh[i] = jet.HFHadronEnergyFraction();
    jet_hfe[i] = jet.HFEMEnergyFraction();
    jet_lep[i] = jet.electronEnergyFraction() + jet.muonEnergyFraction();
  }

//------------ Fill Tree ------------//

  tree->Fill();
}


// ------------ method called once each job just after ending the event loop  ------------
void OffsetTreeMaker::endJob() {

    root_file->Write();
    root_file->Close();
}

int OffsetTreeMaker::getEtaIndex(float eta){

  for (int i=0; i != ETA_BINS; ++i){
    if (etabins[i] <= eta && eta < etabins[i+1]) return i;
  }
  if (eta == etabins[ETA_BINS]) return ETA_BINS-1;
  else return -1;
}


OffsetTreeMaker::Flavor OffsetTreeMaker::getFlavor(int id)
{
    if (id == 211) //h
        return chm;     //initially matched charged hadron
    else if (id == 11) //e
        return lep;
    else if (id == 13) //mu
        return lep;
    else if (id == 22) //gamma
        return ne;
    else if (id == 130) //h0
        return nh;
    else if (id == 1) //h_HF
        return hfh;
    else if (id == 2) //egamma_HF
        return hfe;
    else
        return X;
}

//define this as a plug-in
DEFINE_FWK_MODULE(OffsetTreeMaker);
