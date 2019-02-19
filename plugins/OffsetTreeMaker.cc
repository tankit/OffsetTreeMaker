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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Common/interface/Ref.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
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
const int MAXNPV = 50;
const int MAXJETS = 4;


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
    TH2F *h2_GME_gen, *h2_finnereta_gen; // for GenParticle
    TRandom3* rand;

    int nEta;
    float energy[ETA_BINS], eRMS[ETA_BINS], et[ETA_BINS], etMED[ETA_BINS], etMEAN[ETA_BINS];
    float et_gme[ETA_BINS_GME][PHI_BINS_GME];
    float et_gme_gen[ETA_BINS_GME][PHI_BINS_GME];
    UChar_t ch_et_gme[ETA_BINS_GME][PHI_BINS_GME];
    UChar_t ch_et_gme_gen[ETA_BINS_GME][PHI_BINS_GME];
    UChar_t f[numFlavors][ETA_BINS];  //energy fraction by flavor
    float et_gen[ETA_BINS], etMED_gen[ETA_BINS], etMEAN_gen[ETA_BINS];

    ULong64_t event;
    int run, lumi, bx;
    float mu;
    float rho, rhoC0, rhoCC;

    int nPVall, nPV;
    float pv_ndof[MAXNPV], pv_z[MAXNPV], pv_rho[MAXNPV];

    float ht;
    int nJets;
    float jet_eta[MAXJETS], jet_phi[MAXJETS], jet_pt[MAXJETS], jet_area[MAXJETS];
    float jet_ch[MAXJETS], jet_nh[MAXJETS], jet_ne[MAXJETS], jet_hfh[MAXJETS], jet_hfe[MAXJETS], jet_lep[MAXJETS];

    vector<int> pf_type;
    vector<float> pf_pt, pf_eta, pf_phi, pf_et;
    vector<int> particle_type;
    vector<float> particle_et, particle_pt, particle_eta, particle_phi;

    TString RootFileName_;
    string puFileName_;
    int numSkip_;
    bool isMC_, writeCands_;

    edm::EDGetTokenT< vector<reco::Vertex> > pvTag_;
//    edm::EDGetTokenT< vector<reco::Track> > trackTag_;
    edm::EDGetTokenT< vector<PileupSummaryInfo> > muTag_;
    edm::EDGetTokenT< vector<pat::PackedCandidate> > pfTag_;
    edm::EDGetTokenT< vector<pat::PackedGenParticle> > genTag_;
    edm::EDGetTokenT< vector<reco::GenParticle>> genparticlesTag_;
    edm::EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT<double> rhoC0Tag_;
    edm::EDGetTokenT<double> rhoCCTag_;
    edm::EDGetTokenT< vector<pat::Jet> > pfJetTag_;
};

OffsetTreeMaker::OffsetTreeMaker(const edm::ParameterSet& iConfig)
{
  numSkip_ = iConfig.getParameter<int> ("numSkip");
  RootFileName_ = iConfig.getParameter<string>("RootFileName");
  puFileName_ = iConfig.getParameter<string>("puFileName");
  isMC_ = iConfig.getParameter<bool>("isMC");
  writeCands_ = iConfig.getParameter<bool>("writeCands");
  pvTag_ = consumes< vector<reco::Vertex> >( iConfig.getParameter<edm::InputTag>("pvTag") );
//  trackTag_ = consumes< vector<reco::Track> >( iConfig.getParameter<edm::InputTag>("trackTag") );
  muTag_ = consumes< vector<PileupSummaryInfo> >( iConfig.getParameter<edm::InputTag>("muTag") );
  pfTag_ = consumes< vector<pat::PackedCandidate> >( iConfig.getParameter<edm::InputTag>("pfTag") );
  genTag_ = consumes< vector<pat::PackedGenParticle> >( iConfig.getParameter<edm::InputTag>("genTag") );
  genparticlesTag_ = consumes< vector<reco::GenParticle> >( iConfig.getParameter<edm::InputTag>("GenParticles") ),
  rhoTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoTag") );
  rhoC0Tag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoC0Tag") );
  rhoCCTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoCCTag") );
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
  h2_GME_gen = new TH2F("GME_2D_gen", "GME_2D_yPhi_xEta", ETA_BINS_GME, -5.0, 5.0, PHI_BINS_GME , -1*M_PI , M_PI);
  h2_finnereta_gen = new TH2F("finnereta_ET_gen", "yPhi_xEta", ETA_BINS, etabins, PHI_BINS_GME, phibins);

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
    tree->Branch("particle_type", "std::vector<int>",   &particle_type);
    tree->Branch("particle_pt",   "std::vector<float>", &particle_pt);
    tree->Branch("particle_eta",  "std::vector<float>", &particle_eta);
    tree->Branch("particle_phi",  "std::vector<float>", &particle_phi);
    tree->Branch("particle_et",   "std::vector<float>", &particle_et);
  }

  tree->Branch("mu", &mu, "mu/F");

  tree->Branch("rho",   &rho,   "rho/F");
  tree->Branch("rhoC0", &rhoC0, "rhoC0/F");
  tree->Branch("rhoCC", &rhoCC, "rhoCC/F");

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
  tree->Branch("et_gme_gen",   et_gme_gen, "et_gme_gen[18][11]/F");
  tree->Branch("ch_et_gme",       ch_et_gme,     "ch_et_gme[18][11]/b");
  tree->Branch("ch_et_gme_gen",   ch_et_gme_gen, "ch_et_gme_gen[18][11]/b");
  tree->Branch("et_gen",     et_gen,     "et[nEta]/F");
  tree->Branch("etMED_gen",  etMED_gen,  "etMED_gen[nEta]/F");
  tree->Branch("etMEAN_gen", etMEAN_gen, "etMEAN_gen[nEta]/F");

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

  counter++;
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
    if (mu==0) return;
  }

//------------ Primary Vertices ------------//

  edm::Handle< vector<reco::Vertex> > primaryVertices;
  iEvent.getByToken(pvTag_, primaryVertices);

  nPVall = primaryVertices->size();
  nPV = 0;

  for (int i = 0; i != nPVall; ++i){
    reco::Vertex pv = primaryVertices->at(i);

    pv_ndof[i] = pv.ndof();
    pv_z[i] = pv.z();
    pv_rho[i] = pv.position().rho();

    if( !pv.isFake() && pv_ndof[i] > 4 && pv_z[i] <= 24 && pv_rho[i] <= 2 )
      nPV++;
  }

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

//------------ Gen Particles -----------//

  //edm::Handle<std::vector<reco::GenParticle>> particles;
  //iEvent.getByToken(genparticlesTag_, particles);

  edm::Handle<std::vector<pat::PackedGenParticle>> particles;
  iEvent.getByToken(genTag_, particles);

  memset(et_gen, 0, sizeof(et_gen));
  memset(etMED_gen, 0, sizeof(etMED_gen));
  memset(etMEAN_gen, 0, sizeof(etMEAN_gen));

  h2_GME_gen->Reset(); h2_finnereta_gen->Reset();

  particle_type.clear(); particle_et.clear(); particle_pt.clear(); particle_eta.clear(); particle_phi.clear();

  //vector<reco::GenParticle>::const_iterator i_particle, endparticle = particles->end();
  vector<pat::PackedGenParticle>::const_iterator i_particle, endparticle = particles->end();
  for (i_particle = particles->begin(); i_particle != endparticle; ++i_particle) {
    int etaIndex = getEtaIndex( i_particle->eta() );
    if(writeCands_) {
      particle_type.push_back( i_particle->pdgId() );
      particle_et.push_back( i_particle->et() );
      particle_pt.push_back( i_particle->pt() );
      particle_eta.push_back( i_particle->eta() );
      particle_phi.push_back( i_particle->phi() );
    }
    et_gen[etaIndex] += i_particle->et();
    h2_GME_gen->Fill(i_particle->eta(),i_particle->phi(),i_particle->et());
    h2_finnereta_gen->Fill(i_particle->eta(),i_particle->phi(),i_particle->et());
  }

  for (int ieta = 1; ieta != (ETA_BINS+1); ++ieta){
    vector<double> x_gen; double et_sum_gen = 0;
    for (int iphi = 1; iphi != PHI_BINS_GME+1; ++iphi){
      x_gen.push_back(h2_finnereta_gen->GetBinContent(ieta, iphi));
      et_sum_gen += h2_finnereta_gen->GetBinContent(ieta, iphi);
    }
    sort(x_gen.begin(),x_gen.end());
    float median_gen = x_gen[5]; // for 11 phi bins
    etMED_gen[ieta-1] = median_gen;
    etMEAN_gen[ieta-1] = float(et_sum_gen)/(11);
  }

//------------ PF Particles ------------//

  edm::Handle< vector<pat::PackedCandidate> > pfCandidates;
  iEvent.getByToken(pfTag_, pfCandidates);

  memset(energy, 0, sizeof(energy));    //reset arrays to zero
  memset(et, 0, sizeof(et));
  memset(etMED, 0, sizeof(etMED));
  memset(etMEAN, 0, sizeof(etMEAN));

  nEta = ETA_BINS;

  float eFlavor[numFlavors][ETA_BINS] = {};
  float e2[ETA_BINS] = {};  //energy squared
  int nPart[ETA_BINS] = {}; //number of particles per eta bin

  pf_type.clear(); pf_pt.clear(); pf_eta.clear(); pf_phi.clear(); pf_et.clear();
  h2_GME->Reset();
  h2_finnereta->Reset();

  vector<pat::PackedCandidate>::const_iterator i_pf, endpf = pfCandidates->end();
  for (i_pf = pfCandidates->begin(); i_pf != endpf; ++i_pf) {

    int etaIndex = getEtaIndex( i_pf->eta() );
    Flavor flavor = getFlavor( i_pf->pdgId() );

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

  for (int ieta = 1; ieta != (ETA_BINS_GME+1); ++ieta){
    for (int iphi = 1; iphi != PHI_BINS_GME+1; ++iphi){
      et_gme[ieta-1][iphi-1] =  h2_GME->GetBinContent(ieta, iphi);
      ch_et_gme[ieta-1][iphi-1] =  char(min(255, int( h2_GME->GetBinContent(ieta, iphi)/0.1)));
      et_gme_gen[ieta-1][iphi-1] =  h2_GME_gen->GetBinContent(ieta, iphi);
      ch_et_gme_gen[ieta-1][iphi-1] =  char(min(255, int( h2_GME_gen->GetBinContent(ieta, iphi)/0.1)));
    }
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

  if (root_file !=0) {

    root_file->Write();
    delete root_file;
    root_file = 0;
  }
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
