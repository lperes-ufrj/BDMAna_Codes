/**
 * @file CheatRecoAtmNu.cc
 * @author Leonardo Peres (leoperes@pos.if.ufrj.br)
 * @brief  Cheater Analyzer to Reconstructed information
 * @version 0.1
 * @date 2021-02-26
 * 
 * @copyright Copyright (c) 2021
 * 
 */

// C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>

// some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "TVector3.h"

// "art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"

// LArSoft, nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
//#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"

using CounterMap_t = std::map<int, unsigned int>;
using KinematicMap_t = std::map<int, std::vector<double>>;

CounterMap_t *InitTreeGen(
    auto &pTree,
    auto &inteventno,
    auto &intrunno,
    auto &intsubrunno,
    auto &nPrimaries,
    auto &nDaughters,
    auto &PFP_PdgCode,
    auto &PFP_Parent,
    auto &PFParticleID,
    auto &nPFParticle,
    auto &nPFParticleVtxDaughters,
    auto &nPFParticleTrkDaughters,
    auto &nPFParticleSwrDaughters,
    auto &VtxStatus,
    auto &VtxID,
    auto &VtxXYZ,
    auto &DaughterTrackLengths,
    auto &DaughterVertexTrack,
    auto &DaughterStartMomentumTrack,
    auto &DaughterEndMomentumTrack,
    auto &intCCNC,
    auto &SwrID,
    auto &SwrDirection,
    auto &SwrDirectionErr,
    auto &SwrShowerStart,
    auto &SwrShowerStartErr,
    auto &PandoraNuIDs)
{

    pTree->Branch("EventNo", &inteventno);
    pTree->Branch("RunNo", &intrunno);
    pTree->Branch("SubRunNo", &intsubrunno);

    pTree->Branch("nPrimaries", &nPrimaries);
    pTree->Branch("nDaughters", &nDaughters);
    pTree->Branch("PFP_PdgCode", &PFP_PdgCode);
    pTree->Branch("PFP_Parent", &PFP_Parent);
    pTree->Branch("PFParticleID", &PFParticleID);
    pTree->Branch("nPFParticle", &nPFParticle);
    pTree->Branch("nPFParticleVtxDaughters", &nPFParticleVtxDaughters);
    pTree->Branch("nPFParticleTrkDaughters", &nPFParticleTrkDaughters);
    pTree->Branch("nPFParticleSwrDaughters", &nPFParticleSwrDaughters);

    pTree->Branch("VtxStatus", &VtxStatus);
    pTree->Branch("VtxID", &VtxID);
    pTree->Branch("VtxXYZ", &VtxXYZ);

    pTree->Branch("DaughterTrackLengths", &DaughterTrackLengths);
    pTree->Branch("DaughterVertexTrack", &DaughterVertexTrack);
    pTree->Branch("DaughterStartMomentumTrack", &DaughterStartMomentumTrack);
    pTree->Branch("DaughterEndMomentumTrack", &DaughterEndMomentumTrack);
    pTree->Branch("intCCNC", &intCCNC);

    pTree->Branch("SwrID", &SwrID);
    pTree->Branch("SwrDirection", &SwrDirection);
    pTree->Branch("SwrDirectionErr", &SwrDirectionErr);
    pTree->Branch("SwrShowerStart", &SwrShowerStart);
    pTree->Branch("SwrShowerStartErr", &SwrShowerStartErr);

    pTree->Branch("PandoraNuIDs", &PandoraNuIDs);
}

void ResetCounters(auto &EventNo,
                   auto &RunNo,
                   auto &SubRunNo,
                   auto &nPrimaries,
                   auto &nDaughters,
                   auto &PFP_PdgCode,
                   auto &PFP_Parent,
                   auto &PFParticleID,
                   auto &nPFParticle,
                   auto &nPFParticleVtxDaughters,
                   auto &nPFParticleTrkDaughters,
                   auto &nPFParticleSwrDaughters,
                   auto &VtxStatus,
                   auto &VtxID,
                   auto &VtxXYZ,
                   auto &DaughterTrackLengths,
                   auto &DaughterVertexTrack,
                   auto &DaughterStartMomentumTrack,
                   auto &DaughterEndMomentumTrack,
                   auto &intCCNC,
                   auto &SwrID,
                   auto &SwrDirection,
                   auto &SwrDirectionErr,
                   auto &SwrShowerStart,
                   auto &SwrShowerStartErr,
                   auto &PandoraNuIDs)
{

    EventNo = 0;
    RunNo = 0;
    SubRunNo = 0;
    nPrimaries = 0;
    nDaughters.resize(0);
    PFP_PdgCode.resize(0);
    PFP_Parent.resize(0);
    PFParticleID.resize(0);
    nPFParticle = 0;
    nPFParticleVtxDaughters = 0;
    nPFParticleTrkDaughters = 0;
    nPFParticleSwrDaughters = 0;
    VtxStatus.resize(0);
    VtxID.resize(0);
    VtxXYZ.resize(0);
    DaughterTrackLengths.resize(0);
    DaughterVertexTrack.resize(0);
    DaughterStartMomentumTrack.resize(0);
    DaughterEndMomentumTrack.resize(0);
    intCCNC.resize(0);
    SwrID.resize(0);
    SwrDirection.resize(0);
    SwrDirectionErr.resize(0);
    SwrShowerStart.resize(0);
    SwrShowerStartErr.resize(0);
    PandoraNuIDs.resize(0);
}

int main(int argc, char **argv)
{

    std::vector<std::string> filenames;
    filenames.push_back(argv[1]);

    std::string PandoraLabel = "pandora";
    art::InputTag fPandoraLabel{PandoraLabel};

    std::string GeantLabel = "largeant";
    art::InputTag fGeantLabel{GeantLabel};

    std::string GeneratorLabel = "generator";
    art::InputTag fGeneratorLabel{GeneratorLabel};

    art::InputTag const fPandoraTrackModuleLabel = {
        fPandoraLabel.label() + "Track",
        fPandoraLabel.instance(),
        fPandoraLabel.process()};

    art::InputTag const fPandoraShowerModuleLabel = {
        fPandoraLabel.label() + "Shower",
        fPandoraLabel.instance(),
        fPandoraLabel.process()};

    TFile *fOut = new TFile("Bg_RecoParticles.root", "RECREATE");
    TTree *fTree = new TTree("RecoParticles", " Reco Info from atmnu sample");

    int inteventno, intrunno, intsubrunno, nPrimaries, nPFParticle;
    std::vector<double> nPFParticleVtxDaughters, nDaughters, nPFParticleTrkDaughters, nPFParticleSwrDaughters, 
    DaughterTrackLengths, DaughterStartMomentumTrack, DaughterEndMomentumTrack, PandoraNuIDs,SwrID;
    std::vector<int> PFParticleID, VtxStatus, VtxID, PFP_PdgCode, PFP_Parent, intCCNC;
    double XYZ[3] = {0.0, 0.0, 0.0};
    std::vector<TVector3> VtxXYZ, DaughterVertexTrack, SwrDirection, SwrDirectionErr,SwrShowerStart,SwrShowerStartErr;

    InitTreeGen(fTree,
                inteventno,
                intrunno,
                intsubrunno,
                nPrimaries,
                nDaughters,
                PFP_PdgCode,
                PFP_Parent,
                PFParticleID,
                nPFParticle,
                nPFParticleVtxDaughters,
                nPFParticleTrkDaughters,
                nPFParticleSwrDaughters,
                VtxStatus,
                VtxID,
                VtxXYZ,
                DaughterTrackLengths,
                DaughterVertexTrack,
                DaughterStartMomentumTrack,
                DaughterEndMomentumTrack,
                intCCNC,
                SwrID,
                SwrDirection,
                SwrDirectionErr,
                SwrShowerStart,
                SwrShowerStartErr,
                PandoraNuIDs);

    for (gallery::Event ev(filenames); !ev.atEnd(); ev.next())
    {

        ResetCounters(inteventno,
                      intrunno,
                      intsubrunno,
                      nPrimaries,
                      nDaughters,
                      PFP_PdgCode,
                      PFP_Parent,
                      PFParticleID,
                      nPFParticle,
                      nPFParticleVtxDaughters,
                      nPFParticleTrkDaughters,
                      nPFParticleSwrDaughters,
                      VtxStatus,
                      VtxID,
                      VtxXYZ,
                      DaughterTrackLengths,
                      DaughterVertexTrack,
                      DaughterStartMomentumTrack,
                      DaughterEndMomentumTrack,
                      intCCNC,
                      SwrID,
                      SwrDirection,
                      SwrDirectionErr,
                      SwrShowerStart,
                      SwrShowerStartErr,
                      PandoraNuIDs);

        auto const &evaux = ev.eventAuxiliary();
        intrunno = evaux.run();
        intsubrunno = evaux.subRun();
        inteventno = evaux.event();

        // EVENT INFO ====================================================

        std::cout << "Processing "
                  << "Run " << ev.eventAuxiliary().run() << ", "
                  << "Event " << ev.eventAuxiliary().event() << std::endl;

        // TRUE INFORMATION ===============================================

        auto const &MCParticleHandle = ev.getValidHandle<std::vector<simb::MCParticle>>(fGeantLabel);
        auto const &MCParticleObjs = *MCParticleHandle;

        auto const &MCTruthHandle = ev.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorLabel);
        auto const &MCTruthObjs = *MCTruthHandle;

        // SAVE TRUE INTERACTION CC OR NC =================================

        for (size_t i = 0; i < MCTruthObjs.size(); i++)
        {
            simb::MCTruth MCTruthObj = MCTruthObjs[i];
            simb::MCNeutrino const &thisNu = MCTruthObj.GetNeutrino();
            intCCNC.push_back(thisNu.CCNC());
            std::cout << "thisNu.CCNC() = " << thisNu.CCNC() << std::endl;
        }

        // RECO INFORMATION ===============================================

        auto const &TrackListHandle = ev.getValidHandle<std::vector<recob::Track>>(fPandoraTrackModuleLabel);
        auto const &PandoraTrackObjs = *TrackListHandle;

        auto const &PfparticleHandle = ev.getValidHandle<std::vector<recob::PFParticle>>(fPandoraLabel);
        auto const &PandoraPfparticleObjs = *PfparticleHandle;

        auto const &VertexHandle = ev.getValidHandle<std::vector<recob::Vertex>>(fPandoraLabel);
        auto const &PandoraVertexObjs = *VertexHandle;

        // ASSOCIATIONS PFP+TRACK, PFP+SHOWER, PFP+VERTEX ==================
        art::FindMany<recob::Track> PandoraPfp_trk(PfparticleHandle, ev, fPandoraTrackModuleLabel);
        art::FindMany<recob::Shower> PandoraPfp_Swr(PfparticleHandle, ev, fPandoraShowerModuleLabel);
        art::FindMany<recob::Vertex> PandoraPfp_vtx(PfparticleHandle, ev, PandoraLabel);

        // SAVE ONLY PFPARTICLE INFORMATION ================================

        nPFParticle = PandoraPfparticleObjs.size();
        //int neutrinoID = 99999;

        for (size_t iPfp = 0; iPfp < PandoraPfparticleObjs.size(); iPfp++)
        {

            PFP_PdgCode.push_back(PandoraPfparticleObjs[iPfp].PdgCode());
            if (PandoraPfparticleObjs[iPfp].IsPrimary())
                nPrimaries++;
            nDaughters.push_back(PandoraPfparticleObjs[iPfp].NumDaughters());
            PFParticleID.push_back(PandoraPfparticleObjs[iPfp].Self());
            PFP_Parent.push_back(PandoraPfparticleObjs[iPfp].Parent());
            if (!(PandoraPfparticleObjs[iPfp].IsPrimary() && (std::abs(PandoraPfparticleObjs[iPfp].PdgCode()) == 14 || std::abs(PandoraPfparticleObjs[iPfp].PdgCode()) == 12)))
                PandoraNuIDs.push_back(PandoraPfparticleObjs[iPfp].Self());
        }


        for (size_t iPfp = 0; iPfp < PandoraPfparticleObjs.size(); iPfp++)
        {

            auto const &Pfp_trk = PandoraPfp_trk.at(iPfp);

            nPFParticleTrkDaughters += Pfp_trk.size();
            if (!Pfp_trk.empty())
            {
                for (size_t i = 0; i < Pfp_trk.size(); i++)
                {
                    auto const thisTrkPfp = Pfp_trk[i];
                    DaughterTrackLengths.push_back(thisTrkPfp->Length());
                    DaughterStartMomentumTrack.push_back(thisTrkPfp->StartMomentum());
                    DaughterEndMomentumTrack.push_back(thisTrkPfp->EndMomentum());
                }
            }
        }

        for (size_t iPfp = 0; iPfp < PandoraPfparticleObjs.size(); iPfp++)
        {

            auto const &Pfp_vtx = PandoraPfp_vtx.at(iPfp);

            nPFParticleVtxDaughters += Pfp_vtx.size();

            for (size_t i = 0; i < Pfp_vtx.size(); i++)
            {
                auto const thisVtxPfp = Pfp_vtx[i];
                VtxStatus.push_back(thisVtxPfp->status());
                VtxID.push_back(thisVtxPfp->ID());
                thisVtxPfp->XYZ(XYZ);
                VtxXYZ.push_back(XYZ);
            }
        }

        for (size_t iPfp = 0; iPfp < PandoraPfparticleObjs.size(); iPfp++)
        {

            auto const &Pfp_swr = PandoraPfp_Swr.at(iPfp);

            nPFParticleSwrDaughters += Pfp_swr.size();

            for (size_t i = 0; i < Pfp_swr.size(); i++)
            {

                auto const thisSwrPfp = Pfp_swr[i];

                SwrID.push_back(thisSwrPfp->ID());
                SwrDirection.push_back(thisSwrPfp->Direction());
                SwrDirectionErr.push_back(thisSwrPfp->DirectionErr());
                SwrShowerStart.push_back(thisSwrPfp->ShowerStart());
                SwrShowerStartErr.push_back(thisSwrPfp->ShowerStartErr());
            }
        }

        fTree->Fill();
    }

    fOut->Write();
    return 0;
}
