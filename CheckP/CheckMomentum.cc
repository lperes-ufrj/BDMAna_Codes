/*

November 9, 2020
Adapted by Leonardo Peres from a series of Yun-Tse's original codes

*/


#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"


// LArSoft, nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

using CounterMap_t = std::map< int, unsigned int >;
using KinematicMap_t = std::map< int, std::vector< double > >;

void InitTreeGen( TTree* pTree, int& Event, int& Run, int& SubRun, double& MomentumDiff, bool& IsConserved) {

    pTree->Branch("eventNo", &Event);
    pTree->Branch("runNo", &Run);
    pTree->Branch("subRunNo", &SubRun);

    pTree->Branch( "MomentumNotConserved", &MomentumDiff);
    pTree->Branch( "IsMomentumConserved", &IsConserved );

}

int main( int argc, char ** argv ) {

    std::vector< std::string > Filenames;
    Filenames.push_back( argv[1] );

    std::string GenLabel = "bdm";
    art::InputTag MCTruthTag { GenLabel };

    int Event;
    int Run;
    int SubRun;
    double MomentumDiff;
    bool IsPConserved;

    TFile *fOut = new TFile( "BDM_CheckMomentum.root", "RECREATE" );
    TTree *fTree = new TTree( "EventMomentum", "Check Momentum" );

    InitTreeGen( fTree, Event, Run, SubRun, MomentumDiff, IsPConserved);
    // keep track of how many events do not conserve 4-momentum
    int n_ev = 0;
    int num_ev = 0;
    bool stop = true;

    for ( gallery::Event ev( Filenames ); stop/*!ev.atEnd()*/; ev.next() ) {
        
        num_ev++;
        if(num_ev>4999) stop=false;
        Event  = ev.eventAuxiliary().event();
        Run    = ev.eventAuxiliary().run();
        SubRun = ev.eventAuxiliary().subRun();
        
        std::cout << "Processing "
                  << "Run " << ev.eventAuxiliary().run() << ", "
                  << "Event " << ev.eventAuxiliary().event() << std::endl;

        auto const& MCTruthHandle = ev.getValidHandle< std::vector< simb::MCTruth > >( MCTruthTag );
        auto const& MCTruthObjs = *MCTruthHandle;

        for ( size_t iMCTruth = 0; iMCTruth < MCTruthObjs.size(); ++iMCTruth ) {

            simb::MCTruth MCTruthObj = MCTruthObjs[iMCTruth];
            int nParticles = MCTruthObj.NParticles();
            TLorentzVector EventMomentum, iP, fP;

            // Check the momentum and energy conservation in an interaction
            for ( int iParticle = 0; iParticle < nParticles; ++iParticle ) {

                const simb::MCParticle& MCParticleObj = MCTruthObj.GetParticle( iParticle );
                //std::cout << "GENIE Particle[" << iParticle << "]: PdgCode: " << MCParticleObj.PdgCode() << ", StatusCode: " << MCParticleObj.StatusCode() << std::endl;
                const TLorentzVector& iMomentum = MCParticleObj.Momentum();
                // Initial state particles
                if ( MCParticleObj.StatusCode() == 0 ) iP += iMomentum;
                // Stable final state particles and nuclear remnants
                else if (  MCParticleObj.StatusCode() == 1 || MCParticleObj.StatusCode() == 15 ) fP += iMomentum;
            } // Loop over MCParticles for each MCTruth

            EventMomentum = fP - iP;
            if ( std::abs( EventMomentum.M2() ) > 1e-10 ) {
                //std::cout << "MCTruth " << iMCTruth << " doesn't conserve the momentum and energy.  The total 4-momentum is ( " << EventMomentum.Px() << ", " << EventMomentum.Py() << ", " << EventMomentum.Pz() << ", " << EventMomentum.E() << " )." << std::endl;
                IsPConserved = false;
                n_ev+=1;
                MomentumDiff = std::abs( EventMomentum.M2() );
            } else {
                IsPConserved = true;
            }

           
        } // Loop over MCTruth
        
        fTree->Fill();
        MomentumDiff=0;
    }
    
    std::cout << "Number of events that did not conserve 4-momentum: " << n_ev << std::endl;
    fOut->Write(0,TObject::kOverwrite);
    return 0;
}
