/*

November 5, 2020
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

//particles pdg code
#define ALL_PART_CODE 0
#define IN_DM_CODE -2000010000
#define OUT_DM_CODE 2000010000
#define PROTONS_CODE 2212
#define NEUTRONS_CODE 2112
#define CH_PIONS_CODE 211
#define NE_PIONS_CODE 111
#define MESONS_CODE 300
#define BARYONS_CODE 3000
#define ARGON40_CODE 1000180400
#define ARGON39_CODE 1000180390
#define CL39_CODE 1000170390
#define GEN_PART_CODE 2000000000
#define VISIBLE_PART_CODE 2000000400
  
using CounterMap_t = std::map< int, unsigned int >;
using KinematicMap_t = std::map< int, std::vector< double > >;

void SetParticleTypes( std::vector< int >& ParticleTypes ) {

    ParticleTypes.push_back( ALL_PART_CODE );  // All particles, including visible and invisible
    ParticleTypes.push_back( PROTONS_CODE );  // protons
    ParticleTypes.push_back( NEUTRONS_CODE );  // neutrons
    ParticleTypes.push_back( CH_PIONS_CODE  );  // charged pions
    ParticleTypes.push_back( NE_PIONS_CODE  );  // neutral pions
    ParticleTypes.push_back( MESONS_CODE  );  // Mesons including K0(311), K+(321), KL(130), but not charged and neutral pions
    ParticleTypes.push_back( BARYONS_CODE );  // Baryons including Lambda(3122), Sigma-(3112), Sigma+(3222), but not protons and neutrons
    ParticleTypes.push_back( ARGON40_CODE );  // argon 40
    ParticleTypes.push_back( ARGON39_CODE );  // argon 39
    ParticleTypes.push_back( CL39_CODE );  // Cl 39
    ParticleTypes.push_back( IN_DM_CODE );  // Incident dark matter
    ParticleTypes.push_back( OUT_DM_CODE ); // Outgoing dark matter
    ParticleTypes.push_back( GEN_PART_CODE );  // GENIE artifacts
    ParticleTypes.push_back( VISIBLE_PART_CODE ); // All the visible particles  ( not GENIE artifacts, nor DMs, Ar40, Ar39, Cl39 and neutrons )
 /* ParticleTypes.push_back( 2000000401 ); // All the visible particles but not neutrons
    ParticleTypes.push_back( 2000000410 ); // Leading particle among all the visible particles and neutrons
    ParticleTypes.push_back( 2000000411 ); // Leading particle among all the visible particles but not neutrons
    ParticleTypes.push_back( 2000000412 ); // Leading proton */

}

// CalculateAngle accepts two vectors (in our case, 3-vectors) and calculates the angle between them
double CalculateAngle( TLorentzVector DMP, TLorentzVector ParticleP ) {

    double angle;  
    if(DMP.Vect().Unit()==ParticleP.Vect().Unit()){
       
       angle = 0.0;     

    } else {

        angle = acos(DMP.Vect().Unit()*ParticleP.Vect().Unit());
       // std::cout << angle <<std::endl;

    }
    return angle;

}

void ResetCounters( CounterMap_t& Multiplicity, KinematicMap_t& Px, KinematicMap_t& Py, KinematicMap_t& Pz, KinematicMap_t& P, KinematicMap_t& E, KinematicMap_t& Angle, KinematicMap_t& InitPosX, KinematicMap_t& InitPosY, KinematicMap_t& InitPosZ, /*KinematicMap_t& EndPosX, KinematicMap_t& EndPosY, KinematicMap_t& EndPosZ,*/ KinematicMap_t& Mass ) {

    std::vector< int > ParticleTypes;
    SetParticleTypes( ParticleTypes );

    for ( size_t iType = 0; iType < ParticleTypes.size(); ++iType ) {
        int type = ParticleTypes[iType];
        Multiplicity[type] = 0;
        Px[type].resize( 0 );
        Py[type].resize( 0 );
        Pz[type].resize( 0 );
        P[type].resize( 0 );
        E[type].resize( 0 );
        Angle[type].resize( 0 );
        InitPosX[type].resize( 0 );
        InitPosY[type].resize( 0 );
        InitPosZ[type].resize( 0 );
        Mass[type].resize( 0 );
     // EndPosX[type].clear();
     // EndPosY[type].clear();
     // EndPosZ[type].clear();
    }
}


void SaveInfo( auto const& TruthObj, int PDGCode, CounterMap_t& MultiplicitySI, KinematicMap_t& PxSI, KinematicMap_t& PySI, KinematicMap_t& PzSI, KinematicMap_t& PSI, KinematicMap_t& ESI, TLorentzVector& DMMomemtumSI, KinematicMap_t& AngleSI, KinematicMap_t& InitPosXSI, KinematicMap_t& InitPosYSI, KinematicMap_t& InitPosZSI, /*KinematicMap_t& EndPosXSI, KinematicMap_t& EndPosYSI, KinematicMap_t& EndPosZSI,*/ KinematicMap_t& MassSI){

    double ang = CalculateAngle( DMMomemtumSI, TruthObj.Momentum(0) );
    MultiplicitySI[PDGCode]++;
                    
    InitPosXSI[PDGCode].push_back(TruthObj.Vx(0));
    InitPosYSI[PDGCode].push_back(TruthObj.Vy(0));
    InitPosZSI[PDGCode].push_back(TruthObj.Vz(0));
    //EndPosXSI[PDGCode].push_back(TruthObj.EndX());
    //EndPosYSI[PDGCode].push_back(TruthObj.EndY());
    //EndPosZSI[PDGCode].push_back(TruthObj.EndZ());

    PxSI[PDGCode].push_back(TruthObj.Px(0));
    PySI[PDGCode].push_back(TruthObj.Py(0));
    PzSI[PDGCode].push_back(TruthObj.Pz(0));

    PSI[PDGCode].push_back(TruthObj.P(0));
    ESI[PDGCode].push_back(TruthObj.E(0));
    MassSI[PDGCode].push_back(TruthObj.Mass());
    AngleSI[PDGCode].push_back( ang );

}


CounterMap_t* InitTreeGen( TTree* pTree, CounterMap_t& MultipGen, KinematicMap_t& Px, KinematicMap_t& Py, KinematicMap_t& Pz, KinematicMap_t& P, KinematicMap_t& E, KinematicMap_t& Angle, KinematicMap_t& InitPosX, KinematicMap_t& InitPosY, KinematicMap_t& InitPosZ, /*KinematicMap_t& EndPosX, KinematicMap_t& EndPosY, KinematicMap_t& EndPosZ, */KinematicMap_t& Mass ) {

  //  CounterMap_t* pMultipGen = new CounterMap_t;
  //  CounterMap_t& MultipGen = (*pMultipGen);
    
    //ResetCounters( MultipGen,  Px, Py, Pz, P, E, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);

    pTree->Branch( "nParticlesGen", &MultipGen[ALL_PART_CODE], "nParticleGen/i" );
    pTree->Branch( "nProtonsGen", &MultipGen[PROTONS_CODE], "nProtonsGen/i" );
    pTree->Branch( "nNeutronsGen", &MultipGen[NEUTRONS_CODE], "nNeutronsGen/i" );
    pTree->Branch( "nPionsGen", &MultipGen[CH_PIONS_CODE], "nPionsGen/i" );
    pTree->Branch( "nPi0sGen", &MultipGen[NE_PIONS_CODE], "nPi0sGen/i" );
    pTree->Branch( "nMesonsGen", &MultipGen[MESONS_CODE], "nMesonsGen/i" );
    pTree->Branch( "nBaryonsGen", &MultipGen[BARYONS_CODE], "nBaryonsGen/i" );
    pTree->Branch( "nAr40Gen", &MultipGen[ARGON40_CODE], "nAr40Gen/i" );
    pTree->Branch( "nAr39Gen", &MultipGen[ARGON39_CODE], "nAr39Gen/i" );
    pTree->Branch( "nCl39Gen", &MultipGen[CL39_CODE], "nCl39Gen/i" );
    pTree->Branch( "nDMsGen", &MultipGen[OUT_DM_CODE], "nDMsGen/i" );
    pTree->Branch( "nGENIEGen", &MultipGen[GEN_PART_CODE], "nGENIEGen/i" );
    pTree->Branch( "nVisibleGen", &MultipGen[VISIBLE_PART_CODE], "nVisibleGen/i" );

    
    // Variables for each final state particle
    pTree->Branch( "AllPartPx", &Px[ALL_PART_CODE] );
    pTree->Branch( "AllPartPy", &Py[ALL_PART_CODE] );
    pTree->Branch( "AllPartPz", &Pz[ALL_PART_CODE] );
    pTree->Branch( "AllPartP", &P[ALL_PART_CODE] );
    pTree->Branch( "AllPartE", &E[ALL_PART_CODE] );
    pTree->Branch( "AllPartAngle", &Angle[ALL_PART_CODE] );
    pTree->Branch( "AllPartInitPosX", &InitPosX[ALL_PART_CODE] );
    pTree->Branch( "AllPartInitPosY", &InitPosY[ALL_PART_CODE] );
    pTree->Branch( "AllPartInitPosZ", &InitPosZ[ALL_PART_CODE] );
    pTree->Branch( "AllPartMass", &Mass[ALL_PART_CODE] );


    pTree->Branch( "VisiblePx", &Px[VISIBLE_PART_CODE] );
    pTree->Branch( "VisiblePy", &Py[VISIBLE_PART_CODE] );
    pTree->Branch( "VisiblePz", &Pz[VISIBLE_PART_CODE] );
    pTree->Branch( "VisibleP", &P[VISIBLE_PART_CODE] );
    pTree->Branch( "VisibleE", &E[VISIBLE_PART_CODE] );
    pTree->Branch( "VisibleAngle", &Angle[VISIBLE_PART_CODE] );
    pTree->Branch( "VisibleInitPosX", &InitPosX[VISIBLE_PART_CODE] );
    pTree->Branch( "VisibleInitPosY", &InitPosY[VISIBLE_PART_CODE] );
    pTree->Branch( "VisibleInitPosZ", &InitPosZ[VISIBLE_PART_CODE] );
    pTree->Branch( "VisibleMass", &Mass[VISIBLE_PART_CODE] );
    

    pTree->Branch( "ProtonPx", &Px[PROTONS_CODE] );
    pTree->Branch( "ProtonPy", &Py[PROTONS_CODE] );
    pTree->Branch( "ProtonPz", &Pz[PROTONS_CODE] );
    pTree->Branch( "ProtonP", &P[PROTONS_CODE] );
    pTree->Branch( "ProtonE", &E[PROTONS_CODE] );
    pTree->Branch( "ProtonAngle", &Angle[PROTONS_CODE] );
    pTree->Branch( "ProtonInitPosX", &InitPosX[PROTONS_CODE] );
    pTree->Branch( "ProtonInitPosY", &InitPosY[PROTONS_CODE] );
    pTree->Branch( "ProtonInitPosZ", &InitPosZ[PROTONS_CODE] );
    pTree->Branch( "ProtonMass", &Mass[PROTONS_CODE] );

    
    pTree->Branch( "NeutronPx", &Px[NEUTRONS_CODE] );
    pTree->Branch( "NeutronPy", &Py[NEUTRONS_CODE] );
    pTree->Branch( "NeutronPz", &Pz[NEUTRONS_CODE] );
    pTree->Branch( "NeutronP", &P[NEUTRONS_CODE] );
    pTree->Branch( "NeutronE", &E[NEUTRONS_CODE] );
    pTree->Branch( "NeutronAngle", &Angle[NEUTRONS_CODE] );
    pTree->Branch( "NeutronInitPosX", &InitPosX[NEUTRONS_CODE] );
    pTree->Branch( "NeutronInitPosY", &InitPosY[NEUTRONS_CODE] );
    pTree->Branch( "NeutronInitPosZ", &InitPosZ[NEUTRONS_CODE] );
    pTree->Branch( "NeutronMass", &Mass[NEUTRONS_CODE] );


    pTree->Branch( "PionPx", &Px[CH_PIONS_CODE] );
    pTree->Branch( "PionPy", &Py[CH_PIONS_CODE] );
    pTree->Branch( "PionPz", &Pz[CH_PIONS_CODE] );
    pTree->Branch( "PionP", &P[CH_PIONS_CODE] );
    pTree->Branch( "PionE", &E[CH_PIONS_CODE] );
    pTree->Branch( "PionAngle", &Angle[CH_PIONS_CODE] );
    pTree->Branch( "PionInitPosX", &InitPosX[CH_PIONS_CODE] );
    pTree->Branch( "PionInitPosY", &InitPosY[CH_PIONS_CODE] );
    pTree->Branch( "PionInitPosZ", &InitPosZ[CH_PIONS_CODE] );
    pTree->Branch( "PionMass", &Mass[CH_PIONS_CODE] );


    pTree->Branch( "Pi0Px", &Px[NE_PIONS_CODE] );
    pTree->Branch( "Pi0Py", &Py[NE_PIONS_CODE] );
    pTree->Branch( "Pi0Pz", &Pz[NE_PIONS_CODE] );
    pTree->Branch( "Pi0P", &P[NE_PIONS_CODE] );
    pTree->Branch( "Pi0E", &E[NE_PIONS_CODE] );
    pTree->Branch( "Pi0Angle", &Angle[NE_PIONS_CODE] );
    pTree->Branch( "Pi0InitPosX", &InitPosX[NE_PIONS_CODE] );
    pTree->Branch( "Pi0InitPosY", &InitPosY[NE_PIONS_CODE] );
    pTree->Branch( "Pi0InitPosZ", &InitPosZ[NE_PIONS_CODE] );
    pTree->Branch( "Pi0Mass", &Mass[NE_PIONS_CODE] );

    
    pTree->Branch( "MesonPx", &Px[MESONS_CODE] );
    pTree->Branch( "MesonPy", &Py[MESONS_CODE] );
    pTree->Branch( "MesonPz", &Pz[MESONS_CODE] );
    pTree->Branch( "MesonE", &E[MESONS_CODE] );
    pTree->Branch( "MesonAngle", &Angle[MESONS_CODE] );
    pTree->Branch( "MesonInitPosX", &InitPosX[MESONS_CODE] );
    pTree->Branch( "MesonInitPosY", &InitPosY[MESONS_CODE] );
    pTree->Branch( "MesonInitPosZ", &InitPosZ[MESONS_CODE] );
    pTree->Branch( "MesonMass", &Mass[MESONS_CODE] );
    

    pTree->Branch( "BaryonPx", &Px[BARYONS_CODE] );
    pTree->Branch( "BaryonPy", &Py[BARYONS_CODE] );
    pTree->Branch( "BaryonPz", &Pz[BARYONS_CODE] );
    pTree->Branch( "BaryonP", &P[BARYONS_CODE] );
    pTree->Branch( "BaryonE", &E[BARYONS_CODE] );
    pTree->Branch( "BaryonAngle", &Angle[BARYONS_CODE] );
    pTree->Branch( "BaryonInitPosX", &InitPosX[BARYONS_CODE] );
    pTree->Branch( "BaryonInitPosY", &InitPosY[BARYONS_CODE] );
    pTree->Branch( "BaryonInitPosZ", &InitPosZ[BARYONS_CODE] );
    pTree->Branch( "BaryonMass", &Mass[BARYONS_CODE] );

    
    pTree->Branch( "Ar40Px", &Px[ARGON40_CODE] );
    pTree->Branch( "Ar40Py", &Py[ARGON40_CODE] );
    pTree->Branch( "Ar40Pz", &Pz[ARGON40_CODE] );
    pTree->Branch( "Ar40P", &P[ARGON40_CODE] );
    pTree->Branch( "Ar40E", &E[ARGON40_CODE] );
    pTree->Branch( "Ar40Angle", &Angle[ARGON40_CODE] );
    pTree->Branch( "Ar40InitPosX", &InitPosX[ARGON40_CODE] );
    pTree->Branch( "Ar40InitPosY", &InitPosY[ARGON40_CODE] );
    pTree->Branch( "Ar40InitPosZ", &InitPosZ[ARGON40_CODE] );
    pTree->Branch( "Ar40Mass", &Mass[ARGON40_CODE] );


    pTree->Branch( "Ar39Px", &Px[ARGON39_CODE] );
    pTree->Branch( "Ar39Py", &Py[ARGON39_CODE] );
    pTree->Branch( "Ar39Pz", &Pz[ARGON39_CODE] );
    pTree->Branch( "Ar39P", &P[ARGON39_CODE] );
    pTree->Branch( "Ar39E", &E[ARGON39_CODE] );
    pTree->Branch( "Ar39Angle", &Angle[ARGON39_CODE] );
    pTree->Branch( "Ar39InitPosX", &InitPosX[ARGON39_CODE] );
    pTree->Branch( "Ar39InitPosY", &InitPosY[ARGON39_CODE] );
    pTree->Branch( "Ar39InitPosZ", &InitPosZ[ARGON39_CODE] );
    pTree->Branch( "Ar39Mass", &Mass[ARGON39_CODE] );


    pTree->Branch( "Cl39Px", &Px[CL39_CODE] );
    pTree->Branch( "Cl39Py", &Py[CL39_CODE] );
    pTree->Branch( "Cl39Pz", &Pz[CL39_CODE] );
    pTree->Branch( "Cl39P", &P[CL39_CODE] );
    pTree->Branch( "Cl39E", &E[CL39_CODE] );
    pTree->Branch( "Cl39Angle", &Angle[CL39_CODE] );
    pTree->Branch( "Cl39InitPosX", &InitPosX[CL39_CODE] );
    pTree->Branch( "Cl39InitPosY", &InitPosY[CL39_CODE] );
    pTree->Branch( "Cl39InitPosZ", &InitPosZ[CL39_CODE] );
    pTree->Branch( "Cl39Mass", &Mass[CL39_CODE] );


    pTree->Branch( "InDMPx", &Px[IN_DM_CODE] );
    pTree->Branch( "InDMPy", &Py[IN_DM_CODE] );
    pTree->Branch( "InDMPz", &Pz[IN_DM_CODE] );
    pTree->Branch( "InDMP", &P[IN_DM_CODE] );
    pTree->Branch( "InDME", &E[IN_DM_CODE] );
    pTree->Branch( "InDMAngle", &Angle[IN_DM_CODE] );
    pTree->Branch( "InDMInitPosX", &InitPosX[IN_DM_CODE] );
    pTree->Branch( "InDMInitPosY", &InitPosY[IN_DM_CODE] );
    pTree->Branch( "InDMInitPosZ", &InitPosZ[IN_DM_CODE] );
    pTree->Branch( "InDMMass", &Mass[IN_DM_CODE] );


    pTree->Branch( "OutDMPx", &Px[OUT_DM_CODE] );
    pTree->Branch( "OutDMPy", &Py[OUT_DM_CODE] );
    pTree->Branch( "OutDMPz", &Pz[OUT_DM_CODE] );
    pTree->Branch( "OutDMP", &P[OUT_DM_CODE] );
    pTree->Branch( "OutDME", &E[OUT_DM_CODE] );
    pTree->Branch( "OutDMAngle", &Angle[OUT_DM_CODE] );
    pTree->Branch( "OutDMInitPosX", &InitPosX[OUT_DM_CODE] );
    pTree->Branch( "OutDMInitPosY", &InitPosY[OUT_DM_CODE] );
    pTree->Branch( "OutDMInitPosZ", &InitPosZ[OUT_DM_CODE] );
    pTree->Branch( "OutDMMass", &Mass[OUT_DM_CODE] );
   

    pTree->Branch( "GENIEPx", &Px[GEN_PART_CODE] );
    pTree->Branch( "GENIEPy", &Py[GEN_PART_CODE] );
    pTree->Branch( "GENIEPz", &Pz[GEN_PART_CODE] );
    pTree->Branch( "GENIEP", &P[GEN_PART_CODE] );
    pTree->Branch( "GENIEE", &E[GEN_PART_CODE] );
    pTree->Branch( "GENIEAngle", &Angle[GEN_PART_CODE] );
    pTree->Branch( "GENIEInitPosX", &InitPosX[GEN_PART_CODE] );
    pTree->Branch( "GENIEInitPosY", &InitPosY[GEN_PART_CODE] );
    pTree->Branch( "GENIEInitPosZ", &InitPosZ[GEN_PART_CODE] );
    pTree->Branch( "GENIEMass", &Mass[GEN_PART_CODE] );


 /*
    pTree->Branch( "ProtonEndPosX", &EndPosX[PROTONS_CODE] );
    pTree->Branch( "NeutronEndPosX", &EndPosX[NEUTRONS_CODE] );
    pTree->Branch( "PionEndPosX", &EndPosX[CH_PIONS_CODE] );
    pTree->Branch( "Pi0EndPosX", &EndPosX[NE_PIONS_CODE] );
    pTree->Branch( "MesonEndPosX", &EndPosX[MESONS_CODE] );
    pTree->Branch( "BaryonEndPosX", &EndPosX[BARYONS_CODE] );
    pTree->Branch( "InDMEndPosX", &EndPosX[IN_DM_CODE] );
    pTree->Branch( "OutDMEndPosX", &EndPosX[OUT_DM_CODE] );
    pTree->Branch( "GENIEEndPosX", &EndPosX[GEN_PART_CODE] );

    pTree->Branch( "ProtonEndPosY", &EndPosY[PROTONS_CODE] );
    pTree->Branch( "NeutronEndPosY", &EndPosY[NEUTRONS_CODE] );
    pTree->Branch( "PionEndPosY", &EndPosY[CH_PIONS_CODE] );
    pTree->Branch( "Pi0EndPosY", &EndPosY[NE_PIONS_CODE] );
    pTree->Branch( "MesonEndPosY", &EndPosY[MESONS_CODE] );
    pTree->Branch( "BaryonEndPosY", &EndPosY[BARYONS_CODE] );
    pTree->Branch( "InDMEndPosY", &EndPosY[IN_DM_CODE] );
    pTree->Branch( "OutDMEndPosY", &EndPosY[OUT_DM_CODE] );
    pTree->Branch( "GENIEEndPosY", &EndPosY[GEN_PART_CODE] );

    pTree->Branch( "ProtonEndPosZ", &EndPosZ[PROTONS_CODE] );
    pTree->Branch( "NeutronEndPosZ", &EndPosZ[NEUTRONS_CODE] );
    pTree->Branch( "PionEndPosZ", &EndPosZ[CH_PIONS_CODE] );
    pTree->Branch( "Pi0EndPosZ", &EndPosZ[NE_PIONS_CODE] );
    pTree->Branch( "MesonEndPosZ", &EndPosZ[MESONS_CODE] );
    pTree->Branch( "BaryonEndPosZ", &EndPosZ[BARYONS_CODE] );
    pTree->Branch( "InDMEndPosZ", &EndPosZ[IN_DM_CODE] );
    pTree->Branch( "OutDMEndPosZ", &EndPosZ[OUT_DM_CODE] );
    pTree->Branch( "GENIEEndPosZ", &EndPosZ[GEN_PART_CODE] );
  */


}


int main( int argc, char ** argv ) {

    std::vector< std::string > Filenames;
    Filenames.push_back( argv[1] );

    std::string GenLabel = "bdm";
    art::InputTag MCTruthTag { GenLabel };

    TFile *fOut = new TFile( "BDM_MCAnalysis.root", "RECREATE" );
    TTree *fTree = new TTree( "MCParticles", "Primary MC Particles" );
   

    KinematicMap_t Px, Py, Pz, P, E, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass;
    CounterMap_t MultiplicityGen;


    InitTreeGen( fTree, MultiplicityGen, Px, Py, Pz, P, E, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass );
   
    int num_ev = 0;
    bool stop = true;
    for ( gallery::Event ev( Filenames ); !ev.atEnd() ; ev.next() ) {

        
        if(num_ev>4999) break;
        num_ev++;
        
        for ( auto& multPair: MultiplicityGen ) multPair.second = 0;
    
        std::cout << "Processing "
                  << "Run " << ev.eventAuxiliary().run() << ", "
                  << "Event " << ev.eventAuxiliary().event() << std::endl;

        auto const& MCTruthHandle = ev.getValidHandle< std::vector< simb::MCTruth > >( MCTruthTag );
        auto const& MCTruthObjs = *MCTruthHandle;

        ResetCounters(MultiplicityGen, Px, Py, Pz, P, E, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
        

        for ( size_t iMCTruth = 0; iMCTruth < MCTruthObjs.size(); ++iMCTruth ) {

            simb::MCTruth MCTruthObj = MCTruthObjs[iMCTruth];
            int nParticles = MCTruthObj.NParticles();

            TLorentzVector DMMomemtum;

            //We need incoming momentum DM first to compare with the other particles
            for ( int iParticle = 0; iParticle < nParticles; ++iParticle ) {

                const simb::MCParticle& MCParticleObjDMIn = MCTruthObj.GetParticle( iParticle );
                int pdgCode = MCParticleObjDMIn.PdgCode();

                if ( pdgCode == OUT_DM_CODE && MCParticleObjDMIn.StatusCode() == 0 ) {

                    DMMomemtum = MCParticleObjDMIn.Momentum(0); 
                    SaveInfo( MCParticleObjDMIn, IN_DM_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
                    
                }
            }

            for ( int iPart = 0; iPart < nParticles; ++iPart ) {

                const simb::MCParticle& MCParticleObj = MCTruthObj.GetParticle( iPart );
                int pdgCode = MCParticleObj.PdgCode();
                
                // Flags for Visible particles
                //bool Meson_flag, Baryon_flag, Proton_flag = false;

              //if(pdgCode == IN_DM_CODE || pdgCode == OUT_DM_CODE){  // CHECKING THINGS...
              /*  if(std::abs(pdgCode) == PROTONS_CODE){
                    std::cout << "test particle found! PdgCode: "<< pdgCode << " status code: " << MCParticleObj.StatusCode() << std::endl;   
                    std::cout << "test Initial Energy: " << MCParticleObj.E(0) << std::endl;
                    std::cout << "test Final Energy: " << MCParticleObj.E(MCParticleObj.NumberTrajectoryPoints()-1) << std::endl;
                    std::cout << "test Initial Postion X: " << MCParticleObj.Vx(0) << std::endl;
                    std::cout << "test Final Postion X:   " << MCParticleObj.Vx(MCParticleObj.NumberTrajectoryPoints()-1) << std::endl;
                    //std::cout << "test Initial Postion Z: " << MCParticleObj.Vz(0) << std::endl;
                    std::cout << std::endl;
                }*/
               // std::cout << "DM Px: " << DMMomemtum.Px() << std::endl;
                SaveInfo( MCParticleObj, ALL_PART_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);

                if ( pdgCode > GEN_PART_CODE && pdgCode != OUT_DM_CODE ) { 
                    
                    SaveInfo( MCParticleObj, GEN_PART_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);

                } else if ( abs(pdgCode) > MESONS_CODE && abs(pdgCode) < 400 || abs(pdgCode) == 130) {

                    SaveInfo( MCParticleObj, MESONS_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
                    
                    SaveInfo( MCParticleObj, VISIBLE_PART_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
                
                } else if ( abs(pdgCode) > BARYONS_CODE && abs(pdgCode) < 4000 ) {

                   SaveInfo( MCParticleObj, BARYONS_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);

                   SaveInfo( MCParticleObj, VISIBLE_PART_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);

                } else if ( pdgCode == OUT_DM_CODE && MCParticleObj.StatusCode() == 1 ) {

                    SaveInfo( MCParticleObj, OUT_DM_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
                   
                } else if ( abs(pdgCode) == PROTONS_CODE ) {

                    SaveInfo( MCParticleObj, PROTONS_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);

                    SaveInfo( MCParticleObj, VISIBLE_PART_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);

                } else if ( abs(pdgCode) == NEUTRONS_CODE ) {

                    SaveInfo( MCParticleObj, NEUTRONS_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
                
                } else if ( abs(pdgCode) == CH_PIONS_CODE ) {

                    SaveInfo( MCParticleObj, CH_PIONS_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
                
                    SaveInfo( MCParticleObj, VISIBLE_PART_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);

                } else if ( abs(pdgCode) == NE_PIONS_CODE ) {

                    SaveInfo( MCParticleObj, NE_PIONS_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
               
                    SaveInfo( MCParticleObj, VISIBLE_PART_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);

                } else if ( abs(pdgCode) == ARGON40_CODE ) {

                    SaveInfo( MCParticleObj, ARGON40_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
               
                } else if ( abs(pdgCode) == ARGON39_CODE ) {

                    SaveInfo( MCParticleObj, ARGON39_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
             
                } else if ( abs(pdgCode) == CL39_CODE ) {

                    SaveInfo( MCParticleObj, CL39_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
                
                } else if ( abs(pdgCode) == GEN_PART_CODE ) {

                    SaveInfo( MCParticleObj, GEN_PART_CODE, MultiplicityGen, Px,  Py, Pz, P, E, DMMomemtum, Angle, InitPosX, InitPosY, InitPosZ, /*EndPosX, EndPosY, EndPosZ,*/ Mass);
                    
                }

            } // Loop over MCParticles for each MCTruth
            
           
        } // Loop over MCTruth 
        
        fTree->Fill();
        
    } // End of an event

    fOut->Write(0,TObject::kOverwrite);
    return 0;
}
