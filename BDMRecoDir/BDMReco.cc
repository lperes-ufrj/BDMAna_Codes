/**
 * @file BDMReco.cc
 * @author Leonardo Peres (leoperes@pos.if.ufrj.br) and Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @brief  Analyzer of Reconstructed information of Boosted Dark Matter
 * @version 0.4
 * @date 2021-04-07
 * 
 * 
 * @copyright Copyright (c) 2021
 * 
 */

/*
 * Main changes:
 * 
 *  * tree data moved into its own structure: relying on a long list of
 *    parameters in a function is error-prone, i.e. it's not trivial to keep
 *    the correct order in the function declaration as in the function call;
 *    the recommended approach is instead purely based on names rather than
 *    positions; how did I get that:
 *      * created a new `struct` data type
 *      * moved the variable definition from the main program into that struct
 *      * moved the two functions (`InitTreeGen()` and `ResetCounters()`)
 *        as member function of that struct; since the names matched, I did
 *        not change anything except than removing all the variables from
 *        the parameter list: the two functions will now always operate on the
 *        variable within the new struct
 *        (exception was three `ResetCounters()` parameter names)
 * 
 * Minor changes:
 *  * replaced `auto` with the actual data type in the function declarations;
 *    `auto` in there is a C++20 feature that may be not supported
 *  * replaced `.resize(0)` with `.clear()` (equivalent, but semantically more
 *    expressive)
 *  * `InitTreeGen()` now "returns" `void` (it was not returning anything
 *    already but it was declared to return a pointer to something)
 *  * for start position and momentum, replaced
 *    `TVector3 p_tmp; p_tmp.SetXYZ(p.X(), p.Y(), p.Z()); p_list.push_back(p_tmp);`
 *    style of statements with the single
 *    `p_list.emplace_back(p.X(), p.Y(), p.Z());`
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
#include <array>

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

//#include "fhiclcpp/ParameterSet.h"

// LArSoft, nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "test/Geometry/geometry_unit_test_base.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

//Geometry includes
#include "larcore/Geometry/Geometry.h"
#include "test/Geometry/geometry_unit_test_base.h"
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// DUNE libraries
#include "dune/Geometry/DuneApaChannelMapAlg.h"
#include "dune/Geometry/GeoObjectSorterAPA.h"

// PARTICLES PDG CODES

#define OUT_DM_CODE 2000010000

//===========================================================================

struct DUNEGeometry10kt_1x2x6_EnvironmentConfiguration
  : public testing::BasicGeometryEnvironmentConfiguration<geo::DuneApaChannelMapAlg>
{
  using base_t
    = testing::BasicGeometryEnvironmentConfiguration<geo::DuneApaChannelMapAlg>;
  
  /// Default constructor
  DUNEGeometry10kt_1x2x6_EnvironmentConfiguration() { DefaultInit(); }
  
  /// Constructor; accepts the name as parameter
  DUNEGeometry10kt_1x2x6_EnvironmentConfiguration(std::string name):
    DUNEGeometry10kt_1x2x6_EnvironmentConfiguration()
    { base_t::SetApplicationName(name); }
  
  
    private:
  void DefaultInit()
    {
      // overwrite the configuration that happened in the base class:
      base_t::SetApplicationName("DUNEapp");
      base_t::SetDefaultGeometryConfiguration(R"(
        SurfaceY:        147828 # underground, 4850 feet to cm (DUNE DocDB-3833)
        Name:             "dune10kt_v1_1x2x6"
        GDML:             "dune10kt_v1_1x2x6.gdml"
        ROOT:             "dune10kt_v1_1x2x6.gdml"
        SortingParameters: { ChannelsPerOpDet: 1 }
        DisableWiresInG4:  true
        )");
    }
}; // class DUNEGeometry10kt_1x2x6_EnvironmentConfiguration

class DUNE10kt_1x2x6_Services
  : public testing::GeometryTesterEnvironment<DUNEGeometry10kt_1x2x6_EnvironmentConfiguration>
{
    public:
  using Config_t = DUNEGeometry10kt_1x2x6_EnvironmentConfiguration;
  using Base_t = testing::GeometryTesterEnvironment<Config_t>;

  DUNE10kt_1x2x6_Services(Config_t config): Base_t(std::move(config), false)
    { Setup(); }
  
    protected:
  virtual std::unique_ptr<geo::GeometryCore> CreateNewGeometry() const override
    {
      /*
       * rewritten because the initialization of the channel mapping is custom
       */
      std::string ProviderParameterSetPath
        = this->Config().GeometryParameterSetPath();

      //
      // create the new geometry service provider
      //
      fhicl::ParameterSet ProviderConfig
        = this->Parameters().template get<fhicl::ParameterSet>
          (ProviderParameterSetPath);
      auto new_geom = std::make_unique<geo::GeometryCore>(ProviderConfig);

      std::string RelativePath
        = ProviderConfig.get< std::string>("RelativePath", "");

      std::string
        GDMLFileName = RelativePath + ProviderConfig.get<std::string>("GDML"),
        ROOTFileName = RelativePath + ProviderConfig.get<std::string>("ROOT");

      // Search all reasonable locations for the geometry file;
      // we see if by any chance art's FW_SEARCH_PATH directory is set and try
      // there;
      // if not, we do expect the path to be complete enough for ROOT to cope.
      cet::search_path sp("FW_SEARCH_PATH");

      std::string ROOTfile;
      if (!sp.find_file(ROOTFileName, ROOTfile)) ROOTfile = ROOTFileName;

      // we really don't care of GDML file, since we are not going to run Geant4
      std::string GDMLfile;
      if (!sp.find_file(GDMLFileName, GDMLfile)) {
        mf::LogWarning("CreateNewGeometry") << "GDML file '"
          << GDMLfile << "' not found.";
      }

      // initialize the geometry with the files we have found
      new_geom->LoadGeometryFile(GDMLfile, ROOTfile);

      //
      // create the new channel map
      //
      fhicl::ParameterSet SortingParameters
        = ProviderConfig.get<fhicl::ParameterSet>("SortingParameters", {});
      auto pChannelMap
        = std::make_unique<geo::DuneApaChannelMapAlg>(SortingParameters);
      pChannelMap->setSorter(*new geo::GeoObjectSorterAPA(SortingParameters));

      // connect the channel map with the geometry, that shares ownsership
      // (we give up ours at the end of this method)
      new_geom->ApplyChannelMap(std::move(pChannelMap));

      return new_geom;
    } // CreateNewGeometry()
  
}; // class DUNE10kt_1x2x6_Services


//===========================================================================

using CounterMap_t = std::map<int, unsigned int>;
using KinematicMap_t = std::map<int, std::vector<double>>;

struct TreeData
{

    int inteventno, intrunno, intsubrunno, nPrimaries, nPFParticle, nPFParticleVtxDaughters, nPFParticleTrkDaughters, nPFParticleSwrDaughters;
    std::vector<double> DaughterTrackLengths, DaughterStartMomentumTrack, DaughterEndMomentumTrack, PandoraNuIDs, SwrID, DaughterTrackKE, DaughterTrackRange, DaughterTrackPitch, DaughterTrackMomentumIfMuon, DaughterTrackMomentumIfProton, DaughterMomentumMultiScatter,InDMMomentum, InDMEnergy, OutDMMomentum, OutDMEnergy;
    std::vector<int> best_plane_pid, nTrk_Cal, nTrk_Pid, PFParticleID, VtxStatus, VtxID, PFP_PdgCode, PFP_Parent, intCCNC, nDaughters, Mode, InteractionType, Swrbest_plane;
    std::vector<TVector3> VtxXYZ, DaughterVertexTrack, SwrDirection, SwrDirectionErr, SwrShowerStart, SwrShowerStartErr;
    std::vector<TVector3> DaughterStartPoint, DaughterStartDirection, SunDirectionFromTrueBDM, DaughterEndPoint, DaughterMultiScatterStartingPoint;
    std::vector<std::vector<double>> SwrEnergy, chi2proton;
    std::vector<std::vector<float>> DaughterTrackdEdx, DaughterTrackResidualRange;
    std::vector<std::vector<int>> track_PID_pdg;
    std::vector<bool> track_isContained, IsVtxPrimary, IsVtxDaughter;
    std::vector<TLorentzVector> PrimaryBDMVertex;
    TVector3 TotalMomEvent;

    void InitTreeGen(TTree *pTree);
    void ResetCounters();

}; // struct TreeData

void TreeData::InitTreeGen(TTree *pTree)
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
    pTree->Branch("IsVtxPrimary", &IsVtxPrimary);
    pTree->Branch("IsVtxDaughter", &IsVtxDaughter);

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

    pTree->Branch("Mode", &Mode);
    pTree->Branch("InteractionType", &InteractionType);
    pTree->Branch("DaughterStartPoint", &DaughterStartPoint);
    pTree->Branch("DaughterEndPoint", &DaughterEndPoint);
    pTree->Branch("DaughterStartDirection", &DaughterStartDirection);

    pTree->Branch("SwrEnergy", &SwrEnergy);
    pTree->Branch("TotalMomEvent", &TotalMomEvent);
    pTree->Branch("SunDirectionFromTrueBDM", &SunDirectionFromTrueBDM);
    pTree->Branch("Swrbest_plane", &Swrbest_plane);

    pTree->Branch("DaughterTrackdEdx", &DaughterTrackdEdx);
    pTree->Branch("DaughterTrackResidualRange", &DaughterTrackResidualRange);

    pTree->Branch("best_plane_pid", &best_plane_pid);
    pTree->Branch("nTrk_Cal", &nTrk_Cal);
    pTree->Branch("nTrk_Pid", &nTrk_Pid);

    pTree->Branch("chi2proton", &chi2proton);
    pTree->Branch("track_PID_pdg", &track_PID_pdg);
    pTree->Branch("track_isContained", &track_isContained);
    pTree->Branch("DaughterTrackKE", &DaughterTrackKE);
    pTree->Branch("DaughterTrackRange", &DaughterTrackRange);
    pTree->Branch("DaughterTrackPitch", &DaughterTrackPitch);

    pTree->Branch("DaughterTrackMomentumIfMuon", &DaughterTrackMomentumIfMuon);
    pTree->Branch("DaughterTrackMomentumIfProton", &DaughterTrackMomentumIfProton);
    pTree->Branch("DaughterMomentumMultiScatter", &DaughterMomentumMultiScatter);
    pTree->Branch("DaughterMultiScatterStartingPoint", &DaughterMultiScatterStartingPoint);
    
    pTree->Branch("PrimaryBDMVertex", &PrimaryBDMVertex);
    pTree->Branch("InDMMomentum", &InDMMomentum);
    pTree->Branch("InDMEnergy", &InDMEnergy);
    pTree->Branch("OutDMMomentum", &OutDMMomentum);
    pTree->Branch("OutDMEnergy", &OutDMEnergy);



} // TreeData::InitTreeGen()

void TreeData::ResetCounters()
{

    inteventno = 0;
    intrunno = 0;
    intrunno = 0;
    nPrimaries = 0;
    nDaughters.clear();
    PFP_PdgCode.clear();
    PFP_Parent.clear();
    PFParticleID.clear();
    nPFParticle = 0;
    nPFParticleVtxDaughters = 0;
    nPFParticleTrkDaughters = 0;
    nPFParticleSwrDaughters = 0;
    VtxStatus.clear();
    VtxID.clear();
    VtxXYZ.clear();
    DaughterTrackLengths.clear();
    DaughterVertexTrack.clear();
    DaughterStartMomentumTrack.clear();
    DaughterEndMomentumTrack.clear();
    intCCNC.clear();
    SwrID.clear();
    SwrDirection.clear();
    SwrDirectionErr.clear();
    SwrShowerStart.clear();
    SwrShowerStartErr.clear();
    PandoraNuIDs.clear();
    Mode.clear();
    InteractionType.clear();
    DaughterStartPoint.clear();
    DaughterStartDirection.clear();
    SwrEnergy.clear();
    TotalMomEvent.SetXYZ(0, 0, 0);
    SunDirectionFromTrueBDM.clear();
    Swrbest_plane.clear();
    DaughterTrackdEdx.clear();
    DaughterTrackResidualRange.clear();
    best_plane_pid.clear();
    nTrk_Cal.clear();
    nTrk_Pid.clear();
    track_PID_pdg.clear();
    track_isContained.clear();
    DaughterTrackKE.clear();
    DaughterTrackRange.clear();
    DaughterTrackPitch.clear();
    chi2proton.clear();
    DaughterEndPoint.clear();
    IsVtxPrimary.clear();
    IsVtxDaughter.clear();
    DaughterTrackMomentumIfMuon.clear();
    DaughterTrackMomentumIfProton.clear();
    DaughterMomentumMultiScatter.clear();
    DaughterMultiScatterStartingPoint.clear();
    PrimaryBDMVertex.clear();
    InDMMomentum.clear();
    InDMEnergy.clear();
    OutDMMomentum.clear();
    OutDMEnergy.clear();

} // TreeData::ResetCounters()

//Detector Limits =========================
float fFidVolXmin = 0;
float fFidVolXmax = 0;
float fFidVolYmin = 0;
float fFidVolYmax = 0;
float fFidVolZmin = 0;
float fFidVolZmax = 0;

bool insideFV(geo::Point_t const &vertex)
{

    double const x = vertex.X();
    double const y = vertex.Y();
    double const z = vertex.Z();

    return x > fFidVolXmin && x < fFidVolXmax &&
           y > fFidVolYmin && y < fFidVolYmax &&
           z > fFidVolZmin && z < fFidVolZmax;

} // insideFV()

//geo::GeometryCore const* geom;
//lar::providerFrom<geo::Geometry>() = geom ;

void GeoLimits(geo::GeometryCore const& geom,float fFidVolCutX, float fFidVolCutY, float fFidVolCutZ ) {

    // Define histogram boundaries (cm).
    // For now only draw cryostat=0.
    double minx = 1e9;
    double maxx = -1e9;
    double miny = 1e9;
    double maxy = -1e9;
    double minz = 1e9;
    double maxz = -1e9;
    for (size_t i = 0; i < geom.NTPC(); ++i)
    {
        double local[3] = {0., 0., 0.};
        double world[3] = {0., 0., 0.};
        const geo::TPCGeo &tpc = geom.TPC(i);
        tpc.LocalToWorld(local, world);
        if (minx > world[0] - geom.DetHalfWidth(i))
            minx = world[0] - geom.DetHalfWidth(i);
        if (maxx < world[0] + geom.DetHalfWidth(i))
            maxx = world[0] + geom.DetHalfWidth(i);
        if (miny > world[1] - geom.DetHalfHeight(i))
            miny = world[1] - geom.DetHalfHeight(i);
        if (maxy < world[1] + geom.DetHalfHeight(i))
            maxy = world[1] + geom.DetHalfHeight(i);
        if (minz > world[2] - geom.DetLength(i) / 2.)
            minz = world[2] - geom.DetLength(i) / 2.;
        if (maxz < world[2] + geom.DetLength(i) / 2.)
            maxz = world[2] + geom.DetLength(i) / 2.;
    }

    fFidVolXmin = minx + fFidVolCutX;
    fFidVolXmax = maxx - fFidVolCutX;
    fFidVolYmin = miny + fFidVolCutY;
    fFidVolYmax = maxy - fFidVolCutY;
    fFidVolZmin = minz + fFidVolCutZ;
    fFidVolZmax = maxz - fFidVolCutZ;
} // GeoLimits()


int main(int argc, char **argv)
{
    
    std::vector<std::string> filenames;
    std::string const InputSpec = argv[1];
    if (InputSpec.length() > 5U && InputSpec.substr(InputSpec.length() - 5U) == ".root") {
      filenames.push_back(InputSpec);
    }
    else {
      //std::string ListOfFiles;
      std::ifstream ListOfFiles(InputSpec);

      std::copy(std::istream_iterator<std::string>(ListOfFiles),
                std::istream_iterator<std::string>(),
                std::back_inserter(filenames));
    }

    for (size_t i = 0; i < filenames.size(); i++)
    {
        std::cout << "file name at " << i << " is " << filenames[i] << std::endl;
    }

    std::string PandoraLabel = "pandora";
    art::InputTag fPandoraLabel{PandoraLabel};

    std::string GeantLabel = "largeant";
    art::InputTag fGeantLabel{GeantLabel};

    std::string GeneratorLabel = "bdm";
    art::InputTag fGeneratorLabel{GeneratorLabel};

    art::InputTag const fPandoraTrackModuleLabel = {
        fPandoraLabel.label() + "Track",
        fPandoraLabel.instance(),
        fPandoraLabel.process()};

    art::InputTag const fPandoraShowerModuleLabel = {
        fPandoraLabel.label() + "Shower",
        fPandoraLabel.instance(),
        fPandoraLabel.process()};

    art::InputTag const fPandoraCaloModuleLabel = {
        fPandoraLabel.label() + "calo",
        fPandoraLabel.instance(),
        fPandoraLabel.process()};

    art::InputTag const fPandoraPIDModuleLabel = {
        fPandoraLabel.label() + "pid",
        fPandoraLabel.instance(),
        fPandoraLabel.process()};

    TFile *fOut = new TFile("BDM_RecoParticles.root", "RECREATE");
    TTree *fTree = new TTree("RecoParticles", " Reco Info from BDM sample");

    // summon the geometry, using hard-coded configuration (`DefaultInit()`):
    DUNEGeometry10kt_1x2x6_EnvironmentConfiguration config("BDMreco");
    DUNE10kt_1x2x6_Services DUNEenvironment(config);

    geo::GeometryCore const& geom
      = *(DUNEenvironment.Provider<geo::GeometryCore>());
    
    double XYZ[3] = {0.0, 0.0, 0.0};
    // double DaughterStartPoint_tmp[3] = {0.0, 0.0, 0.0};
    // double DaughterStartDirection_tmp[3] = {0.0, 0.0, 0.0};
    int n_Events = 0;

    // VARIABLES TO SAVE IN THE TREES ===============================================================================
    TreeData td;
    td.InitTreeGen(fTree);

    trkf::TrackMomentumCalculator trkm;

    GeoLimits(geom, 10, 10, 10); // GET DETECTOR GEOMETRY LIMITS exclude 10cm in each axis===============================================================

    for (gallery::Event ev(filenames); !ev.atEnd(); ev.next())
    {
        n_Events++;

        //double TotalMomX = 0;
        //double TotalMomY = 0;
        //double TotalMomZ = 0;

        td.ResetCounters();

        auto const &evaux = ev.eventAuxiliary();
        td.intrunno = evaux.run();
        td.intsubrunno = evaux.subRun();
        td.inteventno = evaux.event();

        // EVENT INFO ============================================================================

        std::cout << "Processing "
                  << "Run " << ev.eventAuxiliary().run() << ", "
                  << "Event " << ev.eventAuxiliary().event() << std::endl;

        // TRUE INFORMATION =======================================================================

        auto const &MCParticleHandle = ev.getValidHandle<std::vector<simb::MCParticle>>(fGeantLabel);
        auto const &MCParticleObjs = *MCParticleHandle;

        auto const &MCTruthHandle = ev.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorLabel);
        auto const &MCTruthObjs = *MCTruthHandle;

        // SAVE TRUE INTERACTION INFOR (GENERATOR) FROM SIGNAL  // 0=CC 1=NC =================================

        for (size_t i = 0; i < MCTruthObjs.size(); i++)
        {
            simb::MCTruth MCTruthObj = MCTruthObjs[i];
            simb::MCNeutrino const &thisNu = MCTruthObj.GetNeutrino();
            td.intCCNC.push_back(thisNu.CCNC()); // 0=CC 1=NC
            td.Mode.push_back(thisNu.Mode());
            td.InteractionType.push_back(thisNu.InteractionType());
            //std::cout << "thisNu.CCNC() = " << thisNu.CCNC() << std::endl;

            int nParticles = MCTruthObj.NParticles();

            TLorentzVector DMMomentum;
            TLorentzVector DMInteracPosition;

            for (int iParticle = 0; iParticle < nParticles; ++iParticle)
            {

                const simb::MCParticle &MCParticleObjDMIn = MCTruthObj.GetParticle(iParticle);
                int pdgCode = MCParticleObjDMIn.PdgCode();

                if (pdgCode == OUT_DM_CODE && MCParticleObjDMIn.StatusCode() == 0)
                {
                    DMMomentum = MCParticleObjDMIn.Momentum(0);
                    DMInteracPosition = MCParticleObjDMIn.Position(0);
                    td.SunDirectionFromTrueBDM.push_back(DMMomentum.Vect().Unit()); // SAVE SUN DIRECTION !!!
                    td.PrimaryBDMVertex.push_back(DMInteracPosition);
                    td.InDMMomentum.push_back(MCParticleObjDMIn.P(0));
                    td.InDMEnergy.push_back(MCParticleObjDMIn.E(0));

                } else if (pdgCode == OUT_DM_CODE && (MCParticleObjDMIn.StatusCode() == 1 || MCParticleObjDMIn.StatusCode() == 15 ) )
                {
                    td.OutDMMomentum.push_back(MCParticleObjDMIn.P(0));
                    td.OutDMEnergy.push_back(MCParticleObjDMIn.E(0));                 
                }
                
            }
        }
   
        
        // RECO INFORMATION ==================================================================

        auto const &TrackListHandle = ev.getValidHandle<std::vector<recob::Track>>(fPandoraTrackModuleLabel);
        auto const &PandoraTrackObjs = *TrackListHandle;

        auto const &PfparticleHandle = ev.getValidHandle<std::vector<recob::PFParticle>>(fPandoraLabel);
        auto const &PandoraPfparticleObjs = *PfparticleHandle;

        auto const &VertexHandle = ev.getValidHandle<std::vector<recob::Vertex>>(fPandoraLabel);
        auto const &PandoraVertexObjs = *VertexHandle;

        
        // ASSOCIATIONS PFP+TRACK, PFP+SHOWER, PFP+VERTEX, TRACK+CALO, TRACK+PID ==================
        art::FindManyP<recob::Track> PandoraPfp_trk(PfparticleHandle, ev, fPandoraTrackModuleLabel);
        art::FindManyP<recob::Shower> PandoraPfp_Swr(PfparticleHandle, ev, fPandoraShowerModuleLabel);
        art::FindManyP<recob::Vertex> PandoraPfp_vtx(PfparticleHandle, ev, PandoraLabel);
        art::FindManyP<anab::Calorimetry> PandoraTrk_Cal(TrackListHandle, ev, fPandoraCaloModuleLabel);
        art::FindManyP<anab::ParticleID> PandoraTrk_Pid(TrackListHandle, ev, fPandoraPIDModuleLabel);

        // SAVE ONLY PFPARTICLE INFORMATION ================================

        td.nPFParticle = PandoraPfparticleObjs.size();
        //int neutrinoID = 99999;

        for (size_t iPfp = 0; iPfp < PandoraPfparticleObjs.size(); iPfp++)
        {

            td.PFP_PdgCode.push_back(PandoraPfparticleObjs[iPfp].PdgCode());
            if (PandoraPfparticleObjs[iPfp].IsPrimary())
                td.nPrimaries++;
            td.nDaughters.push_back(PandoraPfparticleObjs[iPfp].NumDaughters());
            td.PFParticleID.push_back(PandoraPfparticleObjs[iPfp].Self());
            td.PFP_Parent.push_back(PandoraPfparticleObjs[iPfp].Parent());
            if (PandoraPfparticleObjs[iPfp].IsPrimary() && (std::abs(PandoraPfparticleObjs[iPfp].PdgCode()) == 14 || std::abs(PandoraPfparticleObjs[iPfp].PdgCode()) == 12))
                td.PandoraNuIDs.push_back(PandoraPfparticleObjs[iPfp].Self());
        }

        //std::cout << "PandoraNuIDs.size() = " << PandoraNuIDs.size() << std::endl;

        if (td.PandoraNuIDs.size() == 0)
            continue;

        //FROM NOW ON, ONLY SAVE DAUGHTERS INFORMATION OF PANDORA NEUTRINOS PRIMARIES

        for (size_t iPfp = 0; iPfp < PandoraPfparticleObjs.size(); iPfp++)
        {

            //Check if it is a primary vertex or a daughter vertex
            if (PandoraPfparticleObjs[iPfp].Self() == td.PandoraNuIDs[0])
            {
                td.IsVtxPrimary.push_back(true);
            }
            else
            {
                td.IsVtxPrimary.push_back(false);
            }
            if (PandoraPfparticleObjs[iPfp].Parent() == td.PandoraNuIDs[0])
            {
                td.IsVtxDaughter.push_back(true);
            }
            else
            {
                td.IsVtxDaughter.push_back(false);
            }

            // VERTEXES RECO INFO ==================================================================

            auto const &Pfp_vtx = PandoraPfp_vtx.at(iPfp);

            td.nPFParticleVtxDaughters += Pfp_vtx.size();

            if (!Pfp_vtx.empty())
            {
                for (size_t i = 0; i < Pfp_vtx.size(); i++)
                {
                    auto const thisVtxPfp = Pfp_vtx[i];
                    td.VtxStatus.push_back(thisVtxPfp->status());
                    td.VtxID.push_back(thisVtxPfp->ID());
                    thisVtxPfp->XYZ(XYZ);
                    td.VtxXYZ.push_back(XYZ);
                }
            }

            if (PandoraPfparticleObjs[iPfp].Parent() != td.PandoraNuIDs[0]) // Only daughters from now on
                continue;

            // TRACKS RECO INFO ==================================================================

            auto const &Pfp_trk = PandoraPfp_trk.at(iPfp);

            td.nPFParticleTrkDaughters += Pfp_trk.size();

            if (!Pfp_trk.empty())
            {
                for (size_t i = 0; i < Pfp_trk.size(); i++)
                {
                    auto const thisTrkPfp = Pfp_trk[i];
                    td.DaughterTrackLengths.push_back(thisTrkPfp->Length());
                    td.DaughterStartMomentumTrack.push_back(thisTrkPfp->StartMomentum());
                    td.DaughterEndMomentumTrack.push_back(thisTrkPfp->EndMomentum());

                    td.DaughterStartPoint.emplace_back(thisTrkPfp->Start().X(), thisTrkPfp->Start().Y(), thisTrkPfp->Start().Z());
                    td.DaughterEndPoint.emplace_back(thisTrkPfp->End().X(), thisTrkPfp->End().Y(), thisTrkPfp->End().Z());
                    td.DaughterStartDirection.emplace_back(thisTrkPfp->StartDirection().X(), thisTrkPfp->StartDirection().Y(), thisTrkPfp->StartDirection().Z());         

                    double Length = thisTrkPfp->Length();
                    //std::cout << "Track Length is: " << thisTrkPfp->Length() << std::endl;
                    td.DaughterMomentumMultiScatter.push_back(trkm.GetMomentumMultiScatterChi2(thisTrkPfp));
                    td.DaughterMultiScatterStartingPoint.push_back(trkm.GetMultiScatterStartingPoint(thisTrkPfp));
                    td.DaughterTrackMomentumIfMuon.push_back(trkm.GetTrackMomentum(Length, 13));
                    td.DaughterTrackMomentumIfProton.push_back(trkm.GetTrackMomentum(Length, 2212));

                    // ASSOCIATION OF TRACKS WITH CALORIMETRY ===> PFP+TRK+CAL ===============================================

                    bool tmp_isContained = insideFV(thisTrkPfp->End());
                    td.track_isContained.push_back(tmp_isContained);

                    auto const &Trk_cal = PandoraTrk_Cal.at(i);
                    if (Trk_cal.empty())
                        continue;

                    //PIDAcal(Trk_cal, PIDA, )
                    td.nTrk_Cal.push_back(Trk_cal.size());
                    for (size_t p = 0; p < Trk_cal.size(); p++)
                    {
                        auto const thisCal = Trk_cal[p];
                        if (!thisCal->PlaneID().isValid)
                            continue;

                        int planenum = thisCal->PlaneID().Plane; // Get the plane number

                        if (planenum != 2)
                            continue; // Only look at the collection plane.

                        td.DaughterTrackdEdx.push_back(thisCal->dEdx());
                        td.DaughterTrackResidualRange.push_back(thisCal->ResidualRange());
                        td.DaughterTrackKE.push_back(thisCal->KineticEnergy());
                        td.DaughterTrackRange.push_back(thisCal->Range());
                        td.DaughterTrackPitch.push_back(thisCal->TrkPitchC());
                    }

                    auto const &Trk_Pid = PandoraTrk_Pid.at(i);

                    td.nTrk_Pid.push_back(Trk_Pid.size());

                    if (Trk_Pid.empty())
                        continue;

                    int ndfbest_plane = -1;
                    int best_plane = -1;
                    std::vector<int> tmp_pdg;
                    std::vector<double> tmp_chi2proton;

                    for (size_t p = 0; p < Trk_Pid.size(); p++)
                    {
                        auto const thisPid = Trk_Pid[p];

                        int ndfplane = thisPid->Ndf();

                        if (ndfplane > ndfbest_plane)
                        {
                            ndfbest_plane = ndfplane;
                            best_plane = p;
                        }
                    
                        tmp_pdg.push_back(thisPid->Pdg());
                        tmp_chi2proton.push_back(thisPid->Chi2Proton());
                        td.track_PID_pdg.push_back(tmp_pdg);
                        td.chi2proton.push_back(tmp_chi2proton);
                    }
                    
                    td.best_plane_pid.push_back(best_plane);
                    
                }
            }

            // SHOWERS RECO INFO ====================================================================

            auto const &Pfp_swr = PandoraPfp_Swr.at(iPfp);

            td.nPFParticleSwrDaughters += Pfp_swr.size();

            if (!Pfp_swr.empty())
            {
                for (size_t i = 0; i < Pfp_swr.size(); i++)
                {

                    auto const thisSwrPfp = Pfp_swr[i];

                    td.SwrID.push_back(thisSwrPfp->ID());
                    td.SwrDirection.push_back(thisSwrPfp->Direction());
                    td.SwrDirectionErr.push_back(thisSwrPfp->DirectionErr());
                    td.SwrShowerStart.push_back(thisSwrPfp->ShowerStart());
                    td.SwrShowerStartErr.push_back(thisSwrPfp->ShowerStartErr());
                    td.Swrbest_plane.push_back(thisSwrPfp->best_plane());
                    td.SwrEnergy.push_back(thisSwrPfp->Energy());

                }
            }
        }

        fTree->Fill();
    }
    std::cout << "The number of events analyzed was " << n_Events << "." << std::endl;
    fOut->Write();
    return 0;
}
