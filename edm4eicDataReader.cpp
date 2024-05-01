#define edm4eicDataReader_cpp

#include "fastjet/ClusterSequence.hh"

#include "edm4hep/MCParticleData.h"
#include "edm4eic/ReconstructedParticleData.h"
#include "edm4eic/MCRecoParticleAssociationData.h"
#include "edm4eic/InclusiveKinematicsData.h"
#include "podio/ObjectID.h"

#include "edm4eicDataReader.hh"

using genParticleContainer = ROOT::VecOps::RVec<edm4hep::MCParticleData>;
using recParticleContainer = ROOT::VecOps::RVec<edm4eic::ReconstructedParticleData>;
using mc2RecoMatching = ROOT::VecOps::RVec<edm4eic::MCRecoParticleAssociationData>;
using inclEventKinematics = ROOT::VecOps::RVec<edm4eic::InclusiveKinematicsData>;
using scatBeamElectronInfo = ROOT::VecOps::RVec<podio::ObjectID>;

void makeJets(unsigned int iSlot, scatBeamElectronInfo& _scat_beam_e, recParticleContainer& genParticles, recParticleContainer& recParticles, mc2RecoMatching& mc2Reco){
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);
    fastjet::Selector constituentSelector = fastjet::SelectorAbsRapMax(3.5) && fastjet::SelectorPtMin(conPtMin) && fastjet::SelectorPtMax(conPtMax);
    fastjet::Selector jetSelector = fastjet::SelectorAbsRapMax(2.5) && fastjet::SelectorEMin(5) && fastjet::SelectorPtMin(1);

    std::vector<fastjet::PseudoJet> recConstituents;
    for(int iRec = 0; iRec < recParticles.size(); iRec++){
        fastjet::PseudoJet recPartVec(recParticles[iRec].momentum.x, recParticles[iRec].momentum.y, recParticles[iRec].momentum.z, recParticles[iRec].energy);
        recPartVec.set_user_index(iRec);
        recConstituents.push_back(recPartVec);
    }
    fastjet::ClusterSequence csReco(constituentSelector(recConstituents), jetDef);
    std::vector<fastjet::PseudoJet> recoJets = jetSelector(csReco.inclusive_jets(1));

    std::vector<fastjet::PseudoJet> genConstituents;
    for(int iGen = 0; iGen < genParticles.size(); iGen++){
        //if(genParticles[iGen].generatorStatus != 1)continue;
        //double energy = std::sqrt(std::pow(genParticles[iGen].momentum.x, 2) + std::pow(genParticles[iGen].momentum.y, 2) + std::pow(genParticles[iGen].momentum.z, 2) + std::pow(genParticles[iGen].mass, 2));
        fastjet::PseudoJet genPartVec(genParticles[iGen].momentum.x, genParticles[iGen].momentum.y, genParticles[iGen].momentum.z, genParticles[iGen].energy);
        genPartVec.set_user_index(iGen);
        genConstituents.push_back(genPartVec);
    }
    fastjet::ClusterSequence csGen(constituentSelector(genConstituents), jetDef);
    std::vector<fastjet::PseudoJet> genJets = jetSelector(csGen.inclusive_jets(1));

    auto scatRecoElectron = recParticles[_scat_beam_e[0].index];
    if(scatRecoElectron.PDG != 11){
        logs[iSlot] << "makeJets: Scattered electron is not an electron ?! instead, PDG = "<<scatRecoElectron.PDG<< std::endl;
    }
    fastjet::PseudoJet recoScatElectronVec = fastjet::PseudoJet(scatRecoElectron.momentum.x, scatRecoElectron.momentum.y, scatRecoElectron.momentum.z, scatRecoElectron.energy);
    std::unordered_map<int, int> reco2GenJetMatchMap;
    int iReco = -1;
    for(auto& recoJet : recoJets){iReco++;
        if(getDeltaR(recoJet, recoScatElectronVec) < 1.0) continue;

        hist1DMaps.at("h1RecoPt" ).at(iSlot)->Fill(recoJet.pt());
        hist1DMaps.at("h1RecoEta").at(iSlot)->Fill(recoJet.eta());
        hist1DMaps.at("h1RecoPhi").at(iSlot)->Fill(recoJet.phi_std());
        for(auto& constit : recoJet.constituents()){
            hist2DMaps.at("h2RecoProfilevPt").at(iSlot)->Fill(recoJet.pt(), getDeltaR(constit, recoJet), constit.pt()/recoJet.pt());
        }
        hist2DMaps.at("h2RecoNConvPt").at(iSlot)->Fill(recoJet.pt(), recoJet.constituents().size());

        //Matching starts here
        int iMatched = -1;
        float minDR = 1000;
        int iGen = -1;
        for(auto& genJet : genJets){
            iGen++;
            double dR = getDeltaR(genJet, recoJet);
            if(dR > 0.2) continue;
            if(dR < minDR){
                minDR = dR;
                iMatched = iGen;
            }
        }if(iMatched == -1) continue;

        reco2GenJetMatchMap[iReco] = iMatched;
    }

    int iGen = -1;
    for(auto& genJet : genJets){iGen++;

        hist1DMaps.at("h1MCPt" ).at(iSlot)->Fill(genJet.pt());
        hist1DMaps.at("h1MCEta").at(iSlot)->Fill(genJet.eta());
        hist1DMaps.at("h1MCPhi").at(iSlot)->Fill(genJet.phi_std());
        for(auto& constit : genJet.constituents()){
            hist2DMaps.at("h2MCProfilevPt").at(iSlot)->Fill(genJet.pt(), getDeltaR(constit, genJet), constit.pt()/genJet.pt());
        }
        hist2DMaps.at("h2MCNConvPt").at(iSlot)->Fill(genJet.pt(), genJet.constituents().size());

        if(!doBijectiveMatching)continue;
        //Matching starts here
        int iMatched = -1;
        float minDR = 1000;
        int iReco = -1;
        for(auto& recoJet : recoJets){
            iReco++;
            double dR = getDeltaR(genJet, recoJet);
            if(dR > 0.2) continue;
            if(dR < minDR){
                minDR = dR;
                iMatched = iGen;
            }
        }if(iMatched == -1) continue;
        if(reco2GenJetMatchMap.find(iMatched) == reco2GenJetMatchMap.end()){
            reco2GenJetMatchMap.erase(iReco);
            continue;
        }else if(reco2GenJetMatchMap[iMatched] != iGen){
            reco2GenJetMatchMap.erase(iReco);
            continue;
        }
    }

    for(auto& match : reco2GenJetMatchMap){
        int _iReco = match.first;
        int _iGen = match.second;
        hist2DMaps.at("h2Pt").at(iSlot)->Fill(genJets[_iGen].pt(), recoJets[_iReco].pt());
        hist2DMaps.at("h2Eta").at(iSlot)->Fill(genJets[_iGen].eta(), recoJets[_iReco].eta());
        hist2DMaps.at("h2Phi").at(iSlot)->Fill(genJets[_iGen].phi_std(), recoJets[_iReco].phi_std());
        hist2DMaps.at("h2ResvPt").at(iSlot)->Fill(genJets[_iGen].pt(), (recoJets[_iReco].pt() - genJets[_iGen].pt())/genJets[_iGen].pt());
    }
    //hist2DMaps.at("h2NMatchvPt").at(iSlot)->Fill(recJetVec.Pt(), nMatched);
    //hist2DMaps.at("h2NMatchFracvPt").at(iSlot)->Fill(recJetVec.Pt(), nMatched/double(nRec));
}

void readJets(unsigned int iSlot, scatBeamElectronInfo& _scat_beam_e, recParticleContainer& genParticles, recParticleContainer& recParticles, mc2RecoMatching& mc2Reco, recParticleContainer& genJets, recParticleContainer& recJets) {
    //auto _histMap = histMap.GetAtSlot(iSlot);
    auto scatRecoElectron = recParticles[_scat_beam_e[0].index];
    if(scatRecoElectron.PDG != 11){
        logs[iSlot] << "readJets: Scattered electron is not an electron ?! instead, PDG = "<<scatRecoElectron.PDG<< std::endl;
    }
    ROOT::Math::PxPyPzEVector recoScatElectronVec(scatRecoElectron.momentum.x, scatRecoElectron.momentum.y, scatRecoElectron.momentum.z, scatRecoElectron.energy);
    //std::cout << recoScatElectronVec.Pt() <<" "<< recoScatElectronVec.Eta() << " " << recoScatElectronVec.Phi()<<" "<<recoScatElectron.PDG<< std::endl;
    for(int iReco = 0; iReco < recJets.size(); iReco++){
        if(recJets[iReco].type != 0) continue;
        if(recJets[iReco].energy < 5) continue;
        edm4eic::ReconstructedParticleData recJet = recJets[iReco];
        ROOT::Math::PxPyPzEVector recJetVec(recJet.momentum.x, recJet.momentum.y, recJet.momentum.z, recJet.energy);
        if(getDeltaR(recoScatElectronVec, recJetVec) < 1.0) continue;
        if(std::abs(recJetVec.Eta()) > 2.5) continue;
        int iMatched = -1;
        ROOT::Math::PxPyPzEVector matchedJetVec;
        float minDR = 1000;
        for(int iGen = 0; iGen < genJets.size(); iGen++){
            if(genJets[iGen].type != 0) continue;
            if(genJets[iGen].energy < 5) continue;
            edm4eic::ReconstructedParticleData genJet = genJets[iGen];
            ROOT::Math::PxPyPzEVector genJetVec(genJet.momentum.x, genJet.momentum.y, genJet.momentum.z, genJet.energy);
            if(std::abs(genJetVec.Eta()) > 2.5) continue;
            double dR = getDeltaR(genJetVec, recJetVec);
            if(dR > 0.2) continue;
            if(dR < minDR){
                minDR = dR;
                iMatched = iReco;
                matchedJetVec = genJetVec;
            }
        }if(iMatched == -1) continue;

       unsigned int nMatched = 0;
       for(int iMatch = 0; iMatch < mc2Reco.size(); iMatch++){
            int iMatchReco = mc2Reco[iMatch].recID;
            ROOT::Math::PxPyPzEVector recMatchVec(recParticles[iMatchReco].momentum.x, recParticles[iMatchReco].momentum.y, recParticles[iMatchReco].momentum.z, recParticles[iMatchReco].energy);
            if(getDeltaR(recMatchVec, recJetVec) < 1.0) nMatched++;
       }

       unsigned int nRec = 0;
        for(int iRec = 0; iRec < recParticles.size(); iRec++){
            ROOT::Math::PxPyPzEVector recPartVec(recParticles[iRec].momentum.x, recParticles[iRec].momentum.y, recParticles[iRec].momentum.z, recParticles[iRec].energy);
            if(getDeltaR(recPartVec, recJetVec) < 1.0) nRec++;      
        }

        hist2DMaps.at("h2Pt")    .at(iSlot)->Fill(matchedJetVec.Pt() , recJetVec.Pt());
        hist2DMaps.at("h2Eta")   .at(iSlot)->Fill(matchedJetVec.Eta(), recJetVec.Eta());
        hist2DMaps.at("h2Phi")   .at(iSlot)->Fill(matchedJetVec.Phi(), recJetVec.Phi());
        hist2DMaps.at("h2ResvPt").at(iSlot)->Fill(matchedJetVec.Pt(), (recJetVec.Pt() - matchedJetVec.Pt())/matchedJetVec.Pt());
        //hist2DMaps.at("h2NMatchvPt").at(iSlot)->Fill(recJetVec.Pt(), nMatched);
        //hist2DMaps.at("h2NMatchFracvPt").at(iSlot)->Fill(recJetVec.Pt(), nMatched/double(nRec));
    }
}

void readEDM4EICEvents(unsigned int iSlot){
    static std::string treeName = "events";
    static std::vector<std::string> columnNames = {"_InclusiveKinematicsSigma_scat", "GeneratedParticles", "ReconstructedParticles", "ReconstructedParticleAssociations", "GeneratedJets", "ReconstructedJets"};

    ROOT::RDataFrame reader(treeName, fileLists[iSlot]);

    auto _analyze = [iSlot](scatBeamElectronInfo& _scat_beam_e, recParticleContainer& genParticles, recParticleContainer& recParticles, mc2RecoMatching& mc2Reco, recParticleContainer& genJets, recParticleContainer& recJets) {
        static unsigned int nEvents = 0;
        if(!doClusterJets)readJets(iSlot, _scat_beam_e, genParticles, recParticles, mc2Reco, genJets, recJets);
        else makeJets(iSlot, _scat_beam_e, genParticles, recParticles, mc2Reco);
        if(++nEvents % 10000 == 0) logs[iSlot] << "Processed " << nEvents << " events" << std::endl;
    };

    reader.Foreach(_analyze, columnNames);
}

int main(int, char*[]) {
    ROOT::EnableThreadSafety();

    if(!doClusterJets){
        R = 1.0;
        conPtMin = 0.2;
        conPtMax = 100;
    }

    std::vector<std::thread> threads;
    //reader.Describe().Print();
    std::ifstream inFileList("/home/tanmaypani/eic/macros/testFiles.list");
    std::string fileName = "";
    if(!inFileList.is_open()){
        std::cerr << "Error: could not open file list" << std::endl;
        return 1;
    }
    unsigned long nFiles = 0;
    while(std::getline(inFileList, fileName)){
        unsigned int iSlot = nFiles % nThreads;
        if(!logs[iSlot].is_open())logs[iSlot].open("log_" + std::to_string(iSlot) + ".txt");
        logs[iSlot]<<"Adding file #"<<nFiles+1<<" : "<<fileName<<std::endl;
        fileLists[iSlot].push_back(fileName);
        nFiles++;
    }
        std::string cutsString = Form("R = %0.1f, %0.1f < p_{T, const.} < %0.1f", R, conPtMin, conPtMax);

    for(unsigned int iSlot = 0; iSlot < nThreads; iSlot++){
        std::cout << "Thread " << iSlot << " has " << fileLists[iSlot].size() << " files" << std::endl;

        hist1DMaps["h1MCPt" ].push_back(std::make_shared<TH1D>(Form("h1MCPt_%i" , iSlot), Form("h1MCPt , %s ;p_{T, jet}^{truth}(GeV/c); #frac{dN}{dp_{T}}", cutsString.c_str()), jetPtBins.size()-1, jetPtBins.data()));
        hist1DMaps["h1MCEta"].push_back(std::make_shared<TH1D>(Form("h1MCEta_%i", iSlot), Form("h1MCEta, %s ;#eta_{jet}^{truth}       ; #frac{dN}{d#eta} ", cutsString.c_str()), jetEtaBins.size()-1, jetEtaBins.data()));
        hist1DMaps["h1MCPhi"].push_back(std::make_shared<TH1D>(Form("h1MCPhi_%i", iSlot), Form("h1MCPhi, %s ;#phi_{jet}^{truth}       ; #frac{dN}{d#phi} ", cutsString.c_str()), jetPhiBins.size()-1, jetPhiBins.data()));
        hist2DMaps["h2MCProfilevPt"].push_back(std::make_shared<TH2D>(Form("h2MCProfilevPt_%i", iSlot), Form("h2MCProfilevPt, %s; p_{T, jet}^{truth}(GeV/c) ; #Delta R"    , cutsString.c_str()), jetPtBins.size()-1, jetPtBins.data(), unsigned(R/0.05), 0, R));
        hist2DMaps["h2MCNConvPt"   ].push_back(std::make_shared<TH2D>(Form("h2MCNConvPt_%i"   , iSlot), Form("h2MCNConvPt   , %s; p_{T, jet}^{truth}(GeV/c) ; N_{constit.}", cutsString.c_str()), jetPtBins.size()-1, jetPtBins.data(), 30, -0.5, 29.5)        );

        hist1DMaps["h1RecoPt" ].push_back(std::make_shared<TH1D>(Form("h1RecoPt_%i" , iSlot), Form("h1RecoPt , %s;p_{T, jet}^{reco}(GeV/c); #frac{dN}{dp_{T}}", cutsString.c_str()), jetPtBins.size()-1, jetPtBins.data()));
        hist1DMaps["h1RecoEta"].push_back(std::make_shared<TH1D>(Form("h1RecoEta_%i", iSlot), Form("h1RecoEta, %s;#eta_{jet}^{reco}       ; #frac{dN}{d#eta} ", cutsString.c_str()), jetEtaBins.size()-1, jetEtaBins.data()));
        hist1DMaps["h1RecoPhi"].push_back(std::make_shared<TH1D>(Form("h1RecoPhi_%i", iSlot), Form("h1RecoPhi, %s;#phi_{jet}^{reco}       ; #frac{dN}{d#phi} ", cutsString.c_str()), jetPhiBins.size()-1, jetPhiBins.data()));
        hist2DMaps["h2RecoProfilevPt"].push_back(std::make_shared<TH2D>(Form("h2RecoProfilevPt_%i", iSlot), Form("h2RecoProfilevPt, %s; p_{T, jet}^{reco}(GeV/c) ; #Delta R"    , cutsString.c_str()), jetPtBins.size()-1, jetPtBins.data(), unsigned(R/0.05), 0, R));
        hist2DMaps["h2RecoNConvPt"   ].push_back(std::make_shared<TH2D>(Form("h2RecoNConvPt_%i"   , iSlot), Form("h2RecoNConvPt   , %s; p_{T, jet}^{reco}(GeV/c) ; N_{constit.}", cutsString.c_str()), jetPtBins.size()-1, jetPtBins.data(), 30, -0.5, 29.5)        );

        hist2DMaps["h2Pt"           ].push_back(std::make_shared<TH2D>(Form("h2Pt_%i"           , iSlot), Form("h2Pt        , %s; p_{T, jet}^{truth}(GeV/c); p_{T, jet}^{reco}(GeV/c)"                          , cutsString.c_str()), 150, 0, 150, 150, 0, 150)           );
        hist2DMaps["h2Eta"          ].push_back(std::make_shared<TH2D>(Form("h2Eta_%i"          , iSlot), Form("h2Eta       , %s; #eta_{jet}^{truth}       ; #eta_{jet}^{reco}"                                 , cutsString.c_str()), 50, -2.5, 2.5, 50, -2.5, 2.5)       );
        hist2DMaps["h2Phi"          ].push_back(std::make_shared<TH2D>(Form("h2Phi_%i"          , iSlot), Form("h2Phi       , %s; #phi_{jet}^{truth}       ; #phi_{jet}^{reco}"                                 , cutsString.c_str()), 60, -3.14, 3.14, 60, -3.14, 3.14)   );
        hist2DMaps["h2ResvPt"       ].push_back(std::make_shared<TH2D>(Form("h2ResvPt_%i"       , iSlot), Form("h2ResvPt    , %s; p_{T, jet}^{truth}(GeV/c); #Delta(p_{T, jet})(reco, truth)/p_{T, jet}^{truth}", cutsString.c_str()), jetPtBins.size()-1, jetPtBins.data(), 20, -1, 1)            );

        //hist2DMaps["h2NMatchvPt"    ].push_back(std::make_shared<TH2D>(Form("h2NMatchvPt_%i"    , iSlot), Form("h2NMatchvPt, R = %0.1f, %0.1f < p_{T, const.} < %0.1f ;p_{T, jet}^{reco}(GeV/c);NMatchedRecoParticles", R, conPtMin, conPtMax), 150, 0, 150, 50, -0.5, 49.5)                       );
        //hist2DMaps["h2NMatchFracvPt"].push_back(std::make_shared<TH2D>(Form("h2NMatchFracvPt_%i", iSlot), Form("h2NMatchFracvPt, R = %0.1f, %0.1f < p_{T, const.} < %0.1f ;p_{T, jet}^{reco}(GeV/c);NMatchedRecoParticles/NRecoParticles", R, conPtMin, conPtMax), 150, 0, 150, 20, 0, 1)          );
        //hist2DMaps["h2Pt_clustered"        ].push_back(std::make_shared<TH2D>(Form("h2Pt_clustered_%i"       , iSlot), Form("h2Pt, R = %0.1f, %0.1f < p_{T, const.} < %0.1f ;p_{T, jet}^{truth}(GeV/c);p_{T, jet}^{reco}(GeV/c)", R, conPtMin, conPtMax), 150, 0, 150, 150, 0, 150)                      );
        //hist2DMaps["h2Eta_clustered"       ].push_back(std::make_shared<TH2D>(Form("h2Eta_clustered_%i"      , iSlot), Form("h2Eta, R = %0.1f, %0.1f < p_{T, const.} < %0.1f ;#eta_{jet}^{truth};#eta_{jet}^{reco}", R, conPtMin, conPtMax), 50, -2.5, 2.5, 50, -2.5, 2.5)                               );
        //hist2DMaps["h2Phi_clustered"       ].push_back(std::make_shared<TH2D>(Form("h2Phi_clustered_%i"      , iSlot), Form("h2Phi, R = %0.1f, %0.1f < p_{T, const.} < %0.1f ;#phi_{jet}^{truth};#phi_{jet}^{reco}", R, conPtMin, conPtMax), 60, -3.14, 3.14, 60, -3.14, 3.14)                           );
        //hist2DMaps["h2ResvPt_clustered"    ].push_back(std::make_shared<TH2D>(Form("h2ResvPt_clustered_%i"   , iSlot), Form("h2ResvPt, R = %0.1f, %0.1f < p_{T, const.} < %0.1f ;p_{T, jet}^{truth}(GeV/c);#Delta(p_{T, jet})/p_{T, jet}^{truth}", R, conPtMin, conPtMax), 150, 0, 150, 100, -1, 1)      );
        
        threads.emplace_back(readEDM4EICEvents, iSlot);
    }

    for(auto& thread : threads)thread.join();

    std::cout<<"Merging histograms ..."<<std::endl;
    for(auto& hist : hist1DMaps){
        TList list;
        unsigned int nTotEntries = 0;
        std::cout << "__Merging threads for: " << hist.first << std::endl;
        for(auto& h : hist.second){
            std::cout<<"____Entries in thread: "<<h->GetEntries()<<std::endl;
            nTotEntries += h->GetEntries(); 
            if(h == hist.second[0])continue;
            list.Add(h.get());
        }
        hist.second[0]->Merge(&list);
        hist.second[0]->SetName(hist.first.c_str());
        std::cout << "__Merged " << nTotEntries << " entries" << std::endl;
        std::cout << "__Entries in merged histogram: " << hist.second[0]->GetEntries() << std::endl;
        std::cout << "_________________________________________________________"<<std::endl;
    }

    for(auto& hist : hist2DMaps){
        TList list;
        unsigned int nTotEntries = 0;
        std::cout << "__Merging threads for: " << hist.first << std::endl;
        for(auto& h : hist.second){
            std::cout<<"____Entries in thread: "<<h->GetEntries()<<std::endl;
            nTotEntries += h->GetEntries(); 
            if(h == hist.second[0])continue;
            list.Add(h.get());
        }
        hist.second[0]->Merge(&list);
        hist.second[0]->SetName(hist.first.c_str());
        std::cout << "__Merged " << nTotEntries << " entries" << std::endl;
        std::cout << "__Entries in merged histogram: " << hist.second[0]->GetEntries() << std::endl;
        std::cout << "_________________________________________________________"<<std::endl;
    }
    
    std::string histFileName = Form("/home/tanmaypani/eic/macros/jets_edm4eic_R%0.1f_conPt%0.1f-%0.1f", R, conPtMin, conPtMax);
    if(doBijectiveMatching) histFileName += "_bijectiveMatching";
    histFileName += ".hist.root";

    std::cout << "Writing histograms to " << histFileName << std::endl;
    TFile *fout = TFile::Open(histFileName.c_str(), "RECREATE");
    for(auto& hist : hist1DMaps){
        hist.second[0]->Write();
    }
    for(auto& hist : hist2DMaps){
        hist.second[0]->Write();
    }
    fout->Close();
    hist2DMaps.clear();
    hist1DMaps.clear();
    std::cout << "Done!" << std::endl;
    return 0;
}