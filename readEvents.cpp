#include <ROOT/RDataFrame.hxx>
#include "Math/Vector4D.h"

#include <iostream>
#include <string>
#include <vector>

TH2D* h2Pt = new TH2D("h2Pt", "h2Pt", 150, 0, 150, 150, 0, 150);
TH2D* h2ResvPt = new TH2D("h2ResvPt", "h2ResvPt", 150, 0, 150, 100, -1, 1);
TH2D* h2Eta = new TH2D("h2Eta", "h2Eta", 50, -2.5, 2.5, 50, -2.5, 2.5);
TH2D* h2Phi = new TH2D("h2Phi", "h2Phi", 60, -3.14, 3.14, 60, -3.14, 3.14);

void analyze(ROOT::RVecI& _scat_e_reco, ROOT::RVecI& _scat_e_truth,
             ROOT::RVecF& genPx, ROOT::RVecF& genPy, ROOT::RVecF& genPz, ROOT::RVecI& genPDG, ROOT::RVecI& genStatus, ROOT::RVecF& genE,
             ROOT::RVecF& recPx, ROOT::RVecF& recPy, ROOT::RVecF& recPz, ROOT::RVecF& recE, ROOT::RVecI& recType, ROOT::RVecI& recPDG,
             ROOT::RVecU& simId, ROOT::RVecU& recId,
             ROOT::RVecF& genJetPx, ROOT::RVecF& genJetPy, ROOT::RVecF& genJetPz, ROOT::RVecF& genJetE, ROOT::RVecI& genJetType, ROOT::RVecI& genJetPDG,  
             ROOT::RVecF& recJetPx, ROOT::RVecF& recJetPy, ROOT::RVecF& recJetPz, ROOT::RVecF& recJetE, ROOT::RVecI& recJetType, ROOT::RVecI& recJetPDG){
    static unsigned int nEvents = 0;
    int iRecoScatElectron = _scat_e_reco[0];
    int iTruthScatElectron = _scat_e_truth[0];
    ROOT::Math::PxPyPzEVector recoScatElectron(recPx[iRecoScatElectron], recPy[iRecoScatElectron], recPz[iRecoScatElectron], recE[iRecoScatElectron]);
    ROOT::Math::PxPyPzEVector truthScatElectron(genPx[iTruthScatElectron], genPy[iTruthScatElectron], genPz[iTruthScatElectron], genE[iTruthScatElectron]);
    std::cout << truthScatElectron.Pt() <<" "<< truthScatElectron.Eta() << " " << truthScatElectron.Phi() <<" "<<genPDG[iTruthScatElectron]<< std::endl;
    //std::cout << recoScatElectron.Pt() <<" "<< recoScatElectron.Eta() << " " << recoScatElectron.Phi() <<" "<<recPDG[iRecoScatElectron]<< std::endl;
    for(unsigned int iRec = 0; iRec < recJetPx.size(); iRec++){
        if(recJetType[iRec] != 0)continue;
        if(recJetE[iRec] < 5)continue;
        //if(recPDG[iRec] != iGen)std::cout << "id mismatch" << std::endl;
        ROOT::Math::PxPyPzEVector recJetVec(recJetPx[iRec], recJetPy[iRec], recJetPz[iRec], recJetE[iRec]);
        if(ROOT::Math::VectorUtil::DeltaR(recoScatElectron, recJetVec) < 1.0) continue;
        if(std::abs(recJetVec.Eta()) > 2.5) continue;
        int iMatched = -1;
        ROOT::Math::PxPyPzEVector matchedJetVec;        
        float minDR = 1000;
        for(unsigned int iGen = 0; iGen < genPx.size(); iGen++){
            if(genJetType[iGen] != 0)continue;
            if(genJetE[iGen] < 5)continue;
            //if(genPDG[iGen] != iGen)std::cout << "id mismatch" << std::endl;
            ROOT::Math::PxPyPzEVector genJetVec(genJetPx[iGen], genJetPy[iGen], genJetPz[iGen], genJetE[iGen]);
            if(std::abs(genJetVec.Eta()) > 2.5) continue;
            double dR = ROOT::Math::VectorUtil::DeltaR(genJetVec, recJetVec);     
            if(dR > 0.2) continue;
            if(dR < minDR){
                minDR = dR;
                iMatched = iRec;
                matchedJetVec = recJetVec;
            }
        }if(iMatched == -1) continue;
        h2Pt->Fill (matchedJetVec.Pt() , recJetVec.Pt());
        h2Eta->Fill(matchedJetVec.Eta(), recJetVec.Pt());
        h2Phi->Fill(matchedJetVec.Phi(), recJetVec.Pt());
        h2ResvPt->Fill(matchedJetVec.Pt(), (recJetVec.Pt() - matchedJetVec.Pt())/matchedJetVec.Pt());
    }
    if(++nEvents % 10000 == 0) std::cout << "Processed " << nEvents << " events" << std::endl;
} 

void readEvents(){
    std::string input = "/home/tanmaypani/eic/data/campaign-23.12.0-craterlake/minQ2=1000/*.root";
    std::string treeName = "events";
    ROOT::RDataFrame reader(treeName, input);

    std::vector<std::string> allBranches = reader.GetColumnNames();

    std::vector<std::string> particleCols = {"MCParticles", "ReconstructedParticles"};
    std::vector<std::string> jetCols =  {"GeneratedJets", "ReconstructedJets"};
    std::vector<std::string> momCols = {"momentum.x", "momentum.y", "momentum.z", "energy", "type", "PDG"};
    std::vector<std::string> allCols = {};

    allCols.push_back("_InclusiveKinematicsSigma_scat.index");
    allCols.push_back("_InclusiveKinematicsTruth_scat.index");

    for(auto particleCol : particleCols){
        for(auto momCol : momCols){
            if(particleCol == "MCParticles"){
                if(momCol == "energy")continue;
                if(momCol == "type")continue;
            }
            allCols.push_back(particleCol + "." + momCol);
        }
        if(particleCol == "MCParticles"){
            allCols.push_back("MCParticles.generatorStatus");
            allCols.push_back("MCParticles.time");
        }
    }

    allCols.push_back("ReconstructedParticleAssociations.simID");
    allCols.push_back("ReconstructedParticleAssociations.recID");

    for(auto jetCol : jetCols){
        for(auto momCol : momCols){
            allCols.push_back(jetCol + "." + momCol);
        }
    }

    reader.Foreach(analyze, allCols);
          
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    gPad->SetLogz();
    h2Pt->Draw("colz");

    TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);
    gPad->SetLogz();
    h2Eta->Draw("colz");

    TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 800, 600);
    gPad->SetLogz();
    h2Phi->Draw("colz");

    TCanvas* canvas4 = new TCanvas("canvas4", "canvas4", 800, 600);
    gPad->SetLogz();
    h2ResvPt->Draw("colz");
}