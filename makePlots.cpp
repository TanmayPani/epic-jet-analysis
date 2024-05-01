std::unordered_map<std::string, TH1D*> mch1Map;
std::unordered_map<std::string, TH1D*> recoh1Map;

std::unordered_map<std::string, TH2D*> h2Map_old;
std::unordered_map<std::string, TH2D*> h2Map;
std::unordered_map<std::string, TH2D*> mch2Map;
std::unordered_map<std::string, TH2D*> recoh2Map;

void makePlots(){
    auto inFile = TFile::Open("jets_edm4eic_R1.0_conPt0.2-100.0_bijectiveMatching.hist.root", "READ");
    std::cout<<"Reading from file: "<<inFile->GetName()<<std::endl;
    TIter inFileKeys(inFile->GetListOfKeys());
    TKey* key;
    while((key = (TKey*)inFileKeys())){
        auto obj = key->ReadObj();
        std::string objName = obj->GetName();
        std::cout<<"_____Reading object: "<<objName<<std::endl;
        if(obj->IsA()->InheritsFrom(TH1D::Class())){
            if(objName.find("h1Reco") != std::string::npos)recoh1Map[objName] = (TH1D*)obj;
            else if(objName.find("h1MC") != std::string::npos)mch1Map[objName] = (TH1D*)obj;
        }
        else if(obj->IsA()->InheritsFrom(TH2D::Class())){
            if(objName.find("h2Reco") == std::string::npos && objName.find("h2MC") == std::string::npos)
                h2Map_old[objName] = (TH2D*)obj;
            else if(objName.find("h2Reco") != std::string::npos)recoh2Map[objName] = (TH2D*)obj;
            else if(objName.find("h2MC") != std::string::npos)mch2Map[objName] = (TH2D*)obj;
        }
    }

    auto inFile_old = TFile::Open("jets_edm4eic_R1.0_conPt0.2-100.0.hist.root", "READ");
    std::cout<<"Reading from file: "<<inFile_old->GetName()<<std::endl;
    TIter inFileKeys_old(inFile_old->GetListOfKeys());
    TKey* key_old;
    while((key_old = (TKey*)inFileKeys_old())){
        auto obj_old = key_old->ReadObj();
        std::string objName = obj_old->GetName();
        std::cout<<"_____Reading object: "<<objName<<std::endl;
        if(obj_old->IsA()->InheritsFrom(TH2D::Class())){
            if(objName.find("h2Reco") == std::string::npos && objName.find("h2MC") == std::string::npos)
                h2Map[objName] = (TH2D*)obj_old;
        }
    }

    std::unordered_map<std::string, TCanvas*> canvases;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    for(auto& h : h2Map){
        canvases[h.first] = new TCanvas(h.first.c_str(), h.first.c_str(), 1200, 600);
        canvases[h.first]->Divide(2, 1, 0, 0);
        canvases[h.first]->cd(1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.2);
        h.second->Draw("COL");
        canvases[h.first]->cd(2);
        gPad->SetLogz();
        gPad->SetRightMargin(0.2);
        TLatex *tex = new TLatex(0.4, 0.3, "With bijective mapping");
        tex->SetNDC();
        tex->SetTextSize(0.04);
        tex->SetTextColor(kRed);
        h2Map_old[h.first]->Draw("COLZ");    
        tex->Draw();
    }

    for(auto& h : mch1Map){
        canvases[h.first] = new TCanvas(h.first.c_str(), h.first.c_str(), 600, 600);
        std::string recoH1Name = h.first;
        recoH1Name.replace(h.first.find("MC"), 2, "Reco");
        gPad->SetLogy();
        h.second->Scale(1./h.second->Integral(), "width");
        h.second->SetLineColor(kRed+1);
        recoh1Map[recoH1Name]->Scale(1./recoh1Map[recoH1Name]->Integral(), "width");
        recoh1Map[recoH1Name]->SetLineColor(kBlue+1);
        h.second->Draw("EHIST");
        recoh1Map[recoH1Name]->Draw("EHIST SAME");
        TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
        leg->SetTextSize(0.03);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h.second, "Truth", "l");
        leg->AddEntry(recoh1Map[recoH1Name], "Reco", "l");
        leg->Draw();
    }

   // gStyle->SetOptTitle(1);
    for(auto& h : mch2Map){
        std::string recoH2Name = h.first;
        recoH2Name.replace(h.first.find("MC"), 2, "Reco");
        canvases[h.first] = new TCanvas(h.first.c_str(), h.first.c_str(), 1200, 600);
        canvases[h.first]->Divide(2, 1, 0, 0);
        canvases[h.first]->cd(1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.2);
        //gPad->SetTopMargin(0.08);
        h.second->Draw("COL");
        canvases[h.first]->cd(2);
        gPad->SetLogz();
        gPad->SetRightMargin(0.2);
        //gPad->SetTopMargin(0.08);
        //TLatex *tex = new TLatex(0.4, 0.3, "With bijective mapping");
        //tex->SetNDC();
        //tex->SetTextSize(0.04);
        //tex->SetTextColor(kRed);
        recoh2Map[recoH2Name]->Draw("COLZ");    
        //tex->Draw();
    }

    for(auto& c : canvases){
        c.second->Print((c.first+".pdf").c_str());
    }

}