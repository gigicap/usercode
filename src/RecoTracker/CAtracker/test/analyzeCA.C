void analyzeCA(TString filename1, TString filename2){

    int n123_1 , n234_1 , n345_1;
    int n123_2 , n234_2 , n345_2;
        
    
    gStyle->SetOptStat("");
    
    cout<<"Opening files"<<endl;
    
    TFile *f1 = new TFile(filename1,"READ");
    TFile *f2 = new TFile(filename2,"READ");

    cout<<"Getting tree"<<endl;
    
    TTree *t1 = f1->Get("CASeedingStep/evtree");
    TH1I *h123_1 = new TH1I("h123_1","Number of reconstructed triplets (Layer conf = 123)",100,0.,150);
    TH1I *h234_1 = new TH1I("h234_1","Number of reconstructed triplets (Layer conf = 234)",100,0.,200);
    TH1I *h345_1 = new TH1I("h345_1","Number of reconstructed triplets (Layer conf = 345)",100,0.,500);
    
    h123_1->SetLineColor(kRed);
    h234_1->SetLineColor(kRed);
    h345_1->SetLineColor(kRed);

    
    t1->SetBranchAddress("prodCells123", &n123_1);
    t1->SetBranchAddress("prodCells234", &n234_1);
    t1->SetBranchAddress("prodCells345", &n345_1);
    
    Int_t nentries = t1->GetEntries();
    
    for(int i = 0; i<nentries; i++){
        t1->GetEntry(i);
        
        h123_1->Fill(n123_1);
        h234_1->Fill(n234_1);
        h345_1->Fill(n345_1);

    }
    

    
    TTree *t2 = f2->Get("CASeedingStep/evtree");
    TH1I *h123_2 = new TH1I("h123_2","Number of reconstructed triplets (Layer conf = 123)",100,0.,150);
    TH1I *h234_2 = new TH1I("h234_2","Number of reconstructed triplets (Layer conf = 234)",100,0.,200);
    TH1I *h345_2 = new TH1I("h345_2","Number of reconstructed triplets (Layer conf = 345)",100,0.,500);
    
    h123_2->SetLineColor(kBlue);
    h234_2->SetLineColor(kBlue);
    h345_2->SetLineColor(kBlue);


    TLegend *l1 = new TLegend(0.4,0.5,0.7,0.7);
    l1->AddEntry(h123_1,"w/o ClusterShapeFilter","l");
    l1->AddEntry(h234_2,"w ClusterShapeFilter","l");
    
    

    t2->SetBranchAddress("prodCells123", &n123_2);
    t2->SetBranchAddress("prodCells234", &n234_2);
    t2->SetBranchAddress("prodCells345", &n345_2);
    
    Int_t nentries = t2->GetEntries();
    
    for(int i = 0; i<nentries; i++){
        t2->GetEntry(i);
        
        h123_2->Fill(n123_2);
        h234_2->Fill(n234_2);
        h345_2->Fill(n345_2);

    }
    


    new TCanvas;
    h123_1->Draw();
    h123_2->Draw("same");
    l1->Draw("same");

    new TCanvas;
    h234_1->Draw();
    h234_2->Draw("same");
    l1->Draw("same");

    new TCanvas;
    h345_1->Draw();
    h345_2->Draw("same");
    l1->Draw("same");    
    
    cout<<"setting"<<endl;
    
    return;
}

