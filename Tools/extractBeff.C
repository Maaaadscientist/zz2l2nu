void extractBeff(TString inputFile) {
  gStyle->SetOptStat(0);
  TFile *f= TFile::Open(inputFile);

  TH2D *btagEff_den_bjet = (TH2D *) f->Get("btagEff_den_bjet");
  TH2D *btagEff_num_bjet_tagged = (TH2D *) f->Get("btagEff_num_bjet");


  TH2D *btagEff_den_cjet = (TH2D *) f->Get("btagEff_den_cjet");
  TH2D *btagEff_num_cjet_tagged = (TH2D *) f->Get("btagEff_num_cjet");


  TH2D *btagEff_den_udsgjet = (TH2D *) f->Get("btagEff_den_udsgjet");
  TH2D *btagEff_num_udsgjet_tagged = (TH2D *) f->Get("btagEff_num_udsgjet");



  TH2D *b_all = (TH2D *)btagEff_den_bjet->Clone();
  TH2D *b_tagged = (TH2D *)btagEff_num_bjet_tagged->Clone();
  b_tagged->Divide(b_all);


  TH2D *c_all = (TH2D *)btagEff_den_cjet->Clone();
  TH2D *c_tagged = (TH2D *)btagEff_num_cjet_tagged->Clone();
  c_tagged->Divide(c_all);


  TH2D *udsg_all = (TH2D *)btagEff_den_udsgjet->Clone();
  TH2D *udsg_tagged = (TH2D *)btagEff_num_udsgjet_tagged->Clone();
  udsg_tagged->Divide(udsg_all);



  TCanvas *c1 = new TCanvas ("c1","c1",800,800);
  c1->cd();
  c1->SetLogx();
  b_tagged->SetTitle("btag efficiency (b jets)");
  b_tagged->Draw("colztexte");
  c1->Print("data/efficiencyTables/btag_efficiency_plots/btag-b.pdf");


  TCanvas *c2 = new TCanvas ("c2","c2",800,800);
  c2->cd();
  c2->SetLogx();
  c_tagged->SetTitle("btag efficiency (c jets)");
  c_tagged->Draw("colztexte");
  c2->Print("data/efficiencyTables/btag_efficiency_plots/btag-c.pdf");


  TCanvas *c3 = new TCanvas ("c3","c3",800,800);
  c3->cd();
  c3->SetLogx();
  udsg_tagged->SetTitle("btag efficiency (udsg jets)");
  udsg_tagged->Draw("colztexte");
  c3->Print("data/efficiencyTables/btag_efficiency_plots/btag-udsg.pdf");



  ofstream bjet;
  bjet.open("data/efficiencyTables/btag-b.txt");
  ofstream cjet;
  cjet.open("data/efficiencyTables/btag-c.txt");
  ofstream udsgjet;
  udsgjet.open("data/efficiencyTables/btag-udsg.txt");

  TH2D *h = (TH2D *)btagEff_den_udsgjet->Clone();

  for(int i=1; i <= h->GetNbinsY(); i++) {
    for(int j=1; j <= h->GetNbinsX(); j++) {
      udsgjet
		<< h->GetYaxis()->GetBinCenter(i) - h->GetYaxis()->GetBinWidth(i)/2
		<< "  "
		<< h->GetYaxis()->GetBinCenter(i) + h->GetYaxis()->GetBinWidth(i)/2
		<< "  ";
      udsgjet
		<< h->GetXaxis()->GetBinCenter(j) - h->GetXaxis()->GetBinWidth(j)/2
		<< "  "
		<< h->GetXaxis()->GetBinCenter(j) + h->GetXaxis()->GetBinWidth(j)/2
		<< "  ";
      udsgjet << udsg_tagged->GetBinContent(j,i) << " 0.  0. " << endl;


      bjet
		<< h->GetYaxis()->GetBinCenter(i) - h->GetYaxis()->GetBinWidth(i)/2
		<< "  "
		<< h->GetYaxis()->GetBinCenter(i) + h->GetYaxis()->GetBinWidth(i)/2
		<< "  ";
      bjet
		<< h->GetXaxis()->GetBinCenter(j) - h->GetXaxis()->GetBinWidth(j)/2
		<< "  "
		<< h->GetXaxis()->GetBinCenter(j) + h->GetXaxis()->GetBinWidth(j)/2
		<< "  ";
      bjet << b_tagged->GetBinContent(j,i) <<" 0.  0. " << endl;


      cjet
		<< h->GetYaxis()->GetBinCenter(i) - h->GetYaxis()->GetBinWidth(i)/2
		<< "  "
		<< h->GetYaxis()->GetBinCenter(i) + h->GetYaxis()->GetBinWidth(i)/2
		<< "  ";
      cjet
		<< h->GetXaxis()->GetBinCenter(j) - h->GetXaxis()->GetBinWidth(j)/2
		<< "  "
		<< h->GetXaxis()->GetBinCenter(j) + h->GetXaxis()->GetBinWidth(j)/2
		<< "  ";
      cjet << c_tagged->GetBinContent(j,i) << " 0.  0. " << endl;
    } 
  }

}
