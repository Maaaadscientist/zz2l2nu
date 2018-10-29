void extractBeff(){
  gStyle->SetOptStat(0);
  TFile *f= TFile::Open("/user/npostiau/CMSSW_8_0_26_patch1_shears_uncertainties/src/shears/HZZ2l2nu/OUTPUTS/testWithNewBonzais/MERGED/outputHZZ_ZZTo2L2Nu.root"); //Change the path here if needed.
  TH2D *btagEff_den_bjet = (TH2D *) f->Get("btagEff_den_bjet");
  TH2D *btagEff_num_bjet_tagged_loose = (TH2D *) f->Get("btagEff_num_bjet_tagged_loose");
  TH2D *btagEff_num_bjet_tagged_medium = (TH2D *) f->Get("btagEff_num_bjet_tagged_medium");
  TH2D *btagEff_num_bjet_tagged_tight = (TH2D *) f->Get("btagEff_num_bjet_tagged_tight");


  TH2D *btagEff_den_cjet = (TH2D *) f->Get("btagEff_den_cjet");
  TH2D *btagEff_num_cjet_tagged_loose = (TH2D *) f->Get("btagEff_num_cjet_tagged_loose");
  TH2D *btagEff_num_cjet_tagged_medium = (TH2D *) f->Get("btagEff_num_cjet_tagged_medium");
  TH2D *btagEff_num_cjet_tagged_tight = (TH2D *) f->Get("btagEff_num_cjet_tagged_tight");


  TH2D *btagEff_den_udsgjet = (TH2D *) f->Get("btagEff_den_udsgjet");
  TH2D *btagEff_num_udsgjet_tagged_loose = (TH2D *) f->Get("btagEff_num_udsgjet_tagged_loose");
  TH2D *btagEff_num_udsgjet_tagged_medium = (TH2D *) f->Get("btagEff_num_udsgjet_tagged_medium");
  TH2D *btagEff_num_udsgjet_tagged_tight = (TH2D *) f->Get("btagEff_num_udsgjet_tagged_tight");




  TH2D *b_all = (TH2D *)btagEff_den_bjet->Clone();
  TH2D *b_tagged_loose = (TH2D *)btagEff_num_bjet_tagged_loose->Clone();
  TH2D *b_tagged_medium = (TH2D *)btagEff_num_bjet_tagged_medium->Clone();
  TH2D *b_tagged_tight = (TH2D *)btagEff_num_bjet_tagged_tight->Clone();

  b_tagged_loose->Divide(b_all);
  b_tagged_medium->Divide(b_all);
  b_tagged_tight->Divide(b_all);

  TH2D *c_all = (TH2D *)btagEff_den_cjet->Clone();
  TH2D *c_tagged_loose = (TH2D *)btagEff_num_cjet_tagged_loose->Clone();
  TH2D *c_tagged_medium = (TH2D *)btagEff_num_cjet_tagged_medium->Clone();
  TH2D *c_tagged_tight = (TH2D *)btagEff_num_cjet_tagged_tight->Clone();

  c_tagged_loose->Divide(c_all);
  c_tagged_medium->Divide(c_all);
  c_tagged_tight->Divide(c_all);

  TH2D *udsg_all = (TH2D *)btagEff_den_udsgjet->Clone();
  TH2D *udsg_tagged_loose = (TH2D *)btagEff_num_udsgjet_tagged_loose->Clone();
  TH2D *udsg_tagged_medium = (TH2D *)btagEff_num_udsgjet_tagged_medium->Clone();
  TH2D *udsg_tagged_tight = (TH2D *)btagEff_num_udsgjet_tagged_tight->Clone();

  udsg_tagged_loose->Divide(udsg_all);
  udsg_tagged_medium->Divide(udsg_all);
  udsg_tagged_tight->Divide(udsg_all);




  TCanvas *c1 = new TCanvas ("c1","c1",800,800);
  c1->cd();
  c1->SetLogx();
  b_tagged_loose->SetTitle("btag (Loose) efficiency (b jets)");
  b_tagged_loose->Draw("colztexte");
  c1->Print("data/efficiencyTables/btag_efficiency_plots/btag-loose-b.pdf");

  TCanvas *c11 = new TCanvas ("c11","c11",800,800);
  c11->cd();
  c11->SetLogx();
  b_tagged_medium->SetTitle("btag (Medium) efficiency (b jets)");
  b_tagged_medium->Draw("colztexte");
  c11->Print("data/efficiencyTables/btag_efficiency_plots/btag-medium-b.pdf");

  TCanvas *c111 = new TCanvas ("c111","c111",800,800);
  c111->cd();
  c111->SetLogx();
  b_tagged_tight->SetTitle("btag (Tight) efficiency (b jets)");
  b_tagged_tight->Draw("colztexte");
  c111->Print("data/efficiencyTables/btag_efficiency_plots/btag-tight-b.pdf");



  TCanvas *c2 = new TCanvas ("c2","c2",800,800);
  c2->cd();
  c2->SetLogx();
  c_tagged_loose->SetTitle("btag (Loose) efficiency (c jets)");
  c_tagged_loose->Draw("colztexte");
  c2->Print("data/efficiencyTables/btag_efficiency_plots/btag-loose-c.pdf");

  TCanvas *c22 = new TCanvas ("c22","c22",800,800);
  c22->cd();
  c22->SetLogx();
  c_tagged_medium->SetTitle("btag (Medium) efficiency (c jets)");
  c_tagged_medium->Draw("colztexte");
  c22->Print("data/efficiencyTables/btag_efficiency_plots/btag-medium-c.pdf");

  TCanvas *c222 = new TCanvas ("c222","c222",800,800);
  c222->cd();
  c222->SetLogx();
  c_tagged_tight->SetTitle("btag (Tight) efficiency (c jets)");
  c_tagged_tight->Draw("colztexte");
  c222->Print("data/efficiencyTables/btag_efficiency_plots/btag-tight-c.pdf");


  TCanvas *c3 = new TCanvas ("c3","c3",800,800);
  c3->cd();
  c3->SetLogx();
  udsg_tagged_loose->SetTitle("btag (Loose) efficiency (udsg jets)");
  udsg_tagged_loose->Draw("colztexte");
  c3->Print("data/efficiencyTables/btag_efficiency_plots/btag-loose-udsg.pdf");

  TCanvas *c33 = new TCanvas ("c33","c33",800,800);
  c33->cd();
  c33->SetLogx();
  udsg_tagged_medium->SetTitle("btag (Medium) efficiency (udsg jets)");
  udsg_tagged_medium->Draw("colztexte");
  c33->Print("data/efficiencyTables/btag_efficiency_plots/btag-medium-udsg.pdf");

  TCanvas *c333 = new TCanvas ("c333","c333",800,800);
  c333->cd();
  c333->SetLogx();
  udsg_tagged_tight->SetTitle("btag (Tight) efficiency (udsg jets)");
  udsg_tagged_tight->Draw("colztexte");
  c333->Print("data/efficiencyTables/btag_efficiency_plots/btag-tight-udsg.pdf");




  ofstream bjetl;
  bjetl.open("data/efficiencyTables/btag-loose-b.txt");
  ofstream cjetl;
  cjetl.open("data/efficiencyTables/btag-loose-c.txt");
  ofstream udsgjetl;
  udsgjetl.open("data/efficiencyTables/btag-loose-udsg.txt");

  ofstream bjetm;
  bjetm.open("data/efficiencyTables/btag-medium-b.txt");
  ofstream cjetm;
  cjetm.open("data/efficiencyTables/btag-medium-c.txt");
  ofstream udsgjetm;
  udsgjetm.open("data/efficiencyTables/btag-medium-udsg.txt");

  ofstream bjett;
  bjett.open("data/efficiencyTables/btag-tight-b.txt");
  ofstream cjett;
  cjett.open("data/efficiencyTables/btag-tight-c.txt");
  ofstream udsgjett;
  udsgjett.open("data/efficiencyTables/btag-tight-udsg.txt");

  TH2D *ccc = (TH2D *)btagEff_den_udsgjet->Clone();


  for(int i=1; i<=ccc->GetNbinsY();i++){
    for(int j=1; j<=ccc->GetNbinsX();j++){
      udsgjett<<ccc->GetYaxis()->GetBinCenter(i)-ccc->GetYaxis()->GetBinWidth(i)/2 << "  "<< ccc->GetYaxis()->GetBinCenter(i)+ccc->GetYaxis()->GetBinWidth(i)/2<<"  ";
      udsgjett<<ccc->GetXaxis()->GetBinCenter(j)-ccc->GetXaxis()->GetBinWidth(j)/2 << "  "<< ccc->GetXaxis()->GetBinCenter(j)+ccc->GetXaxis()->GetBinWidth(j)/2<<"  ";
      udsgjett<<udsg_tagged_tight->GetBinContent(j,i)<<" 0.  0. "<<endl;

      udsgjetm<<ccc->GetYaxis()->GetBinCenter(i)-ccc->GetYaxis()->GetBinWidth(i)/2 << "  "<< ccc->GetYaxis()->GetBinCenter(i)+ccc->GetYaxis()->GetBinWidth(i)/2<<"  ";
      udsgjetm<<ccc->GetXaxis()->GetBinCenter(j)-ccc->GetXaxis()->GetBinWidth(j)/2 << "  "<< ccc->GetXaxis()->GetBinCenter(j)+ccc->GetXaxis()->GetBinWidth(j)/2<<"  ";
      udsgjetm<<udsg_tagged_medium->GetBinContent(j,i)<<" 0.  0. "<<endl;

      udsgjetl<<ccc->GetYaxis()->GetBinCenter(i)-ccc->GetYaxis()->GetBinWidth(i)/2 << "  "<< ccc->GetYaxis()->GetBinCenter(i)+ccc->GetYaxis()->GetBinWidth(i)/2<<"  ";
      udsgjetl<<ccc->GetXaxis()->GetBinCenter(j)-ccc->GetXaxis()->GetBinWidth(j)/2 << "  "<< ccc->GetXaxis()->GetBinCenter(j)+ccc->GetXaxis()->GetBinWidth(j)/2<<"  ";
      udsgjetl<<udsg_tagged_loose->GetBinContent(j,i)<<" 0.  0. "<<endl;


      bjett<<ccc->GetYaxis()->GetBinCenter(i)-ccc->GetYaxis()->GetBinWidth(i)/2 << "  "<< ccc->GetYaxis()->GetBinCenter(i)+ccc->GetYaxis()->GetBinWidth(i)/2<<"  ";
      bjett<<ccc->GetXaxis()->GetBinCenter(j)-ccc->GetXaxis()->GetBinWidth(j)/2 << "  "<< ccc->GetXaxis()->GetBinCenter(j)+ccc->GetXaxis()->GetBinWidth(j)/2<<"  ";
      bjett<<b_tagged_tight->GetBinContent(j,i)<<" 0.  0. "<<endl;

      bjetm<<ccc->GetYaxis()->GetBinCenter(i)-ccc->GetYaxis()->GetBinWidth(i)/2 << "  "<< ccc->GetYaxis()->GetBinCenter(i)+ccc->GetYaxis()->GetBinWidth(i)/2<<"  ";
      bjetm<<ccc->GetXaxis()->GetBinCenter(j)-ccc->GetXaxis()->GetBinWidth(j)/2 << "  "<< ccc->GetXaxis()->GetBinCenter(j)+ccc->GetXaxis()->GetBinWidth(j)/2<<"  ";
      bjetm<<b_tagged_medium->GetBinContent(j,i)<<" 0.  0. "<<endl;

      bjetl<<ccc->GetYaxis()->GetBinCenter(i)-ccc->GetYaxis()->GetBinWidth(i)/2 << "  "<< ccc->GetYaxis()->GetBinCenter(i)+ccc->GetYaxis()->GetBinWidth(i)/2<<"  ";
      bjetl<<ccc->GetXaxis()->GetBinCenter(j)-ccc->GetXaxis()->GetBinWidth(j)/2 << "  "<< ccc->GetXaxis()->GetBinCenter(j)+ccc->GetXaxis()->GetBinWidth(j)/2<<"  ";
      bjetl<<b_tagged_loose->GetBinContent(j,i)<<" 0.  0. "<<endl;

      cjett<<ccc->GetYaxis()->GetBinCenter(i)-ccc->GetYaxis()->GetBinWidth(i)/2 << "  "<< ccc->GetYaxis()->GetBinCenter(i)+ccc->GetYaxis()->GetBinWidth(i)/2<<"  ";
      cjett<<ccc->GetXaxis()->GetBinCenter(j)-ccc->GetXaxis()->GetBinWidth(j)/2 << "  "<< ccc->GetXaxis()->GetBinCenter(j)+ccc->GetXaxis()->GetBinWidth(j)/2<<"  ";
      cjett<<c_tagged_tight->GetBinContent(j,i)<<" 0.  0. "<<endl;

      cjetm<<ccc->GetYaxis()->GetBinCenter(i)-ccc->GetYaxis()->GetBinWidth(i)/2 << "  "<< ccc->GetYaxis()->GetBinCenter(i)+ccc->GetYaxis()->GetBinWidth(i)/2<<"  ";
      cjetm<<ccc->GetXaxis()->GetBinCenter(j)-ccc->GetXaxis()->GetBinWidth(j)/2 << "  "<< ccc->GetXaxis()->GetBinCenter(j)+ccc->GetXaxis()->GetBinWidth(j)/2<<"  ";
      cjetm<<c_tagged_medium->GetBinContent(j,i)<<" 0.  0. "<<endl;

      cjetl<<ccc->GetYaxis()->GetBinCenter(i)-ccc->GetYaxis()->GetBinWidth(i)/2 << "  "<< ccc->GetYaxis()->GetBinCenter(i)+ccc->GetYaxis()->GetBinWidth(i)/2<<"  ";
      cjetl<<ccc->GetXaxis()->GetBinCenter(j)-ccc->GetXaxis()->GetBinWidth(j)/2 << "  "<< ccc->GetXaxis()->GetBinCenter(j)+ccc->GetXaxis()->GetBinWidth(j)/2<<"  ";
      cjetl<<c_tagged_loose->GetBinContent(j,i)<<" 0.  0. "<<endl;



    } 
  }

}
