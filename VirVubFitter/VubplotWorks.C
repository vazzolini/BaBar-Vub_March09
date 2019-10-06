void VubplotWorks(TString thepath, TString thefile, TString theoutpath, Int_t err=-1, Int_t staterr=0) {
  //modify to rescale Bauer for the bug in acceptance
  //err=1 plot BRBR with error bars (default)
  //err=2 plot absolute value or errors
  //err=3 plot relative errors
  //staterr=0 only data statistics (default)
  //staterr=1 data+MC statistics
  //  gROOT->SetStyle("Plain"); 
  gROOT->SetStyle("BABAR"); 
  char filetsC[100];
  char name[100], name2[100], nameeps[100], title[100];
  char buffertC[100];
  char canvtitl[100], leghead[100], nameY[100];
  double pbr, err1;
  double vub1, vub2, errvub1, errvub2;
  int bin;

  //take the partial branching fraction, and total relative error
  //correct for wrong pstarfactor
  //compute vub according to 0.00432 sqrt(PBF/0.002/epsphsp*1.55/taub)

  double tauB = 1.61;
  double fact1 =0.00424;
  double fact2 =0.00445;
  double denom = 0.002;
  double fixT = 1.61;
  double unc = 0.048;
  double brxsl = 0.1083;
  //  double csi = 128.3;
  double csi = 0.00779;
  double scfct=pow(4.74/4.7,5.);
  //  double scfct=1;

  //  double G15[16]  = {0,0,0,0,0,0,0,0.15,0.185,0.185,0.17,0.15,0.12,0.10,0.07,0.05};
  double G15[16]  = {0,0,0,0,0,0,0,0,0,0.185,0.17,0.15,0.12,0.10,0.07,0.05};
  double G17[16]  = {0,0,0,0,0,0.34,0.325,0.3,0.27,0.23,0.195,0.16,0.13,0.10,0.07,0.05};
  double G186[16] = {0,0,0,0,0,0.415,0.375,0.33,0.29,0.245,0.205,0.165,0.13,0.10,0.07,0.05};
  double R15[16]  = {0,0,0,0,0,0,0,0,0,-0.135,-0.1025,-0.0725,-0.045,-0.0225,-0.004,0};
  double R17[16]  = {0,0,0,0,0,-0.1125,-0.1,-0.078,-0.06,-0.04,-0.027,-0.015,-0.005,0,0,0};
  double R186[16] = {0,0,0,0,0,-0.055,-0.042,-0.03,-0.02,-0.0125,-0.006,-0.002,0,0,0,0};

  double eG15[8] = {0,0,0,0,0,0.241661,0.331059,0.536027};
  double eG17[8] = {0,0,0,0.17,0.181797,0.232594,0.328261,0.536027};
  double eG186[8]= {0,0,0,0.146287,0.179234,0.232594,0.333654,0.536027};


  //new numbers from Bauer

  double bG17[8]  = {0., 0., 0.379, 0.343, 0.281, 0.209, 0.140, 0.08};
  //shape fcn variation (according to what?)
  double dbG17[8] = {0., 0., 0.341, 0.312, 0.259, 0.197, 0.137, 0.08};
  double ubG17[8] = {0., 0., 0.408, 0.365, 0.295, 0.215, 0.141, 0.08};
  //relative errors
  double roG17[8] = {0., 0., 0.25, 0.18, 0.15, 0.16, 0.18, 0.21};
  double rmG17[8] = {0., 0., 0.06, 0.07, 0.08, 0.11, 0.16, 0.28}; 
 
  double rupG17[8], rdnG17[8]; 
  double aupG17[8], adnG17[8];

  for(int i=0; i<8; i++) {
    int indx=i*2;
    if(bG17[i] != 0){
      rupG17[i] = sqrt(pow((ubG17[i]-bG17[i])/bG17[i],2)+pow(roG17[i],2)+pow(rmG17[i],2));
      rdnG17[i] = sqrt(pow((dbG17[i]-bG17[i])/bG17[i],2)+pow(roG17[i],2)+pow(rmG17[i],2));
      aupG17[i] = rupG17[i]*bG17[i];
      adnG17[i] = rdnG17[i]*bG17[i];
      //print out to compare any differences
      cout << "BLL original:     " << indx << "  " << G17[indx]*scfct  << " +" << 
	eG17[i]*scfct*G17[indx] << " -" << eG17[i]*scfct*G17[indx] << endl;
      cout << "BLL from Bauer :  " << indx << "  " << bG17[i]  << " +" << aupG17[i] << " -" << adnG17[i] << endl;
      //CB CHANGE HERE TO THE NEW NUMBERS PROVIDED BY BAUER
      G17[indx]=bG17[i];
      eG17[i]=0.5*(rupG17[i]+rdnG17[i]);
      cout << "now using for Vub : " << indx << " " << G17[indx] << " +/- " << eG17[i] << endl;
      scfct = 1.;
    }
  }

  //CLEO points

  double eff15[16]  = {0.537675, 0.515206, 0.491637, 0.464306, 0.434074, 0.403882, 0.370748, 0.336415, 0.30042,
		       0.264906, 0.228211, 0.193597, 0.159044, 0.127571, 0.0984194,0.0721088};
  double eff17[16]  = {0.67617,  0.641757, 0.605362, 0.565946, 0.522629, 0.477851, 0.430992, 0.383713, 0.336975, 
		       0.290836, 0.244558, 0.202481, 0.162845, 0.128451, 0.0984194,0.0721088};
  double eff186[16] = {0.767567, 0.723249, 0.67601,  0.626311, 0.572849, 0.518307, 0.462845, 0.407063, 0.352841, 
		       0.30068,  0.24938,  0.204162, 0.162985, 0.128451, 0.0984194, 0.0721088};

  double edfn15p[8] = {0.209452, 0.201863, 0.197533, 0.181355, 0.16153, 0.147316, 0.132007, 0.16177};
  double edfn15m[8] = {-0.24064,-0.242683,-0.245107,-0.242252,-0.241771,-0.233951,-0.22905,-0.229124};

  double edfn17p[8] = {0.165436, 0.154011, 0.140397, 0.128995, 0.114536, 0.112127, 0.116275, 0.161973};
  double edfn17m[8] = {-0.235634,-0.232499,-0.231255,-0.223188,-0.218839,-0.208354,-0.212286,-0.228921};

  double edfn186p[8] = {0.130955, 0.117357, 0.105347, 0.0951093, 0.0890969, 0.0984896, 0.11593, 0.161973};
  double edfn186m[8] = {-0.220457,-0.215095,-0.208946,-0.200025,-0.195122,-0.188209,-0.209895,-0.228921};


  //BELLE points

//   double eff15[16] = {0.422911, 0.408634, 0.39322, 0.374437, 0.352961, 0.330947, 0.307178, 0.280997, 
// 		       0.253121, 0.224427, 0.195414, 0.164526, 0.135394, 0.107378, 0.0801595, 0.0585843};

//   double eff17[16] = {0.570887, 0.545862, 0.518505, 0.487717, 0.452223, 0.41661, 0.38002, 0.340199, 
// 		       0.299721, 0.258664, 0.217807, 0.177248, 0.140857, 0.108694, 0.0802193, 0.0585843}; 

//   double eff186[16] = {0.683729, 0.647637, 0.608495, 0.56656, 0.518365, 0.471107, 0.42319, 0.372522, 
// 			0.322812, 0.273619, 0.225424, 0.179821, 0.141196, 0.108694, 0.0802193, 0.0585843}; 

//   double edfn15p[8] = {0.0980968, 0.0984567, 0.0986731, 0.097414, 0.0919901, 0.0846722, 0.0755308, 0.0840579};
//   double edfn15m[8] = {-0.110725, -0.107117, -0.102505, -0.107017, -0.106972, -0.106458, -0.0973049, -0.0865729};

//   double edfn17p[8] = {0.0969043, 0.0936212, 0.0901938, 0.0827782, 0.0751565, 0.0698312, 0.0645315, 0.0834982};
//   double edfn17m[8] = {-0.100249, -0.100004, -0.0985491, -0.0961657, -0.0929852, -0.0874882, -0.081612, -0.0862592};

//   double edfn186p[8] = {0.0809665, 0.076463, 0.0711819, 0.0620341, 0.055078, 0.0571218, 0.0633879, 0.0834982};
//   double edfn186m[8] = {-0.099985, -0.0978889, -0.0932955, -0.089035, -0.082958, -0.0796461, -0.0802844, -0.0862592};


  //old points
  //  double eff15[16]  = {0.603798, 0.575883, 0.545139, 0.513021, 0.477495, 0.440415, 0.403535, 0.364562, 0.32527, 
  // 		       0.285202, 0.244994, 0.206619, 0.1693, 0.135966, 0.105959, 0.0787224};
  
  //   double eff17[16]  = {0.736655, 0.695391, 0.651477, 0.607005, 0.557154, 0.506764, 0.457092, 0.406523, 
  // 		       0.355915, 0.307039, 0.259041, 0.213872, 0.172189, 0.136424, 0.105959, 0.0787224};
  
  //   double eff186[16] = {0.815139, 0.764909, 0.712307, 0.657714, 0.598617, 0.540059, 0.482456, 0.424615, 
  // 		       0.36777,  0.313415, 0.262368, 0.214987, 0.172269, 0.136424, 0.105959, 0.0787224};
  
  TCanvas  *c1 = new TCanvas("canv","canv",800,600);
  c1->cd(1);
  
  for(int ij=staterr; ij<=staterr; ij++){
    int lines = 0;
    if(err<0){
      sprintf(name,"%s%i","PBR",ij);  sprintf(title, "PBR"); TH1D* h = new TH1D(name, title, 8, -1,15);
      sprintf(name,"%s%i","PBR2",ij);  sprintf(title, "PBR2"); TH1D* h2 = new TH1D(name, title, 8, -0.9,15.1); 
    }else{
      sprintf(name,"%s%i","PBR",ij);  sprintf(title, "PBR"); TH1D* h = new TH1D(name, title, 15, 0.5+0.1*ij,15.5+0.1*ij);
      sprintf(name,"%s%i","PBR2",ij);  sprintf(title, "PBR2"); TH1D* h2 = new TH1D(name, title, 15, 0.5+0.1*ij,15.5+0.1*ij); 
    }
    if(ij==1){
    sprintf(filetsC,"%s%s%i",thepath.Data(),thefile.Data(),15);
    sprintf(leghead,"Mx < 1.5");
    }else if(ij==2){
    sprintf(filetsC,"%s%s%i",thepath.Data(),thefile.Data(),17);
    sprintf(leghead,"Mx < 1.7");
    }else if(ij==3){
    sprintf(filetsC,"%s%s%i",thepath.Data(),thefile.Data(),186);
    sprintf(leghead,"Mx < 1.86");
    }
    cout << "Open file  " << filetsC << ij << endl;
    ifstream tsC(filetsC);
    sprintf(name,"%s%i","PBR",ij);
    sprintf(name2,"%s%i","PBR2",ij);
    while (tsC.getline(buffertC, 200, '\n')) {
      sscanf(buffertC, "%lf %lf",&pbr,&err1);
      pbr/=1000.;
      //	cout << "control" <<  endl;
      int indx=2*lines;
      if(ij==1){
	//defazio neubert
	vub1=fact1*sqrt(pbr/denom/eff15[indx]*fixT/tauB)*1000.;
	errvub1=sqrt(pow(err1/2,2)+pow(unc,2)+pow(0.5*(edfn15p[lines]-edfn15m[lines]),2))*vub1;
	if(G15[indx]!=0){
	  vub2=sqrt(pbr*csi/G15[indx]/scfct)*1000;	    
	  errvub2=sqrt(pow(err1,2)+pow(eG15[lines],2))*vub2/2;
	}else{
	  vub2 = 0;
	  errvub2 = 0;
	}
      }else if(ij==2){
	//defazio neubert
	vub1=fact1*sqrt(pbr/denom/eff17[indx]*fixT/tauB)*1000;
	errvub1=sqrt(pow(err1/2,2)+pow(unc,2)+pow(0.5*(edfn17p[lines]-edfn17m[lines]),2))*vub1;
	if(G17[indx]!=0){
	  vub2=sqrt(pbr*csi/G17[indx]/scfct)*1000;	    
	  errvub2=sqrt(pow(err1,2)+pow(eG17[lines],2))*vub2/2;
	}else{
	  vub2 = 0;
	  errvub2 = 0;
	}
      }else if(ij==3){
	//defazio neubert
	vub1=fact1*sqrt(pbr/denom/eff186[indx]*fixT/tauB)*1000;
	errvub1=sqrt(pow(err1/2,2)+pow(unc,2)+pow(0.5*(edfn186p[lines]-edfn186m[lines]),2))*vub1;
	if(G17[indx]!=0){
	  vub2=sqrt(pbr*csi/G186[indx]/scfct)*1000;	    
	  errvub2=sqrt(pow(err1,2)+pow(eG186[lines],2))*vub2/2;
	}else{
	  vub2 = 0;
	  errvub2 = 0;
	}
      }
      
      lines++;
      cout << "PBR=          " << pbr << "+/-" << pbr*err1 << endl;
      cout << "|V_ub| DFN =   " << vub1 << "+/-" << errvub1 << endl; 
      cout << "|V_ub| BLL =   " << vub2 << "+/-" << errvub2 << endl; 
      sprintf(nameY,"|V_{ub}| (10^{-3})");
      ((TH1D*)gDirectory->Get(name))->SetBinContent(lines,vub1);
      ((TH1D*)gDirectory->Get(name2))->SetBinContent(lines,vub2);
      ((TH1D*)gDirectory->Get(name))->SetBinError(lines,errvub1);
      ((TH1D*)gDirectory->Get(name2))->SetBinError(lines,errvub2);
      //	cout << "PBR " << pbr << endl;
    }
    ((TH1D*)gDirectory->Get(name))->SetMinimum(0.);
    ((TH1D*)gDirectory->Get(name2))->SetMinimum(0.);
    //    ((TH1D*)gDirectory->Get(name))->SetTitleOffset(1.3,"Y");
    //    ((TH1D*)gDirectory->Get(name))->SetTitleOffset(1.3,"X");
    //    ((TH1D*)gDirectory->Get(name2))->SetTitleOffset(1.3,"Y");
    //    ((TH1D*)gDirectory->Get(name2))->SetTitleOffset(1.3,"X");
    ((TH1D*)gDirectory->Get(name))->SetXTitle("q^{2}_{cut} (GeV^{2}/c^{4})");
    ((TH1D*)gDirectory->Get(name2))->SetXTitle("q^{2}_{cut} (GeV^{2}/c^{4})");
    ((TH1D*)gDirectory->Get(name))->SetYTitle(nameY);
    ((TH1D*)gDirectory->Get(name2))->SetYTitle(nameY);
    ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1);  ((TH1D*)gDirectory->Get(name))->SetLineWidth(2);
    ((TH1D*)gDirectory->Get(name2))->SetMarkerSize(1.2);  ((TH1D*)gDirectory->Get(name2))->SetLineWidth(2);
    if(TMath::Abs(err)==2){
      ((TH1D*)gDirectory->Get(name))->SetMaximum(8);
      ((TH1D*)gDirectory->Get(name2))->SetMaximum(8);
    }else if(TMath::Abs(err)==1){
      ((TH1D*)gDirectory->Get(name))->SetMaximum(8);
      ((TH1D*)gDirectory->Get(name2))->SetMaximum(8);
    }else if(TMath::Abs(err)==3){
      ((TH1D*)gDirectory->Get(name))->SetMaximum(0.26);
      ((TH1D*)gDirectory->Get(name2))->SetMaximum(0.26);
    }
    
    //plot filled histograms
    gStyle->SetOptStat(0);  gStyle->SetOptTitle(0);
    
    TLegendEntry *legge; 
    TLegend *leg;
    leg = new TLegend(0.7,0.2,0.9,0.45);
    //  leg = new TLegend(0.65,0.15,0.85,0.30);
    leg->SetBorderSize(0.); leg->SetTextSize(0.07);
    //make text size 0.08 to improve reading 
    leg->SetFillStyle(4000);  leg->SetFillColor(0); 
    
    ((TH1D*)gDirectory->Get(name2))->SetMarkerStyle(20); ((TH1D*)gDirectory->Get(name2))->SetMarkerColor(1);
    ((TH1D*)gDirectory->Get(name2))->SetLineColor(1);
    ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(24); ((TH1D*)gDirectory->Get(name))->SetMarkerColor(2);
    ((TH1D*)gDirectory->Get(name))->SetLineColor(2);

    /*
    ((TH1D*)gDirectory->Get("PBR2"))->SetMarkerStyle(21); ((TH1D*)gDirectory->Get("PBR2"))->SetMarkerColor(2);
    ((TH1D*)gDirectory->Get("PBR22"))->SetMarkerStyle(25); ((TH1D*)gDirectory->Get("PBR22"))->SetMarkerColor(2);
    ((TH1D*)gDirectory->Get("PBR22"))->SetLineColor(2);
    ((TH1D*)gDirectory->Get("PBR2"))->SetLineColor(2);
    
    ((TH1D*)gDirectory->Get("PBR1"))->SetMarkerStyle(20); ((TH1D*)gDirectory->Get("PBR1"))->SetMarkerColor(1);
    ((TH1D*)gDirectory->Get("PBR21"))->SetMarkerStyle(24); ((TH1D*)gDirectory->Get("PBR21"))->SetMarkerColor(1);
    ((TH1D*)gDirectory->Get("PBR21"))->SetLineColor(1);
    ((TH1D*)gDirectory->Get("PBR1"))->SetLineColor(1);
    */

    ((TH1D*)gDirectory->Get(name2))->Draw("p");
    ((TH1D*)gDirectory->Get(name))->Draw("psame");

    legge = leg->AddEntry(((TH1D*)gDirectory->Get(name2)),  "BLL", "p"); 
    legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)),  "DFN", "p"); 
    sprintf(leghead,"(b)");
    legge = leg->SetHeader(leghead);
    leg->Draw();
    
    gPad->Update();
    sprintf(nameeps,"%s%s%s%i%i%s",theoutpath.Data(),thefile.Data(),"_s_",TMath::Abs(err),staterr,".eps");
    c1->Print(nameeps);
    sprintf(nameeps,"%s%s%s%i%i%s",theoutpath.Data(),thefile.Data(),"_s_",TMath::Abs(err),staterr,".gif");
    c1->Print(nameeps);
  }
}
