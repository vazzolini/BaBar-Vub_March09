
void redAll() {

  redFitTree("/nfs/farm/babar/AWG6/Recoil/newprod/root/anaQA-prl/reduced/csx-vubnre_small.root", 
                                          "/u/ec/ursl/d/root/b2ufit/prl/rcsx-vubnre_small.root");

  redFitTree("/nfs/farm/babar/AWG6/Recoil/newprod/root/anaQA-prl/reduced/csx-vubmix_small.root", 
                                          "/u/ec/ursl/d/root/b2ufit/prl/rcsx-vubmix_small.root");

  redFitTree("/nfs/farm/babar/AWG6/Recoil/newprod/root/anaQA-prl/reduced/csx-genb-new_small.root", 
                                          "/u/ec/ursl/d/root/b2ufit/prl/rcsx-genb-new_small.root");

  redFitTree("/nfs/farm/babar/AWG6/Recoil/newprod/root/anaQA-prl/reduced/csx-data_small.root", 
                                          "/u/ec/ursl/d/root/b2ufit/prl/rcsx-data_small.root");

}


void redFitTree(const char *filename, const char *newname) {

  cout << "Reducing " << filename << endl;

   // -- Get old file, old tree and set top branch address
   TFile *oldfile = new TFile(filename);
   TTree *oldtree = (TTree*)oldfile->Get("events");

   // -- Disable all branches ...
   oldtree->SetBranchStatus("*",0);

   // -- ... and switch on those you'd like to write out into the new tree
   oldtree->SetBranchStatus("mes",1);
   oldtree->SetBranchStatus("intpur",1);
   oldtree->SetBranchStatus("brecoflav", 1); 

   oldtree->SetBranchStatus("mxhadfit",1);
   oldtree->SetBranchStatus("q2fit",1);
   oldtree->SetBranchStatus("brecocharge",1);
   oldtree->SetBranchStatus("xcharge",1);
   oldtree->SetBranchStatus("lcharge",1);
   oldtree->SetBranchStatus("mm2",1);
   oldtree->SetBranchStatus("nle",1);
   oldtree->SetBranchStatus("nel",1);
   oldtree->SetBranchStatus("nmu",1);
   oldtree->SetBranchStatus("pcms",1);
   oldtree->SetBranchStatus("tlab",1);
   oldtree->SetBranchStatus("plab",1);
   oldtree->SetBranchStatus("nchg",1);
   oldtree->SetBranchStatus("nneu",1);
   oldtree->SetBranchStatus("nks",1);
   oldtree->SetBranchStatus("nkp",1);
   oldtree->SetBranchStatus("wdeltam",1);

//    oldtree->SetBranchStatus("totweightTrkMult",1);
//    oldtree->SetBranchStatus("totweightNutMult",1);
//    oldtree->SetBranchStatus("totweight",1);


   oldtree->SetBranchStatus("vub",1);
   oldtree->SetBranchStatus("vcb",1);
   oldtree->SetBranchStatus("other",1);
   oldtree->SetBranchStatus("kplus",1);


   oldtree->SetBranchStatus("ecmsgen",1);
   oldtree->SetBranchStatus("mxhadgen",1);
   oldtree->SetBranchStatus("q2Gen",1);
   oldtree->SetBranchStatus("Gvxbtyp",1);
   oldtree->SetBranchStatus("GfDpi",1);
   oldtree->SetBranchStatus("GfDk",1);
   oldtree->SetBranchStatus("GfDks",1);
   oldtree->SetBranchStatus("GfDpiz",1);
   oldtree->SetBranchStatus("GfDlep",1);
			     
   //-- Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile(newname,"recreate");
   TTree *newtree = oldtree->CloneTree();

   newtree->Print();
   newfile->Write();
   delete oldfile;
   delete newfile;
}
