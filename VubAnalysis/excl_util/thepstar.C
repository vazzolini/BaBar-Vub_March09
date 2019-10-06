void thepstar(int leptype, int mode){

  TChain *chainvcb = new TChain("events");
  TChain *chainvub = new TChain("events");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-b0cock-2000.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-b0cock-2001.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-b0cock-2002.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-bpcock-2000.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-bpcock-2001.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-bpcock-2002.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-bpcock-2002.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2000a.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2000b.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2000c.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2001a.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2001b.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2001c.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2001d.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2001e.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2001f.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2002a.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2002b.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2002c.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2002d.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbch-new-2002e.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2000a.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2000b.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2000c.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2001a.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2001b.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2001c.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2001d.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2001e.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2002a.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2002b.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2002c.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2002d.root");
  chainvcb->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-genbnu-new-2002e.root");

  chainvub->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-reso-new.root");
  chainvub->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-brevubmix-old2001.root");
  chainvub->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-eta.root");
  chainvub->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-etap.root");
  chainvub->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-a0.root");
  chainvub->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-vubmix-new2000.root");
  chainvub->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-vubmix-new2001.root");
  chainvub->Add("/nfs/babar/recoil/Vub_excl/anaQA-prod_042204/root/anaQA-excl00/csx-vubmix-new2002.root");

  pstarfactor  pippo;
  pippo.Bookhist(mode);
  pippo.Init(chainvcb);
  pippo.Loop(leptype, mode, 0);
  pippo.Init(chainvub);
  pippo.Loop(leptype, mode, 1);
  pippo.FitMes();
  pippo.Finalize(mode);
}
