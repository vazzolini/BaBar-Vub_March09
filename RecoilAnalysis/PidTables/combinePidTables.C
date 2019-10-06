
// Origin of the PidTables (01/14/03)
// ----------------------------------
// Thorsten and Francesco had provided consistent and the latest tables, respectively. 
//
//   cp /nfs/farm/babar/AWG/PID/tables_new/2002-r10/electrons+/e.emcrad.*.urs 2002/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons+/e.emcrad.*.urs 2001/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons+/e.emcrad.*.urs 2000/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2002-r10/electrons-/e.emcrad.*.urs 2002/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons-/e.emcrad.*.urs 2001/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons-/e.emcrad.*.urs 2000/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables/2002-r10/electrons+/mu.mumug2.* 2002/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons+/mu.mumug2.* 2001/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons+/mu.mumug2.* 2000/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables/2002-r10/electrons-/mu.mumug2.* 2002/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons-/mu.mumug2.* 2001/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons-/mu.mumug2.* 2000/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2002-r10/electrons+/k.dstar.*.urs 2002/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons+/k.dstar.*.urs 2001/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons+/k.dstar.*.urs 2000/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2002-r10/electrons-/k.dstar.*.urs 2002/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons-/k.dstar.*.urs 2001/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons-/k.dstar.*.urs 2000/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2002-r10/electrons+/pi.ksTau31.*.urs 2002/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons+/pi.ksTau31.*.urs 2001/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons+/pi.ksTau31.*.urs 2000/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2002-r10/electrons-/pi.ksTau31.*.urs 2002/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons-/pi.ksTau31.*.urs 2001/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons-/pi.ksTau31.*.urs 2000/electrons-/

//   cp /nfs/farm/babar/AWG/PID/tables_new/2002-r10/electrons+/p*.urs 2002/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons+/p*.urs 2001/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons+/p*.urs 2000/electrons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2002-r10/electrons-/p*.urs 2002/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/electrons-/p*.urs 2001/electrons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/electrons-/p*.urs 2000/electrons-/

//   cp /nfs/farm/babar/AWG/PID/tables/2002-r10/muons+/* 2002/muons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/muons+/* 2001/muons+/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/muons+/* 2000/muons+/
//   cp /nfs/farm/babar/AWG/PID/tables/2002-r10/muons-/* 2002/muons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/muons-/* 2001/muons-/
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/muons-/* 2000/muons-/

//   cp /nfs/farm/babar/AWG/PID/tables/2002-r10/kaons+/*Micro* 2002/kaons+
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/kaons+/*Micro* 2001/kaons+
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/kaons+/*Micro* 2000/kaons+
//   cp /nfs/farm/babar/AWG/PID/tables/2002-r10/kaons-/*Micro* 2002/kaons-
//   cp /nfs/farm/babar/AWG/PID/tables_new/2001-r10/kaons-/*Micro* 2001/kaons-
//   cp /nfs/farm/babar/AWG/PID/tables_new/2000-r10/kaons-/*Micro* 2000/kaons-




void electrons() {

  // -- Electron efficiency
  cout << endl;
  PidTable a1("/u/ec/ursl/d/PidTables/2002/electrons+/e.emcrad.LH.tight.urs");
  PidTable b1("/u/ec/ursl/d/PidTables/2001/electrons+/e.emcrad.LH.tight.urs");
  PidTable c1("/u/ec/ursl/d/PidTables/2000/electrons+/e.emcrad.LH.tight.urs");
  a1.combine(b1); 
  a1.combine(c1); 
  a1.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons+/e.emcrad.LH.tight.urs");
 
  cout << endl;
  PidTable a2("/u/ec/ursl/d/PidTables/2002/electrons-/e.emcrad.LH.tight.urs");
  PidTable b2("/u/ec/ursl/d/PidTables/2001/electrons-/e.emcrad.LH.tight.urs");
  PidTable c2("/u/ec/ursl/d/PidTables/2000/electrons-/e.emcrad.LH.tight.urs");
  a2.combine(b2); 
  a2.combine(c2); 
  a2.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons-/e.emcrad.LH.tight.urs");

  // -- Pion misid
  cout << endl;
  PidTable a3("/u/ec/ursl/d/PidTables/2002/electrons+/pi.ksTau31.LH.tight.urs");
  PidTable b3("/u/ec/ursl/d/PidTables/2001/electrons+/pi.ksTau31.LH.tight.urs");
  PidTable c3("/u/ec/ursl/d/PidTables/2000/electrons+/pi.ksTau31.LH.tight.urs");
  a3.combine(b3); 
  a3.combine(c3); 
  a3.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons+/pi.ksTau31.LH.tight.urs");

  cout << endl;
  PidTable a4("/u/ec/ursl/d/PidTables/2002/electrons-/pi.ksTau31.LH.tight.urs");
  PidTable b4("/u/ec/ursl/d/PidTables/2001/electrons-/pi.ksTau31.LH.tight.urs");
  PidTable c4("/u/ec/ursl/d/PidTables/2000/electrons-/pi.ksTau31.LH.tight.urs");
  a4.combine(b4); 
  a4.combine(c4); 
  a4.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons-/pi.ksTau31.LH.tight.urs");


  // -- Kaon misid
  cout << endl;
  PidTable a5("/u/ec/ursl/d/PidTables/2002/electrons+/k.dstar.LH.tight.urs");
  PidTable b5("/u/ec/ursl/d/PidTables/2001/electrons+/k.dstar.LH.tight.urs");
  PidTable c5("/u/ec/ursl/d/PidTables/2000/electrons+/k.dstar.LH.tight.urs");
  a5.combine(b5); 
  a5.combine(c5); 
  a5.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons+/k.dstar.LH.tight.urs");

  cout << endl;
  PidTable a6("/u/ec/ursl/d/PidTables/2002/electrons-/k.dstar.LH.tight.urs");
  PidTable b6("/u/ec/ursl/d/PidTables/2001/electrons-/k.dstar.LH.tight.urs");
  PidTable c6("/u/ec/ursl/d/PidTables/2000/electrons-/k.dstar.LH.tight.urs");
  a6.combine(b6); 
  a6.combine(c6); 
  a6.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons-/k.dstar.LH.tight.urs");


  // -- Proton misid
  cout << endl;
  PidTable a7("/u/ec/ursl/d/PidTables/2002/electrons+/p.lambdaproton.LH.tight.urs");
  PidTable b7("/u/ec/ursl/d/PidTables/2001/electrons+/p.lambdaproton.LH.tight.urs");
  PidTable c7("/u/ec/ursl/d/PidTables/2000/electrons+/p.lambdaproton.LH.tight.urs");
  a7.combine(b7); 
  a7.combine(c7); 
  a7.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons+/p.lambdaproton.LH.tight.urs");

  cout << endl;
  PidTable a8("/u/ec/ursl/d/PidTables/2002/electrons-/p.lambdaproton.LH.tight.urs");
  PidTable b8("/u/ec/ursl/d/PidTables/2001/electrons-/p.lambdaproton.LH.tight.urs");
  PidTable c8("/u/ec/ursl/d/PidTables/2000/electrons-/p.lambdaproton.LH.tight.urs");
  a8.combine(b8); 
  a8.combine(c8); 
  a8.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons-/p.lambdaproton.LH.tight.urs");


  // -- Muon misid
  cout << endl;
  PidTable a9("/u/ec/ursl/d/PidTables/2002/electrons+/mu.mumug2.LH.Tight");
  PidTable b9("/u/ec/ursl/d/PidTables/2001/electrons+/mu.mumug2.LH.Tight");
  PidTable c9("/u/ec/ursl/d/PidTables/2000/electrons+/mu.mumug2.LH.Tight");
  a9.combine(b9); 
  a9.combine(c9); 
  a9.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons+/mu.mumug2.LH.Tight");

  cout << endl;
  PidTable a10("/u/ec/ursl/d/PidTables/2002/electrons-/mu.mumug2.LH.Tight");
  PidTable b10("/u/ec/ursl/d/PidTables/2001/electrons-/mu.mumug2.LH.Tight");
  PidTable c10("/u/ec/ursl/d/PidTables/2000/electrons-/mu.mumug2.LH.Tight");
  a10.combine(b10); 
  a10.combine(c10); 
  a10.dumpToFile("/u/ec/ursl/d/PidTables/200x/electrons-/mu.mumug2.LH.Tight");

}  


// ----------------------------------------------------------------------
void muons() {

  // -- Muon efficiency
  cout << endl;
  PidTable a1("/u/ec/ursl/d/PidTables/2002/muons+/mu.eemumu.Micro.Tight");
  PidTable b1("/u/ec/ursl/d/PidTables/2001/muons+/mu.eemumu.Micro.Tight");
  PidTable c1("/u/ec/ursl/d/PidTables/2000/muons+/mu.eemumu.Micro.Tight");
  a1.combine(b1); 
  a1.combine(c1); 
  a1.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons+/mu.eemumu.Micro.Tight");

  cout << endl;
  PidTable a2("/u/ec/ursl/d/PidTables/2002/muons-/mu.eemumu.Micro.Tight");
  PidTable b2("/u/ec/ursl/d/PidTables/2001/muons-/mu.eemumu.Micro.Tight");
  PidTable c2("/u/ec/ursl/d/PidTables/2000/muons-/mu.eemumu.Micro.Tight");
  a2.combine(b2); 
  a2.combine(c2); 
  a2.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons-/mu.eemumu.Micro.Tight");

  // -- Pion misid
  cout << endl;
  PidTable a3("/u/ec/ursl/d/PidTables/2002/muons+/pi.Dstar.Micro.Tight");
  PidTable b3("/u/ec/ursl/d/PidTables/2001/muons+/pi.Dstar.Micro.Tight");
  PidTable c3("/u/ec/ursl/d/PidTables/2000/muons+/pi.Dstar.Micro.Tight");
  a3.combine(b3); 
  a3.combine(c3); 
  a3.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons+/pi.Dstar.Micro.Tight");

  cout << endl;
  PidTable a4("/u/ec/ursl/d/PidTables/2002/muons-/pi.Dstar.Micro.Tight");
  PidTable b4("/u/ec/ursl/d/PidTables/2001/muons-/pi.Dstar.Micro.Tight");
  PidTable c4("/u/ec/ursl/d/PidTables/2000/muons-/pi.Dstar.Micro.Tight");
  a4.combine(b4); 
  a4.combine(c4); 
  a4.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons-/pi.Dstar.Micro.Tight");

  // -- Kaon misid
  cout << endl;
  PidTable a5("/u/ec/ursl/d/PidTables/2002/muons+/k.Dstar.Micro.Tight");
  PidTable b5("/u/ec/ursl/d/PidTables/2001/muons+/k.Dstar.Micro.Tight");
  PidTable c5("/u/ec/ursl/d/PidTables/2000/muons+/k.Dstar.Micro.Tight");
  a5.combine(b5); 
  a5.combine(c5); 
  a5.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons+/k.Dstar.Micro.Tight");

  cout << endl;
  PidTable a6("/u/ec/ursl/d/PidTables/2002/muons-/k.Dstar.Micro.Tight");
  PidTable b6("/u/ec/ursl/d/PidTables/2001/muons-/k.Dstar.Micro.Tight");
  PidTable c6("/u/ec/ursl/d/PidTables/2000/muons-/k.Dstar.Micro.Tight");
  a6.combine(b6); 
  a6.combine(c6); 
  a6.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons-/k.Dstar.Micro.Tight");

  // -- Proton misid
  cout << endl;
  PidTable a7("/u/ec/ursl/d/PidTables/2002/muons+/p.Lambda.Micro.Tight");
  PidTable b7("/u/ec/ursl/d/PidTables/2001/muons+/p.Lambda.Micro.Tight");
  PidTable c7("/u/ec/ursl/d/PidTables/2000/muons+/p.Lambda.Micro.Tight");
  a7.combine(b7); 
  a7.combine(c7); 
  a7.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons+/p.Lambda.Micro.Tight");

  cout << endl;
  PidTable a8("/u/ec/ursl/d/PidTables/2002/muons-/p.Lambda.Micro.Tight");
  PidTable b8("/u/ec/ursl/d/PidTables/2001/muons-/p.Lambda.Micro.Tight");
  PidTable c8("/u/ec/ursl/d/PidTables/2000/muons-/p.Lambda.Micro.Tight");
  a8.combine(b8); 
  a8.combine(c8); 
  a8.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons-/p.Lambda.Micro.Tight");

  // -- Dummy
  cout << endl;
  PidTable a9("/u/ec/ursl/d/PidTables/2002/muons+/mu.eemumu.Micro.Tight");
  PidTable b9("/u/ec/ursl/d/PidTables/2001/muons+/mu.eemumu.Micro.Tight");
  PidTable c9("/u/ec/ursl/d/PidTables/2000/muons+/mu.eemumu.Micro.Tight");
  a9.combine(b9); 
  a9.combine(c9); 
  a9.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons+/e.bla.Micro.Tight");

  cout << endl;
  PidTable a10("/u/ec/ursl/d/PidTables/2002/muons-/mu.eemumu.Micro.Tight");
  PidTable b10("/u/ec/ursl/d/PidTables/2001/muons-/mu.eemumu.Micro.Tight");
  PidTable c10("/u/ec/ursl/d/PidTables/2000/muons-/mu.eemumu.Micro.Tight");
  a10.combine(b10); 
  a10.combine(c10); 
  a10.dumpToFile("/u/ec/ursl/d/PidTables/200x/muons-/e.bla.Micro.Tight");

}  


// ----------------------------------------------------------------------
void kaons() {

  // -- Kaon efficiency
  cout << endl;
  PidTable a1("/u/ec/ursl/d/PidTables/2002/kaons+/k.Dstar.Micro.Tight");
  PidTable b1("/u/ec/ursl/d/PidTables/2001/kaons+/k.Dstar.Micro.Tight");
  PidTable c1("/u/ec/ursl/d/PidTables/2000/kaons+/k.Dstar.Micro.Tight");
  a1.combine(b1); 
  a1.combine(c1); 
  a1.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons+/k.Dstar.Micro.Tight");

  cout << endl;
  PidTable a2("/u/ec/ursl/d/PidTables/2002/kaons-/k.Dstar.Micro.Tight");
  PidTable b2("/u/ec/ursl/d/PidTables/2001/kaons-/k.Dstar.Micro.Tight");
  PidTable c2("/u/ec/ursl/d/PidTables/2000/kaons-/k.Dstar.Micro.Tight");
  a2.combine(b2); 
  a2.combine(c2); 
  a2.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons-/k.Dstar.Micro.Tight");


  // -- Pion misid
  cout << endl;
  PidTable a3("/u/ec/ursl/d/PidTables/2002/kaons+/pi.Dstar.Micro.Tight");
  PidTable b3("/u/ec/ursl/d/PidTables/2001/kaons+/pi.Dstar.Micro.Tight");
  PidTable c3("/u/ec/ursl/d/PidTables/2000/kaons+/pi.Dstar.Micro.Tight");
  a3.combine(b3); 
  a3.combine(c3); 
  a3.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons+/pi.Dstar.Micro.Tight");

  cout << endl;
  PidTable a4("/u/ec/ursl/d/PidTables/2002/kaons-/pi.Dstar.Micro.Tight");
  PidTable b4("/u/ec/ursl/d/PidTables/2001/kaons-/pi.Dstar.Micro.Tight");
  PidTable c4("/u/ec/ursl/d/PidTables/2000/kaons-/pi.Dstar.Micro.Tight");
  a4.combine(b4); 
  a4.combine(c4); 
  a4.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons-/pi.Dstar.Micro.Tight");


  // -- Electron misid
  cout << endl;
  PidTable a5("/u/ec/ursl/d/PidTables/2002/kaons+/e.emcrad.Micro.Tight");
  PidTable b5("/u/ec/ursl/d/PidTables/2001/kaons+/e.emcrad.Micro.Tight");
  PidTable c5("/u/ec/ursl/d/PidTables/2000/kaons+/e.emcrad.Micro.Tight");
  a5.combine(b5); 
  a5.combine(c5); 
  a5.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons+/e.emcrad.Micro.Tight");

  cout << endl;
  PidTable a6("/u/ec/ursl/d/PidTables/2002/kaons-/e.emcrad.Micro.Tight");
  PidTable b6("/u/ec/ursl/d/PidTables/2001/kaons-/e.emcrad.Micro.Tight");
  PidTable c6("/u/ec/ursl/d/PidTables/2000/kaons-/e.emcrad.Micro.Tight");
  a6.combine(b6); 
  a6.combine(c6); 
  a6.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons-/e.emcrad.Micro.Tight");


  // -- Proton misid
  cout << endl;
  PidTable a7("/u/ec/ursl/d/PidTables/2002/kaons+/p.Lambda.Micro.Tight");
  PidTable b7("/u/ec/ursl/d/PidTables/2001/kaons+/p.Lambda.Micro.Tight");
  PidTable c7("/u/ec/ursl/d/PidTables/2000/kaons+/p.Lambda.Micro.Tight");
  a7.combine(b7); 
  a7.combine(c7); 
  a7.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons+/p.Lambda.Micro.Tight");

  cout << endl;
  PidTable a8("/u/ec/ursl/d/PidTables/2002/kaons-/p.Lambda.Micro.Tight");
  PidTable b8("/u/ec/ursl/d/PidTables/2001/kaons-/p.Lambda.Micro.Tight");
  PidTable c8("/u/ec/ursl/d/PidTables/2000/kaons-/p.Lambda.Micro.Tight");
  a8.combine(b8); 
  a8.combine(c8); 
  a8.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons-/p.Lambda.Micro.Tight");


  // -- Muon misid
  cout << endl;
  PidTable a9("/u/ec/ursl/d/PidTables/2002/kaons+/mu.eemumu.Micro.Tight");
  PidTable b9("/u/ec/ursl/d/PidTables/2001/kaons+/mu.eemumu.Micro.Tight");
  PidTable c9("/u/ec/ursl/d/PidTables/2000/kaons+/mu.eemumu.Micro.Tight");
  a9.combine(b9); 
  a9.combine(c9); 
  a9.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons+/mu.eemumu.Micro.Tight");

  cout << endl;
  PidTable a10("/u/ec/ursl/d/PidTables/2002/kaons-/mu.eemumu.Micro.Tight");
  PidTable b10("/u/ec/ursl/d/PidTables/2001/kaons-/mu.eemumu.Micro.Tight");
  PidTable c10("/u/ec/ursl/d/PidTables/2000/kaons-/mu.eemumu.Micro.Tight");
  a10.combine(b10);
  a10.combine(c10);
  a10.dumpToFile("/u/ec/ursl/d/PidTables/200x/kaons-/mu.eemumu.Micro.Tight");

}  
