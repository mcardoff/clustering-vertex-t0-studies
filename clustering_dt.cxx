#include <RtypesCore.h>
#include <TROOT.h>
#include "clustering_constants.h"
#include "event_processing.h"
#include "plotting_utilities.h"

using namespace MyUtl;

void clustering_dt() {
  gStyle->SetOptStat(0);

  TChain chain ("ntuple");
  setupChain(chain, "./ntuple/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  ROOT::EnableImplicitMT(); // uses all cores

  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);

  // HGTD Times
  std::map<Score, AnalysisObj> hgtdtimesMap;
  hgtdtimesMap.emplace(Score::HGTD  , AnalysisObj("hgtdtimes", "HGTD Times", Score::HGTD)  );
  hgtdtimesMap.emplace(Score::TRKPTZ, AnalysisObj("hgtdtimes", "HGTD Times", Score::TRKPTZ));
  hgtdtimesMap.emplace(Score::PASS  , AnalysisObj("hgtdtimes", "HGTD Times", Score::PASS)  );

  gErrorIgnoreLevel = kWarning;
  std::cout << "Starting Event Loop\n";
  bool progress = true;
  Long64_t nEvent = chain.GetEntries();
  Long64_t evtMax = nEvent+1;
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t thisEvnt = chain.GetReadEntry() - chain.GetChainOffset(); // +1 bc its 0 indexed
    Long64_t readNum = chain.GetReadEntry()+1;

    // if (hgtdtimesMap[Score::TRKPTZ].get("pufrac").effTotal->Integral() > 10e3)
    //   break;
    
    if (progress && readNum % 1000 == 0)
      std::cout << "Progress: " << readNum << "/" << nEvent << "\n";
    auto passHgtd = processEventData(&branch, false, true, false, hgtdtimesMap);
  }

  auto maps = {&hgtdtimesMap};

  for (auto &map: maps)
    for (auto &[k, analysis] : *map)
      analysis.postProcessing();
      
  for (auto &map: maps)
    for (auto &[k, analysis] : *map)
      analysis["hs_track"]->printEfficiencyStats();

  std::cout << "FINISHED PROCESSING\n";

  for (const auto *const KEY: {"fjet", "ftrack", "pu_frac", "hs_track", "pu_track"}) {
    moneyPlot(Form("./hgtd_trkptz_%s.pdf", KEY), KEY, canvas,
	      {&hgtdtimesMap.at(HGTD), &hgtdtimesMap.at(TRKPTZ),
	      });
  }

  std::vector<AnalysisObj*> plts =
    {&hgtdtimesMap.at(HGTD), &hgtdtimesMap.at(TRKPTZ), &hgtdtimesMap.at(PASS)};
  inclusivePlot(Form("./inclusivereso_logscale.pdf"), true, true,
		-400, 400, canvas, plts);

  inclusivePlot(Form("./inclusivereso_linscale.pdf"), false, true,
		-200, 200, canvas, plts);
  
  std::cout << "FINISHED PLOT PRINTING\n";  
  return;
}
