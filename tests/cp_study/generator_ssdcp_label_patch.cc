#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "EvtGenCleanLabels.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace Pythia8;

int main(int argc, char* argv[]) {
    int nEvents = 50000;
    if (argc > 1) {
        nEvents = std::atoi(argv[1]);
    }

    const std::string mainDecayFile = "DECAY_2014_NOLONGLIFE.DEC";
    const std::string userDecayFile = "Bd_JpsiKs_SSDCP.dec";
    const std::string pdtFile = "evt.pdl";
    const std::string outputFile = "cp_study/results_ssdcp_label_patch/events_SSDCP_label_patch.hepmc";
    const std::string truthFile = "cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv";

    Pythia pythia;
    pythia.readString("Beams:eCM = 13600.");
    pythia.readString("HardQCD:all = off");
    pythia.readString("HardQCD:gg2bbbar = on");
    pythia.readString("HardQCD:qqbar2bbbar = on");
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("PhaseSpace:pTHatMin = 8.0");
    pythia.readString("ParticleDecays:mixB = off");
    pythia.readString("511:mayDecay = off");
    pythia.readString("-511:mayDecay = off");
    pythia.init();

    EvtGenDecays evtgenDecays(
        &pythia,
        mainDecayFile,
        pdtFile,
        nullptr,
        nullptr,
        1
    );
    evtgenDecays.readDecayFile(userDecayFile);

    HepMC::Pythia8ToHepMC toHepMC;
    HepMC::IO_GenEvent ascii_io(outputFile.c_str(), std::ios::out);
    std::ofstream truth(truthFile);
    truth << "hepmc_event,py_index,parent_pdg,other_b_pdg,flip,t_mm,t_ps\n";

    int nAccepted = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;

        evtgenDecays.decay();
        for (const auto& row : evtgenDecays.ssdcpTruth) {
            truth << nAccepted << ","
                  << row.pyIndex << ","
                  << row.parentPdg << ","
                  << row.otherBPdg << ","
                  << row.flip << ","
                  << row.tMm << ","
                  << row.tMm / 0.299792458 << "\n";
        }

        HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
        toHepMC.fill_next_event(pythia, hepmcevt);
        ascii_io << hepmcevt;
        delete hepmcevt;

        ++nAccepted;
        if (nAccepted <= 5 || nAccepted % 10000 == 0) {
            std::cout << "Processed event " << nAccepted << std::endl;
        }
    }

    pythia.stat();

    std::cout << "\nFinished." << std::endl;
    std::cout << "Accepted events: " << nAccepted << std::endl;
    std::cout << "Output written to: " << outputFile << std::endl;
    std::cout << "Truth tags written to: " << truthFile << std::endl;
    return 0;
}
