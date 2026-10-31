#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/EvtGen.h"
#include "Pythia8Plugins/HepMC2.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

#include <iostream>
#include <string>

using namespace Pythia8;

int main() {
    const int nEvents = 10000;

    const std::string mainDecayFile = "DECAY_2014_NOLONGLIFE.DEC";
    const std::string userDecayFile = "Bd_JpsiKs_SSDCP.dec";
    const std::string pdtFile = "evt.pdl";
    const std::string outputFile = "cp_study/results_mixing0_test/events_SSDCP_mixing0.hepmc";

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
        0
    );
    evtgenDecays.readDecayFile(userDecayFile);

    HepMC::Pythia8ToHepMC toHepMC;
    HepMC::IO_GenEvent ascii_io(outputFile.c_str(), std::ios::out);

    int nAccepted = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;
        evtgenDecays.decay();
        HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
        toHepMC.fill_next_event(pythia, hepmcevt);
        ascii_io << hepmcevt;
        delete hepmcevt;
        ++nAccepted;
    }

    std::cout << "Accepted events: " << nAccepted << std::endl;
    std::cout << "Output written to: " << outputFile << std::endl;
    return 0;
}
