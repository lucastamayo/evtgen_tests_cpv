#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "Pythia8Plugins/EvtGen.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

#include <iostream>
#include <string>

using namespace Pythia8;

int main() {

    const int nEvents = 50000;

    // -----------------------------
    // Arquivos de entrada/saída
    // -----------------------------
    const std::string mainDecayFile = "DECAY_2014_NOLONGLIFE.DEC";
    const std::string userDecayFile = "Bd_JpsiKs_SSDCP.dec";
    const std::string pdtFile       = "evt.pdl";

    // salva direto no EOS
    const std::string outputFile    = "events_SSDCP2.hepmc";

    // -----------------------------
    // Inicializa Pythia
    // -----------------------------
    Pythia pythia;

    pythia.readString("Beams:eCM = 13600.");

    pythia.readString("HardQCD:all = off");
    pythia.readString("HardQCD:gg2bbbar = on");
    pythia.readString("HardQCD:qqbar2bbbar = on");

    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("PhaseSpace:pTHatMin = 8.0");

    // turning off B mixing in Pythia
    pythia.readString("ParticleDecays:mixB = off");

    // turning off b decays in pythia
    pythia.readString("511:mayDecay = off");
    pythia.readString("-511:mayDecay = off");

    pythia.init();

    // -----------------------------
    // calling evtgen
    // -----------------------------
    // mixing = 1  -> incoherent B mixing
    EvtGenDecays evtgenDecays(
        &pythia,
        mainDecayFile,
        pdtFile,
        nullptr,
        nullptr,
        1
    );

    evtgenDecays.readDecayFile(userDecayFile);

    // -----------------------------
    // to hepmc
    // -----------------------------
    HepMC::Pythia8ToHepMC toHepMC;
    HepMC::IO_GenEvent ascii_io(outputFile.c_str(), std::ios::out);

    int nAccepted = 0;

    // -----------------------------
    // Loop de eventos
    // -----------------------------
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

        if (!pythia.next()) continue;

        // aplica EvtGen no record do Pythia:
        // decay + lifetime + mixing + atualização do event record
        evtgenDecays.decay();

        // converte evento final do Pythia para HepMC
        HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
        toHepMC.fill_next_event(pythia, hepmcevt);

        // grava no arquivo
        ascii_io << hepmcevt;
        delete hepmcevt;

        ++nAccepted;

        if (nAccepted <= 5) {
            std::cout << "Processed event " << nAccepted << std::endl;
        }
    }

    // estatísticas finais
    pythia.stat();

    std::cout << "\nFinished." << std::endl;
    std::cout << "Accepted events: " << nAccepted << std::endl;
    std::cout << "Output written to: " << outputFile << std::endl;

    return 0;
}