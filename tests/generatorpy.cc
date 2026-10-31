#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/EvtGen.h"

#include <iostream>
#include <fstream>
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

    const std::string outputFile    = "events_SSDCP2_pythia.txt";

    std::ofstream out(outputFile);

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

    // mixing feito no EvtGen
    pythia.readString("ParticleDecays:mixB = off");

    // impedir Pythia de decair B0
    pythia.readString("511:mayDecay = off");
    pythia.readString("-511:mayDecay = off");

    pythia.init();

    // -----------------------------
    // EvtGen
    // -----------------------------
    EvtGenDecays evtgenDecays(
        &pythia,
        mainDecayFile,
        pdtFile,
        nullptr,
        nullptr,
        1 // incoherent (pp)
    );

    evtgenDecays.readDecayFile(userDecayFile);

    int nAccepted = 0;

    // -----------------------------
    // Loop de eventos
    // -----------------------------
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

        if (!pythia.next()) continue;

        // aplica EvtGen
        evtgenDecays.decay();

        // =============================
        // SALVAR EVENTO
        // =============================
        out << "\n=====================================\n";
        out << "Event " << iEvent << "\n";

        // imprime o event record completo
        for (int i = 0; i < pythia.event.size(); ++i) {

    const Particle& p = pythia.event[i];

    out << i << " "
        << p.id() << " "
        << p.name() << " "
        << p.status() << " "
        << p.mother1() << " "
        << p.mother2() << " "
        << p.daughter1() << " "
        << p.daughter2() << " "
        << p.px() << " "
        << p.py() << " "
        << p.pz() << " "
        << p.e() << "\n";
}

        ++nAccepted;

        if (nAccepted % 1000 == 0) {
            std::cout << "Processed " << nAccepted << " events\n";
        }
    }

    // -----------------------------
    // Final
    // -----------------------------
    pythia.stat();

    std::cout << "\nFinished.\n";
    std::cout << "Accepted events: " << nAccepted << "\n";
    std::cout << "Output written to: " << outputFile << "\n";

    out.close();

    return 0;
}