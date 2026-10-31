#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#ifdef LABEL_PATCH
#include "EvtGenCleanLabels.h"
#else
#include "Pythia8Plugins/EvtGen.h"
#endif

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace Pythia8;

int main(int argc, char* argv[]) {
    if (argc != 8) {
        std::cerr << "Usage: generate SVS|SSD_CP events output.hepmc truth.csv seed mixing data_directory\n";
        return 2;
    }
    const std::string model = argv[1];
    const int nEvents = std::stoi(argv[2]);
    const int seed = std::stoi(argv[5]);
    const int mixing = std::stoi(argv[6]);
    if ((model != "SVS" && model != "SSD_CP") || nEvents <= 0 || seed < 1 || seed > 900000000 || (mixing != 0 && mixing != 1)) return 2;
    const std::string data = argv[7];
    const std::string mainDecayFile = data + "/DECAY_2014_NOLONGLIFE.DEC";
    const std::string userDecayFile = data + (model == "SVS" ? "/Bd_JpsiKs_SVS.dec" : "/Bd_JpsiKs_SSDCP.dec");
    const std::string pdtFile = data + "/evt.pdl";
    const std::string outputFile = argv[3];
    const std::string truthFile = argv[4];
    for (const auto& path : {mainDecayFile, userDecayFile, pdtFile}) {
        if (!std::ifstream(path)) { std::cerr << "Cannot read " << path << "\n"; return 2; }
    }
    if (std::ifstream(outputFile) || std::ifstream(truthFile)) {
        std::cerr << "Output already exists; choose a new directory.\n"; return 2;
    }
    { std::ofstream check(outputFile); if (!check) return 2; }
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
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));
    if (!pythia.init()) return 3;

    EvtGenDecays evtgenDecays(
        &pythia,
        mainDecayFile,
        pdtFile,
        nullptr,
        nullptr,
        mixing
    );
    evtgenDecays.readDecayFile(userDecayFile);

    HepMC::Pythia8ToHepMC toHepMC;
    HepMC::IO_GenEvent ascii_io(outputFile.c_str(), std::ios::out);
    std::ofstream truth(truthFile);
    if (!truth) return 2;
    truth << "hepmc_event,py_index,parent_pdg,other_b_pdg,flip,t_mm,t_ps\n";

    int nAccepted = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;

        evtgenDecays.decay();
#ifdef LABEL_PATCH
        for (const auto& row : evtgenDecays.ssdcpTruth) {
            truth << nAccepted << ","
                  << row.pyIndex << ","
                  << row.parentPdg << ","
                  << row.otherBPdg << ","
                  << row.flip << ","
                  << row.tMm << ","
                  << row.tMm / 0.299792458 << "\n";
        }

#endif
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
    return nAccepted > 0 ? 0 : 4;
}
