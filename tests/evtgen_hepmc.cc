#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

#include <iostream>

// -----------------------------------------
// Conversão recursiva EvtGen → HepMC
// -----------------------------------------
void fillHepMC(EvtParticle* p, HepMC::GenVertex* vtx,
               HepMC::GenEvent* evt) {

    if (!p) return;

    // -----------------------------
    // Momento
    // -----------------------------
    EvtVector4R p4 = p->getP4();
    HepMC::FourVector mom(p4.get(1), p4.get(2), p4.get(3), p4.get(0));

    int pdg = EvtPDL::getStdHep(p->getId());
    int status = (p->getNDaug() == 0) ? 1 : 2;

    HepMC::GenParticle* hp = new HepMC::GenParticle(mom, pdg, status);

    // -----------------------------
    // 🔴 POSIÇÃO DO VÉRTICE (ESSENCIAL)
    // -----------------------------
    EvtVector4R x4 = p->get4Pos();

    HepMC::FourVector pos(
        x4.get(1),  // x
        x4.get(2),  // y
        x4.get(3),  // z
        x4.get(0)   // t
    );

    // -----------------------------
    // vértice de produção
    // -----------------------------
    if (!vtx) {
        vtx = new HepMC::GenVertex(pos);
        evt->add_vertex(vtx);
    }

    vtx->add_particle_out(hp);

    // -----------------------------
    // vértice de decaimento
    // -----------------------------
    if (p->getNDaug() > 0) {

        HepMC::GenVertex* decayVtx = new HepMC::GenVertex(pos);
        evt->add_vertex(decayVtx);

        decayVtx->add_particle_in(hp);

        for (int i = 0; i < p->getNDaug(); i++) {
            fillHepMC(p->getDaug(i), decayVtx, evt);
        }
    }
}

// -----------------------------------------
// MAIN
// -----------------------------------------
int main() {

    std::string decayFile  = "Bd_JpsiKs_SSDCPB.dec";
    std::string pdtFile    = "evt.pdl";
    std::string outputFile = "events_evtgen.hepmc";

    // -----------------------------
    // Random engine
    // -----------------------------
    EvtRandomEngine* randomEngine = new EvtSimpleRandomEngine();

    // -----------------------------
    // Inicializa EvtGen
    // -----------------------------
    EvtGen evtgen(decayFile.c_str(), pdtFile.c_str(), randomEngine);

    // -----------------------------
    // Upsilon(4S) em repouso
    // -----------------------------
    EvtId upsId = EvtPDL::getId("Upsilon(4S)");
    double mUps = EvtPDL::getMass(upsId);

    EvtVector4R pInit(mUps, 0.0, 0.0, 0.0);

    // -----------------------------
    // HepMC output
    // -----------------------------
    HepMC::IO_GenEvent ascii_io(outputFile.c_str(), std::ios::out);

    // -----------------------------
    // Loop de eventos
    // -----------------------------
    for (int i = 0; i < 100000; i++) {

        // cria Upsilon(4S)
        EvtParticle* ups = EvtParticleFactory::particleFactory(upsId, pInit);

        // decai tudo (mixing + SSD_CP)
        ups->decay();

        // debug
        if (i < 3) {
            std::cout << "\n==== EVENTO " << i << " ====\n";
            ups->printTree();
        }

        // cria evento HepMC
        HepMC::GenEvent* evt = new HepMC::GenEvent();

        // converte árvore
        fillHepMC(ups, nullptr, evt);

        // escreve no arquivo
        ascii_io << evt;

        delete evt;
        ups->deleteTree();
    }

    std::cout << "\nFinished.\n";
    std::cout << "Output: " << outputFile << std::endl;

    return 0;
}