#include <iostream>
#include <fstream>


#include "TTree.h"
#include "TFile.h"
#include "TMath.h"


using namespace TMath;



int FIND_GG()
{
    auto tr_ph = new TChain( "tr_ph" );
    tr_ph->Add( "root://cmd//scan2019/scan2019_tr_ph_fc_e887.5_v8.root" );
    
    std::cout << "BEGIN?\n";
    
    //see -> https://cmd.inp.nsk.su/docdb/files/0003/000381/003/ryzhenenkov_lum2021.pdf
    //select gamma gamma
    const char selection[] = "(nt == 0) && (nph == 2)   && (Abs(Abs(phphi[1] - phphi[0]) - Pi()) < 0.2) && (Abs(phth[1] + phth[0] - Pi()) < 0.6) &&  (  (Pi() - 1.0 > (phth[0] +Pi() - phth[1])/2) && ((phth[0] +Pi() - phth[1])/2 > 1.0)  )   &&   (0.5*emeas < phen[0]) && (0.5*emeas < phen[1]) && (phen[0] < 1.5*emeas) && (phen[1] < 1.5*emeas) ";
    
    std::cout << "BEFORE?\n";
    
    tr_ph->Draw( "evnum:runnum",  selection);
    
    std::cout << "DRAW?\n";

    auto selnum = tr_ph->GetSelectedRows();
    auto evnum = tr_ph->GetV1();
    auto runnum = tr_ph->GetV2();

    std::cout << selnum << std::endl;
    std::cout << "ARRS\n";

    std::ofstream out("scan2019_tr_ph_fc_e887.5_v8.txt");
    
    std::cout << "OPEN\n";

    long counter = 0;
    for(int i=0; i<selnum; i++, counter++) {
        if (i==0) 
        {
            out << runnum[i] << "    :   \"";
            out << evnum[i];
            continue;
        }
        if (runnum[i] != runnum[i-1]) {
            out << "\"     -n     " << counter << "\n\n\n";
            counter = 0;
            out << runnum[i] <<  "   :   \"";
            out << evnum[i];
        } else {
            out <<  ';' <<  evnum[i];
        }
    }
    out << "\"     -n     " << counter << "\n\n\n";
    out.close();

    std::cout << "END\n";
    return 0;
}
