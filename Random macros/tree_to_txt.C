#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
using namespace std;

void tree_to_txt(){

TFile *f1=new TFile("LHC15o_Train330.root"); //"tree_data.root"
TTree *t1;
f1->GetObject("jRecTree",t1); //"myTree"
TBranch *b1 = t1->GetBranch(""dimuonBranch""); //"dimuonBranch"
// TBranch *b2 = t1->GetBranch("jRecY");
// TBranch *b3 = t1->GetBranch("jRecPt");
TLeaf *l1a=b1->GetLeaf("m");
TLeaf *l1b=b1->GetLeaf("y");
TLeaf *l1c=b1->GetLeaf("pt");

Float_t m1;
Float_t y1;
Float_t pt1;

l1a->SetAddress(&m1);
l1b->SetAddress(&y1);
l1c->SetAddress(&pt1);

// b1->SetAddress(&m1);
// b2->SetAddress(&y1);
// b3->SetAddress(&pt1);

Int_t c1=0;
//Get the entries and fill them to hist

  ofstream myfile;
  myfile.open ("data.txt");
 // myfile << "m           y           pt " << endl;


for(Int_t i=0;i<10;i++){ //loop over events
    t1->GetEntry(i);
   myfile << m1 << " " << y1 << " " << pt1<< endl; //write to file
  }
  myfile.close();

}

