#ifndef EVENTINSPECTOR_HH
#define EVENTINSPECTOR_HH


#include "Helper.h"
#include "InputManager.h"
#include <unistd.h>



class EventInspector {
protected:
  InputManager im;
  virtual void ExecuteInLoop() = 0;
  virtual bool Selected() = 0;
  void InitialiseInput(){
    std::cout << "Setting up the input tree" << std::endl;
    im.SetInputFile(InFileName.c_str());
    im.SetInputTree(InTreeName.c_str());
    im.LoadTree();
    std::cout << "Done initialising!" << std::endl;
  };
  
private:
  std::string OutFileName;
  std::string InFileName;
  std::string InTreeName;
  int LoopOverNEvents;
  TCanvas c;
  
protected:
  EventInspector(int argc, char** argv) {
    int opt;
    // Shut GetOpt error messages down (return '?'): 
    extern char *optarg;
    extern int  optopt;

    OutFileName = "Event_Info.pdf";
    InFileName  = "/dune/data/users/plasorak/anahist.root";
    InTreeName  = "snanagaushit/SNSimTree";

    while ( (opt = getopt(argc, argv, "n:o:t:i:")) != -1 ) {  // for each option...
      switch ( opt ) {
      case 'n':
        LoopOverNEvents = atoi(optarg);
        break;
      case 't':
        InTreeName = std::string(optarg);
        break;
      case 'i':
        InFileName = std::string(optarg);
        break;
      case 'o':
        OutFileName = std::string(optarg);
        break;
      case '?':  // unknown option...
        std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
        break;
      }
    }
    InitialiseInput();

    
  };
  
  void SetOutFileName(std::string s) {
    if (OutFileName == "Event_Info.pdf") {
      OutFileName = s;
    } else {
      std::cout << "Not overriding the name of the file since it was changed in the command line." << std::endl;
      std::cout << "The file name is " << OutFileName << std::endl;
    }
  };
  
  ~EventInspector() {
    std::cout << "Closing the file " << OutFileName << std::endl;
    c.Print((OutFileName+"]").c_str());
  };
  
public:
  void Loop() {
    if (LoopOverNEvents>0)
      LoopOverNEvents = std::min(LoopOverNEvents,(int)im.GetEntries());
    else
      LoopOverNEvents = im.GetEntries();
    for (int iEntry=0; iEntry<LoopOverNEvents; ++iEntry) {
      im.GetEntry(iEntry);
      PrintProgress(iEntry,LoopOverNEvents);
      if(!Selected()) continue;
      ExecuteInLoop();
    }
    std::cout << "Opening the pdf file for saving " << OutFileName << std::endl;
    c.Print((OutFileName+"[").c_str());
  };
  
  void AddPlot() {
    c.Print(OutFileName.c_str());
  }

};




#endif
