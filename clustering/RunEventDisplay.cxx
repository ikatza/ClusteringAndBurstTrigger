#include "PhotonTimingDisplay.hh"
#include "PhotonWirePositionDisplay.hh"
#include "WirePositionTimingDisplay.hh"
// #include "WireTimingDisplay.hh"
#include <iostream>


void PrintUsage() {
  std::cout << "Option n is the event number"           << std::endl;
  std::cout << "Option i is the file name"              << std::endl;
  std::cout << "Option t is the tree name"              << std::endl;
  std::cout << "Option e is the event type"             << std::endl;
  std::cout << " -e 0 Non BT"                           << std::endl;
  std::cout << " -e 1 SNnu"                             << std::endl;
  std::cout << " -e 2 APA"                              << std::endl;
  std::cout << " -e 3 CPA"                              << std::endl;
  std::cout << " -e 4 Ar39"                             << std::endl;
  std::cout << " -e 5 Neutron"                          << std::endl;
  std::cout << " -e 6 Krypton"                          << std::endl;
  std::cout << " -e 7 Polonium"                         << std::endl;
  std::cout << " -e 8 Radon"                            << std::endl;
  std::cout << " -e 9 Ar42"                             << std::endl;
  std::cout << "Option d can be:"                       << std::endl;
  std::cout << " -d 0 Photon Timing display"            << std::endl;
  std::cout << " -d 1 Photon and wire position display" << std::endl;
  std::cout << " -d 2 Wire Position and Timing display" << std::endl;
  std::cout << " -d 3 Wire Timing display"              << std::endl;
  exit(1);
};

int main(int argc, char** argv){
  int opt;
  // Shut GetOpt error messages down (return '?'): 
  extern char *optarg;
  extern int  optopt;

  int nEvent=-1;
  int Type=-1;
  int Display=-1;
  std::string FileName="";
  std::string TreeName="";
  while ( (opt = getopt(argc, argv, "n:e:d:i:t:")) != -1 ) {  // for each option...
    switch ( opt ) {
    case 'n':
      nEvent = atoi(optarg);
      break;
    case 't':
      TreeName = optarg;
      break;
    case 'e':
      Type = atoi(optarg);
      break;
    case 'i':
      FileName = optarg;
      break;
    case 'd':
      Display = atoi(optarg);
      break;
    case '?':  // unknown option...
      std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
      PrintUsage();
      break;
    }
  }
  
  switch(Display){
  case 0:{
    PhotonTimingDisplay ptd(FileName,TreeName);
    ptd.DisplayEvent(nEvent,Type);
    return 0;
  }
  case 1:{
    PhotonWirePositionDisplay pwpd(FileName,TreeName);
    pwpd.DisplayEvent(nEvent,Type);
    return 0;
  }
  case 2:{
    WirePositionTimingDisplay wptd(FileName,TreeName);
    wptd.DisplayEvent(nEvent,Type);
    return 0;
  }
  default:
    std::cout << "The requested event display doesn't exit (or hasn't been implemented yet!)"
              << std::endl;
    return 1;
  }


};

