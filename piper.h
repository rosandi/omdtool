#include <iostream>
#include <string>
#include <cstdlib>

#ifndef __PIPER_H__
#define __PIPER_H__

using std::string;
namespace omd {
  
  class Piper {
    FILE* extpipe;

  public:

    Piper(const char *prg=NULL) {
      if(prg) open(prg);
    }
    
    void open(const char* prg) {
      if((extpipe=popen(prg, "w"))==NULL) throw "Cannot open pipe\n";
    }

    virtual ~Piper() {
      pclose(extpipe);
    }
      
    void send(const string cmd) {
      fprintf(extpipe, "%s\n", cmd.c_str());
      fflush(extpipe);
    }

  };
  
}

#endif
