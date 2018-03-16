#include "FastGlauber/FastGlauber.h"
#include "configReader.h"
#include <stdlib.h>

using namespace std;
int main(int argc, char *argv[]){

    
    //read config file 
    string outputFileName, system, type, repulsionDistance; 
    int nevents = 0, energy = 0;
    configReader *reader = new configReader(argv[1]);
    vector< map<string, string> > config = reader->GetField();
    map<string, string>::iterator iter;
    bool flag = false;
    int isDeformed = 0;
    for(int i=0;i<config.size();i++){
        for(iter=config[i].begin();iter!=config[i].end();iter++){
            if((iter->first).find("FastGlauber")!=string::npos || flag){
                flag = true;
                if(iter->first.find("outputName")!=string::npos) outputFileName = iter->second;
                if(iter->first.find("system")!=string::npos) system = iter->second;
                if(iter->first.find("energy")!=string::npos) energy = atoi(iter->second.c_str());
                if(iter->first.find("type")!=string::npos) type = iter->second;
                if(iter->first.find("isDeformed")!=string::npos) isDeformed = atoi(iter->second);
                if(iter->first.find("nevents")!=string::npos) nevents = atoi(iter->second.c_str());
                if(iter->first.find("repulsionDistance")!=string::npos) repulsionDistance = iter->second;
            }
        }
        flag = false;
    }

    FastGlauber* maker = new FastGlauber(outputFileName, type, energy );
    //maker->Print("type");
    //maker->Run(nevents);
    //maker->Finish();



    return 0;
}

