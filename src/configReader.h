#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;
class configReader{

    public:
        configReader(char* input);
        ~configReader();
        vector< map<string, string> > GetField();

    private:

        bool readFile(char* input);
        vector< map<string, string> > Field;
};

inline vector< map<string, string> > configReader::GetField() { return Field;}


configReader::configReader(char* input) : Field(0){
    if(readFile(input)) 
    cout<<"Read config file succeed!"<<endl;
}

bool configReader::readFile(char* input){

    fstream inputFile;
    inputFile.open(input);
    if(!inputFile.is_open()) {
        cout<<"Config File Can not open!"<<endl;
        return false;
    }
    char tmp[1000];
    int n = -1;
    while(!inputFile.eof()){
        inputFile.getline(tmp, 1000);
        string line(tmp);
        if(line.find("[") != string::npos){
            map<string, string> FieldMap; 
            Field.push_back(FieldMap); 
            n++;
            Field[n].insert(pair<string, string>(line, line));
        }
        size_t pos = line.find("=");
        if(pos != string::npos) 
            if(pos != line.size()){
                string s1 = line.substr(0, pos);
                string s2 = line.substr(pos+1);
                Field[n].insert(pair<string, string>(s1, s2));
            }
    }
    return true;
}

#endif
