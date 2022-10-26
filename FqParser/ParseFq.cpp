#include "ParseFq.h"
#include <fstream>
#include <string>
#include <iostream>
#include <array>
#include <map>
#include <stdlib.h>
#include <iterator>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <set>

using std::istream_iterator;
using namespace std;

map<string,string> argsMap;
std::map<std::string, std::string>outFileNameMap;

ParseFq::ParseFq(const map<string,string> &argsMap)
{
    this->argsMap = argsMap;
}

ParseFq::~ParseFq()
{
    //dtor
}

void ParseFq::ReadFq(string fileName,string file_Type){
    std::ofstream outfile;
    //cout << "Given File to load "+fileName <<endl;
    string out_File_Name = fileName;
    const size_t last_slash_idx = out_File_Name.find_last_of("\\/");
    if (std::string::npos != last_slash_idx)
    {
        out_File_Name.erase(0, last_slash_idx + 1);
    }
    
    out_File_Name = argsMap.find("out")->second+"tmp_"+out_File_Name;
    //cout<<out_File_Name<<endl;
    
    outFileNameMap.emplace(file_Type,out_File_Name);
    outfile.open(out_File_Name, std::ios_base::app);
    
    ifstream file(fileName);//("E:\\Data\\NGS-16\\DZA02_S2_L001_R1_001.fastq");
    string str;
    int lineCounter = 0;
    std::array<string,4> readAr;
    int readCounter = 0;
    while (std::getline(file, str))
    {
        if(readCounter<=10000000){
            if(lineCounter==4){
                lineCounter=0;
                readCounter++;
                //cout << readCounter <<endl;
                this->ProcessRead(readAr,outfile);
                readAr.empty();
                readAr[lineCounter] = str;
                lineCounter++;
            }
            else{
                readAr[lineCounter] = str;
                lineCounter++;
            }
        }
    }
    cout << readCounter <<endl;
    file.close();
    outfile.close();
}

void ParseFq::ProcessRead(std::array<string,4> readAr,std::ofstream &outfile){
    std::vector<string> readVect;
    string seq = readAr[1];
    string qualSeq = readAr[3];
    bool doSave = true;
    if(argsMap.count("lcrop")>0){
        string numVal = argsMap.find("lcrop")->second;
        int cLen = atoi(numVal.c_str());
        seq = seq.substr(cLen);
        qualSeq = qualSeq.substr(cLen);
        //cout<<"1"+seq<<endl;
    }
    if(argsMap.count("tcrop")>0){
        string numVal = argsMap.find("tcrop")->first;
        int cLen = atoi(numVal.c_str());
        int lastIndex = seq.length()-cLen;
        if(lastIndex<seq.length()){
            seq = seq.substr(0,lastIndex);
            qualSeq = qualSeq.substr(0,lastIndex);
            //cout<<"2"+seq<<endl;
        }
    }
    if(argsMap.count("leading")>0){
        string numVal = argsMap.find("leading")->second;
        int qual = atoi(numVal.c_str());
        QualityTrim(seq,qualSeq,1,qual);
        //cout<<"3"+seq<<endl;
    }
    if(argsMap.count("trailing")>0){
        string numVal = argsMap.find("trailing")->second;
        int qual = atoi(numVal.c_str());
        QualityTrim(seq,qualSeq,2,qual);
        //cout<<"4"+seq<<endl;
    }
    if(argsMap.count("windowlength")>0){
        string numVal = argsMap.find("windowlength")->second;
        int val_st = numVal.find('-');
        int winSz = atoi(numVal.substr(0,val_st).c_str());
        int qual = atoi(numVal.substr(val_st+1).c_str());
        WindowTrim(seq,qualSeq,winSz,qual);
    }
    if(argsMap.count("minlen")>0){
        string numVal = argsMap.find("minlen")->second;
        int len = atoi(numVal.c_str());
        if(seq.length()<len){
            doSave = false;
        }
        if(doSave){
            outfile << readAr[0]+"\n";
            outfile << seq+"\n";
            outfile << readAr[2]+"\n";
            outfile << qualSeq+"\n";
        }
    }
}

void ParseFq::QualityTrim(string &seq,string &qualSeq,int pos,int qual){
    std::vector<int> qualVect = ConvertCharToInt(qualSeq);
    
    if(pos == 1){
        int maxPosToCheck = seq.length()*1/3;
        int lowestPos = 0;
        for(auto i=0;i<maxPosToCheck && i<qualVect.size();i++){
            int curQual = qualVect.at(i);
            if(curQual<qual){
                lowestPos = i;
            }
        }
        seq = seq.substr(lowestPos);
        qualSeq = qualSeq.substr(lowestPos);
    }else if(pos == 2){
        int maxPosToCheck = seq.length()-(seq.length()*1/3);
        int lowestPos = 0;
        for(auto i=seq.length();i<=maxPosToCheck;i--){
            int curQual = qualVect.at(i);
            if(curQual<qual){
                lowestPos = i;
            }
        }
        if(lowestPos == 0){
            seq = seq;
            qualSeq = qualSeq;
        }else if(lowestPos>=maxPosToCheck){
            seq = seq.substr(0,lowestPos);
            qualSeq = qualSeq.substr(0,lowestPos);
        }
    }
}

void ParseFq::WindowTrim(string &seq,string &qualSeq,int winSize,int qual){
    std::vector<int> qualVect = ConvertCharToInt(qualSeq);
    std::vector<int> winQVect;
    for(int i=0;i<qualVect.size();i++){
        int win_Avg=0,num_bases = 0;
        for(int j=i;j<(j+winSize) && j<qualVect.size();j++){
            win_Avg=win_Avg+(qualVect[j]);
            num_bases++;
        }
        int win_res = (win_Avg/num_bases);
        if(win_res>qual){
            winQVect.push_back(i);
        }
    }
    
    vector< vector<int> > res;
    if(winQVect.size()<qualVect.size()){
        split(winQVect,res);
        
        compar comp;
        std::sort(res.begin(),res.end(),comp);
        
        int vId = 0;
        int st=0,stop = 0;
        for(vector<int> i:res){
            if(i.size()>=1){
                if(vId == 0){
                    st = i.front();
                    stop = i.back();
                }
                vId++;
            }
        }
        if(st>0 && stop<seq.length()){
            seq = seq.substr(st,stop);
            qualSeq = qualSeq.substr(st,stop);
        }
    }
}

std::vector< int > ParseFq:: ConvertCharToInt(string qualSeq){
    std::vector<int> qualVect;
    for(auto c :qualSeq){
        qualVect.push_back(((int)c)-33); //qualMap.find(c)->second
    }
    return qualVect;
}

void ParseFq::split(const std::vector<int> &origin, vector< vector<int> > &result)
{
    result.clear();
    if(origin.empty()) return;
    
    result.resize(1);
    result[0].push_back(origin[0]);
    
    for(size_t i = 1; i < origin.size(); ++i)
    {
        if(origin[i] != origin[i-1] + 1) result.push_back(vector<int>());
        result.back().push_back(origin[i]);
    }
}

void ParseFq::CheckFileReads(string file_to_load,string file_to_check){
    
    ifstream file(file_to_load);//("E:\\Data\\NGS-16\\DZA02_S2_L001_R1_001.fastq");
    string str;
    int lineCounter = 0;
    int readCounter = 0;
    std::map<std::string, std::string> m;
    vector<string>parts;
    while (std::getline(file, str))
    {
        if(lineCounter==4){
            lineCounter=0;
            if(lineCounter == 0){
                parts = splitStr(' ',0,str);
                m.emplace(parts[0],"");
            }
            lineCounter++;
        }
        else{
            if(lineCounter == 0){
                parts = splitStr(' ',0,str);
                m.emplace(parts[0],"");
            }
            lineCounter++;
        }
        parts.clear();
    }
    file.close();
    
    std::ofstream outfile;
    string out_File_Name = file_to_check;
    const size_t last_slash_idx = out_File_Name.find_last_of("\\/");
    if (std::string::npos != last_slash_idx)
    {
        out_File_Name.erase(0, last_slash_idx + 1+4);
    }
    
    out_File_Name = argsMap.find("out")->second+"trimmed_"+out_File_Name;
    //cout<<out_File_Name<<endl;
    outfile.open(out_File_Name, std::ios_base::app);
    
    ifstream fileCheck(file_to_check);
    lineCounter = 0;
    std::array<string,4> readAr;
    readCounter = 0;
    while (std::getline(fileCheck, str))
    {
        if(lineCounter==4){
            lineCounter=0;
            readCounter++;
            parts = splitStr(' ',0,readAr[0]);
            bool read_present =  m.count(parts[0])>0;
            if(read_present){
                outfile << readAr[0]+"\n";
                outfile << readAr[1]+"\n";
                outfile << readAr[2]+"\n";
                outfile << readAr[3]+"\n";
            }
            readAr.empty();
            readAr[lineCounter] = str;
            lineCounter++;
        }
        else{
            readAr[lineCounter] = str;
            lineCounter++;
        }
        parts.clear();
    }
    fileCheck.close();
    outfile.close();
}


vector<string> ParseFq::splitStr(char delim, int rep,string data) {
    vector<string> flds;
    string work = data;
    string buf = "";
    int i = 0;
    while (i < work.length()) {
        if (work[i] != delim)
            buf += work[i];
        else if (rep == 1) {
            flds.push_back(buf);
            buf = "";
        } else if (buf.length() > 0) {
            flds.push_back(buf);
            buf = "";
        }
        i++;
    }
    if (!buf.empty())
        flds.push_back(buf);
    return flds;
}

std::map<std::string, std::string> ParseFq::getOutFileMap(){
    return outFileNameMap;
}
