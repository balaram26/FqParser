#include <iostream>
#include <sstream>
#include <thread>
#include <cstdio>


#include "ParseFq.h"
using namespace std;

std::map<std::string, std::string> mappify2(std::string const& s)
{
    std::map<std::string, std::string> m;
    
    std::string::size_type key_pos = 0;
    std::string::size_type key_end;
    std::string::size_type val_pos;
    std::string::size_type val_end;
    
    while((key_end = s.find('=', key_pos)) != std::string::npos)
    {
        if((val_pos = s.find_first_not_of("= ", key_end)) == std::string::npos)
            break;
        
        val_end = s.find('\t', val_pos);
        m.emplace(s.substr(key_pos, key_end - key_pos), s.substr(val_pos, val_end - val_pos));
        
        key_pos = val_end;
        if(key_pos != std::string::npos)
            ++key_pos;
    }
    return m;
}

void ProcessPE(const map<string,string> &argsMap){
    for(auto const& p: argsMap)
        std::cout << '{' << p.first << " => " << p.second << '}' << '\n';
    
    ParseFq fqParser(argsMap);
    std::thread t1 (&ParseFq::ReadFq,&fqParser,argsMap.find("pe1")->second,"pe1"), t2(&ParseFq::ReadFq,&fqParser,argsMap.find("pe2")->second,"pe2");
    t1.join();
    t2.join();
    cout<<"Trimming completed"<<endl;
    fqParser.CheckFileReads(fqParser.getOutFileMap().find("pe2")->second,fqParser.getOutFileMap().find("pe1")->second);
    fqParser.CheckFileReads(fqParser.getOutFileMap().find("pe1")->second,fqParser.getOutFileMap().find("pe2")->second);
    for(std::pair<string,string> p: fqParser.getOutFileMap())
        remove(p.second.c_str());
    cout<<"Checking completed"<<endl;
}

void ProcessSE(const map<string,string> &argsMap){
    for(auto const& p: argsMap)
        std::cout << '{' << p.first << " => " << p.second << '}' << '\n';
    
    ParseFq fqParser(argsMap);
    std::thread t1 (&ParseFq::ReadFq,&fqParser,argsMap.find("se")->second,"se");
    t1.join();
    cout<<"Trimming completed"<<endl;
    fqParser.CheckFileReads(fqParser.getOutFileMap().find("se")->second,fqParser.getOutFileMap().find("se")->second);
    cout<<"Checking completed"<<endl;
}

void add_default_parms(map<string,string>& argsMap){
    if(argsMap.count("lcrop")<=0){
        argsMap.emplace("lcrop","10");
    }
    if(argsMap.count("tcrop")<=0){
        argsMap.emplace("tcrop","10");
    }
    if(argsMap.count("leading")<=0){
        argsMap.emplace("leading","20");
    }
    if(argsMap.count("trailing")<=0){
        argsMap.emplace("trailing","20");
    }
    if(argsMap.count("windowlength")<=0){
        argsMap.emplace("windowlength","4-20");
    }
    if(argsMap.count("minlen")<=0){
        argsMap.emplace("minlen","36");
    }
}

int main(int argc,char *argv[])
{
    //    string fileName = "E:\\Data\\NGS-16\\DZA02_S2_L001_R1_001.fastq";
    //    string dir = "E:\\Data\\NGS-16\\";
    //    string tFNames = "pe1="+dir+"DZA02_S2_L001_R1_001.fastq pe2="+dir+"DZA02_S2_L001_R2_001.fastq";
    //    string outDir = "E:\\Data\\Results\\";
    //
    //    std::string temout="";
    //    temout = "lcrop=10 tcrop=50 leading=20 trailing=20 windowlength=4-20 minlen=36 "+tFNames+" out="+outDir;
    //    cout<<temout<<endl;
    
    std::ostringstream input_Args;
    for(int i = 1; i < argc; i++)
        input_Args << argv[i] << "\t";
    
    map<string,string> argsMap = mappify2(input_Args.str());//("lcrop=10\ttcrop=50\tleading=20\ttrailing=20\twindowlength=4-20\tminlen=36\t"+tFNames+"\tout="+outDir);
    
    if(argsMap.count("se")>=1 && argsMap.count("out")>=1){
        add_default_parms(argsMap);
        ProcessSE(argsMap);
    }else if(argsMap.count("pe1")>=1 && argsMap.count("pe2")>=1 && argsMap.count("out")>=1){
        add_default_parms(argsMap);
        ProcessPE(argsMap);
    }else{
        cout<<"Input arguments dont match!!\n"<<endl;
        cout<<"FqParser.exe pe1=filename pe2=filename out=output_directory [options]\n"<<endl;
        cout<<"FqParser.exe se=filename out=output_directory [options]\n"<<endl;
        cout<<"Options:\n"<<endl;
        string opt = "";
        opt = "lcrop=10 crop 10 bases from 3 prime end of the read.\n";
        opt += "tcrop=20 crop 10 bases from 5 prime end of the read.\n";
        opt += "leading=20 check quality of bases from 3' end of the read to the center.\n";
        opt += "trailing=20 check quality of bases from 5' end of the read to the center.\n";
        opt += "windowlength=4-20 scroll through a bases with window length of 4 and check average base quality of 20\n";
        opt += "minlen=36 Check minimum read length of 36 bases\n";
        cout<<opt<<endl;
    }
    
    return 0;
}
