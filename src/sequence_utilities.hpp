#ifndef SEQ_UTILS
#define SEQ_UTILS

#include<queue>
#include<tuple>
#include<algorithm>
#include<random>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<unordered_set>

#include "common.hpp"



using namespace std;

namespace SeqUtils{
    string QUALITY2ASCII = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

    string fastaFormatName(const string &name){
        string res;
        for (char c: name){
            if (c != '>' && c != ' ' && c != '\t') res += c;
        }
        return res;
    }

    string fastaFormatRead(const string &name){
        string res;
        for (char c: name){
            if (('a' <= c && 'z' >= c) || ('A' <= c && 'Z' >= c) || c == '-') res += c;
        }
        return res;
    }
}

struct SeqParser{
    bool isFasta = false;
    string seq, quality, L1;
    ifstream fin;

    common::LogInfo LOG;

    const long long STEP_SIZE = 1 << 27;
    long long cnt = 0, lenCnt = 0, threshold = STEP_SIZE;
    
    SeqParser(string fileName, string fileFormat = "auto", int verbose = common::LogInfo::DEFAULT_VERBOSE) : LOG(verbose){
        ifstream ftemp(fileName);
        string temp;
        if (!getline(ftemp, temp)){
            cerr << "Input file path '" << fileName << "' seems empty!\n";
            exit(0);
        }
        if (temp.length() == 0){
            cerr << "Input file path '" << fileName << "' seems empty!\n";
            exit(0);
        }
        if (fileFormat == "fastq" || (fileFormat == "auto" && seemsFastq(temp))) initFastq(fileName);
        else if (fileFormat == "fasta" || (fileFormat == "auto" && seemsFasta(temp))) initFasta(fileName);
        else {
            cerr << "Input file '" << fileName << "' bad format! Only FASTA and FASTQ are supported! (Did you forget unzipping it?)\n";
            exit(0);
        }
    }

    bool nextSeq(){
        if (isFasta) return nextFastaSeq();
        else return nextFastqSeq();
    }

    string getSeq(size_t threshold = 0){
        if (isFasta) return seq;
        if (threshold >= SeqUtils::QUALITY2ASCII.length()) threshold = SeqUtils::QUALITY2ASCII.length() - 1;
        string res;
        for (size_t i = 0; i < seq.size(); i++){
            if (quality[i] >= SeqUtils::QUALITY2ASCII[threshold]) res += seq[i];
            else res += 'N';
        }
        return res;
    }

    static string STRING_REVERSE_COMPLEMENT(const string &seq){
        string res;
        for (int i = seq.size() - 1; i >= 0; i--){
            switch (seq[i]) {
                case 'A': res += 'T'; break;
                case 'C': res += 'G'; break;
                case 'G': res += 'C'; break;
                case 'T': res += 'A'; break;
                default: res += 'N';
            }
        }
        return res;
    }

    static string formatSeq(const string raw){
        string res;
        for (char c: raw){
            switch (c) {
                case 'A': case 'a': res += 'A'; break;
                case 'C': case 'c': res += 'C'; break;
                case 'G': case 'g': res += 'G'; break;
                case 'T': case 't': case 'U': case 'u': res += 'T'; break;
                case ' ': case '\t': case '\r': case '\n': break; 
                default: res += 'N';
            }
        }
        return res;
    }

private:
    static bool seemsFastq(const string &line){
        return line[0] == '@';
    }

    static bool seemsFasta(const string &line){
        return line[0] == '>';
    }

    void initFastq(const string &fileName){
        isFasta = false;
        fin.open(fileName);
    }

    void initFasta(const string &fileName){
        isFasta = true;
        fin.open(fileName);
        if (!getline(fin, L1)){
            cerr << "Input file path '" << fileName << "' seems empty!\n";
            exit(0);
        }
        if (L1.length() == 0) getline(fin, L1);
    }

    bool nextFastaSeq(){
        if (L1 == "") return false;

        string line;
        seq = "";
        while (getline(fin, line)){
            if (line.length() > 0 && line[0] == '>') {
                cnt++;
                lenCnt += seq.length();
                if (lenCnt >= threshold){
                    while (lenCnt >= threshold) threshold += STEP_SIZE;
                    LOG << cnt << " sequences read (" << lenCnt << " BPs); the last sequence is '" + L1 + "'" << endl;
                }

                L1 = line;
                return true;
            }
            seq += formatSeq(line);
        }
        
        cnt++;
        lenCnt += seq.length();
        if (lenCnt >= threshold){
            while (lenCnt >= threshold) threshold += STEP_SIZE;
            LOG << cnt << " sequences read (" << lenCnt << " BPs); the last sequence is '" + L1 + "'" << endl;
        }

        L1 = "";
        return true;
    }

    bool nextFastqSeq(){
        string L2, L3;
        if (!getline(fin, L1)) return false;
        if (L1.length() == 0) getline(fin, L1);
        getline(fin, L2);
        getline(fin, L3);
        getline(fin, quality);
        seq = formatSeq(L2);

        cnt++;
        lenCnt += seq.length();
        if (lenCnt >= threshold){
            while (lenCnt >= threshold) threshold += STEP_SIZE;
            LOG << cnt << " sequences read (" << lenCnt << " BPs); the last sequence is '" + L1 + "'" << endl;
        }
        return true;
    }
};

#endif