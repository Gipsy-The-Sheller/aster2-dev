#include <future>

//waster-lite for menory-save waster analysis

#define DEBUG

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <vector>
#include <deque>
#include <array>
#include <cmath>
#include <bitset>
#include <new>

#include "sequence_utilities.hpp"

using std::string;
using std::tuple;
using std::vector;
using std::array;


constexpr int kSize = 21;
constexpr int kFlankSize = (int)(kSize-1) / 2;
constexpr unsigned long long B_SIZE = 1ULL<<29;
constexpr unsigned long long BT_SIZE = B_SIZE / 8;
constexpr unsigned long long BUFFER_SIZE = 1024 * 1024 * 64;
constexpr char EMPTY = 0b00;
constexpr int ALIGNMENT = 1024;

using BITSET = std::bitset<B_SIZE>;
using EsTablePtr = BITSET* (*)[17];

// char mem[12ULL<<30];
char* mem = nullptr;
char* membitarray[16][3];
constexpr size_t MEM_SIZE = 12ULL << 30;

typedef long long ll;

static_assert(kFlankSize * 2 <= 64, "Sequence too long for 64-bit hash");

struct HashSquare;
struct BinSeq;

struct BinSeq{
    const static int K = 10, SHIFT = 2 * (K-1);
    const static int MASK = (1 << SHIFT) - 1;
    
    std::bitset<1ULL << 30> seq;

    int n=0, m=0;
    ll count = 0, seqlength;

    template<typename T>
    char TwoBit(T address){
        // return ((this->seq[address * 2] ? 1 : 0) << 2) | (this->seq[address * 2 + 1] ? 1 : 0);
        return (seq[address * 2]<<1 | seq[address * 2 + 1]);
    }

    bool Yieldable(){
        return count > seqlength ? false : true;
    }

    tuple<int, int> Yield(){
        // seq <<= 2;
        n = ((n & MASK) << 2) | TwoBit(count);
        m = (m >> 2) | ((3u ^ TwoBit(count + K + 1)) << SHIFT);
        
        count ++;
        return n > m ? std::make_tuple(n, m) : std::make_tuple(m, n);
    }

    BinSeq(const string& sequence) {
        seqlength = sequence.size();
        for (int i=0; i < sequence.size(); i++){

            switch (std::toupper(sequence[i])){
                case 'A':{
                    seq[2 * i] = 0; seq[2 * i + 1] = 0; break;
                } case 'C': {
                    seq[2 * i] = 0; seq[2 * i + 1] = 1; break;
                } case 'G':{
                    seq[2 * i] = 1; seq[2 * i + 1] = 0; break;
                } case 'T':{
                    seq[2 * i] = 1; seq[2 * i + 1] = 1; break;
                }
            }
        }
        for (char i = 0 ; i < K - 1; i++){
            Yield();
        }
    }
};

struct HashSquare {
    ll n, m, c;
    ll k, v, t, b, r;
    bool reverse, discard;
    ll AdjustedFlankHash, AdjustedCoreHash; 

    template <typename T>
    T RevCompHash(T x) {
        constexpr int effective_bits = 2 * kFlankSize;  // 20
        constexpr int pairs = effective_bits / 2;        // 10
        T result = 0;
        for (int i = 0; i < pairs; ++i) {
            T pair = (x >> (effective_bits - 2 - 2*i)) & 0x03;
            T rev = 3 - pair;
            result |= (rev << (effective_bits - 2 - 2*i));
        }
        return result;
    }

    HashSquare(ll n_, ll m_, ll c_) : n(n_), m(m_), c(c_) {
        discard = (n == m);
        if (!discard){
            reverse = (n > m);
            AdjustedFlankHash = !reverse ? (n >> kFlankSize) | m : (RevCompHash<ll>(m) >> kFlankSize) | RevCompHash<ll>(n);
            AdjustedCoreHash = !reverse ? c : RevCompHash<ll>(c);   
            //std::tie(m, n) = !reverse ? std::make_tuple(n, m) : std::make_tuple(RevCompHash<ll>(m), RevCompHash<ll>(n));
            if (n < m) std::tie(m, n) = std::make_tuple(n, m);
            k = n * (n - 1) / 2 + m;
            v = k / 17;
            t = k % 17;
            b = v / 65;
            r = v % 65;         
        }

    }
};

std::pair<ll, ll> recover_n_m(ll t, ll b, ll r) {
    ll k = (b * 65 + r) * 17 + t;
    ll lo = 0, hi = (1LL << 20);
    while (lo < hi) {
        ll mid = (lo + hi + 1) / 2;
        if (mid * (mid - 1) / 2 <= k) lo = mid;
        else hi = mid - 1;
    }
    ll n = lo;
    ll m = k - n * (n - 1) / 2;
    return {n, m};
}

void FilterInputWorker(char **FilterTable, string fileName, int fileorder, char* fileBuffer){
    // std::cerr << "thread: " << fileorder << std::endl;
    std::cerr << std::format("[Task {}] - thread {} - initialized.\n", fileName, fileorder);
    SeqParser* seqfile = new SeqParser(fileName, fileBuffer);
    while (seqfile->nextSeq()) {
        BinSeq* sequence = new BinSeq(seqfile->getSeq(1)); // fix to adapt for FASTQ real files.
        ll n, m;
        while(sequence->Yieldable()){
            std::tie(n, m) = sequence->Yield();
            HashSquare hs(n, m, EMPTY);
            
            if ((hs.discard || hs.t == 16) || (hs.t != fileorder)) continue;
            FilterTable[hs.t][hs.b] = std::min((int)FilterTable[hs.t][hs.b], (int)hs.r);
        }
    }
    std::cerr << std::format("[Task {}] finished.\n", fileName);

}

void CrossStatWorker(char **FilterTable, EsTablePtr EsTable, string fileName, int fileorder, char* fileBuffer){
    std::cerr << std::format("This is task: {}\n", fileorder);
    SeqParser* seqfile = new SeqParser(fileName, fileBuffer);
    while (seqfile->nextSeq()) {
        BinSeq* sequence = new BinSeq(seqfile->getSeq(1));
        ll n, m;
        while(sequence->Yieldable()){
            std::tie(n, m) = sequence->Yield();
            HashSquare hs(n, m, EMPTY);
            if (hs.discard || hs.t == 16) continue;
            if (hs.r == FilterTable[hs.t][hs.b] && hs.r < 63){
                //if (hs.t - (hs.t%4) <= fileorder && hs.t + 4 - (hs.t%4) > fileorder){
                    if (EsTable[fileorder % 4][hs.t]){
                        (*EsTable[fileorder % 4][hs.t])[hs.b] = 1;
                    }                    
                //}
            }
        }
    }
    std::cerr << std::format("[task {}] finished.\n", fileorder);
}

void StatDepoWorker(char **FilterTable, EsTablePtr EsTable, int fileorder){
    std::cerr << std::format("This is task: {}\n", fileorder);
    for(ll b=0; b<B_SIZE; b++){
        if (FilterTable[fileorder][b] != 63){
            int res = 0;
            for (int i=0; i<4; i++){
                res += EsTable[i][fileorder] ? (int)((*EsTable[i][fileorder])[b]) : 0;
            }
            FilterTable[fileorder][b] = (res << 6) | FilterTable[fileorder][b];            
        }

    }
    std::cerr << std::format("\n[task {}]: finished\n", fileorder);
}

// struct WorkFlowc {
//     WorkFlow {

//     }


// }

int main(int argc, char** argv){
    mem = static_cast<char*>(::operator new(MEM_SIZE, std::align_val_t(ALIGNMENT)));
    // TO DELETE: ::operator delete(mem, std::align_val_t(ALIGNMENT));
    if(!mem){
        throw "Failed to allocate memory!";
        exit(-1);
    }

    
    char *memstart = mem, *memEst = &mem[4LL << 30], *memFlt = memstart;

    char *bufferStart = mem + 16 * B_SIZE + 48 * BT_SIZE; // 1024 * 1024 * 64 = 64MBytes; 

    char *FilterTable[17];
    for (int i = 0; i < 16; i++){
        FilterTable[i] = mem + i * B_SIZE;
    }
    FilterTable[16] = nullptr;
    std::fill_n (FilterTable[0], 16 * (B_SIZE), 63);


    if(true){

        std::vector<std::thread> threads;
        for (int i=1; i<17; i++){
            std::string fileName = std::format("simc{}.fa", i);
            threads.emplace_back(FilterInputWorker, 
                FilterTable, std::move(fileName), i-1, bufferStart + (i-1) * BUFFER_SIZE);
        }

        for (auto& t : threads) {
            t.join();
        }
    }

    char *EsStart = mem + 16 * B_SIZE;
    BITSET *EsTable[16][3] = {};
    for (int i = 0; i < 16; i++){
        for (int j = 0; j < 3; j++){
            EsTable[i][j] = new (EsStart + (i * 3 + j) * BT_SIZE) BITSET();
        }
    }

    if(true){
        std::cerr << "Estimating...\n";

        using BITSET = std::bitset<B_SIZE>;
        // BITSET *EsTable[4][4] = {};

        BITSET *blockEsTable[4][17] = {};
        int blockEsTableCnt[17] = {};
        
        int istart = 0; // 4, 8, 12
        for (int istart = 0; istart < 16; istart += 4) {
            for (int i = istart; i < istart + 4; ++i) { // i -> file
                for (int j = istart; j < istart + 4; ++j) { // j -> t
                    if (i != j){
                        blockEsTable[i - istart][j] = EsTable[j][blockEsTableCnt[j]];
                        blockEsTableCnt[j]++;
                    }
                }
            }            
        }
        

        {
            
            std::cerr << "Adding working threads: 16 CrossStatWorkers\n";
            std::vector <std::thread> CSWthreads;
            for (int i=0; i<16; i++){
                CSWthreads.emplace_back(CrossStatWorker, FilterTable, blockEsTable, std::format("simc{}.fa", i+1), i, bufferStart + i*BUFFER_SIZE);
            }

            for (auto &i : CSWthreads) {
                if (i.joinable()){
                    i.join();
                }
            }            
        }

        {
            std::cerr << "Adding working threads: 4 StatDepoWorkers\n";
            std::vector <std::thread> SDWthreads;
            for (int i=0; i<16; i++){
                SDWthreads.emplace_back(StatDepoWorker, FilterTable, blockEsTable, i);
            }

            for (auto &i : SDWthreads) {
                if (i.joinable()){
                    i.join();
                }
            }            
        }

        std::cerr << "Validating results on file 0\n";

        {
            int temp_res[4] = {};
            for (ll i=0; i<B_SIZE; i++) {
                //char: 127 = 0b111111
                //char: 63  = 0b011111
                auto count = int((FilterTable[0][i] >> 6) & 0x03);
                temp_res[count] ++ ;
            }
            cerr << temp_res[0] << ' ' << temp_res[1] << ' ' << temp_res[2] << ' ' << temp_res[3] << endl;
        }

        

    }

    std::cerr << "Finished";

    return 0;
}