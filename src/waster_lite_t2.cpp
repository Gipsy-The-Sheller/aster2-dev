//waster-lite for menory-save waster analysis

#define DEBUG

#include <cstdio>
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

#include "sequence_utilities.hpp"

using std::string;
using std::tuple;
using std::vector;
using std::array;

const int kSize = 21;
const int kFlankSize = (int)(kSize-1) / 2;
const unsigned long long B_SIZE = 1ULL<<29;
const char EMPTY = 0b00;

typedef long long ll;

static_assert(kFlankSize * 2 <= 64, "Sequence too long for 64-bit hash");

struct HashSquare;
struct Filter;
struct BinSeq;

ll SeqHash(const string&);
string Hash2Seq(const ll);
string RevComp(string);
template<typename T> 
T RevCompHash(T);
inline HashSquare ProcessSingleMer(string);

struct BinSeq{
    const static int K = 10, SHIFT = 2 * (K-1);
    const static int MASK = (1 << SHIFT) - 1;
    
    std::bitset<1ULL << 30> seq;

    int n=0, m=0;
    ll count = 0, seqlength;

    template<typename T>
    char TwoBit(T address){
        // return ((this->seq[address * 2] ? 1 : 0) << 2) | (this->seq[address * 2 + 1] ? 1 : 0);
        return 2 * seq[address * 2] + seq[address * 2 + 1];
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

            switch (sequence[i]){
                case 'a':{
                    seq[2 * i] = 0; seq[2 * i + 1] = 0; break;
                } case 'c': {
                    seq[2 * i] = 0; seq[2 * i + 1] = 1; break;
                } case 'g':{
                    seq[2 * i] = 1; seq[2 * i + 1] = 0; break;
                } case 't':{
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

inline HashSquare ProcessSingleMer(string Mer){
    if (Mer.size() != kSize){
        throw std::invalid_argument("Invalid k-mer to process!");
    }
    auto nMer = Mer.substr(0, kFlankSize);
    char cNuc = Mer[kFlankSize];
    auto mMer = Mer.substr(kFlankSize+1, kFlankSize);

    int c;
    ll n, m;

    switch(std::toupper(cNuc)){
        case 'A': c = 0; break;
        case 'C': c = 1; break;
        case 'G': c = 2; break;
        case 'T': c = 3; break;
    }

    n = SeqHash(nMer);
    m = SeqHash(RevComp(mMer));

    HashSquare hs(n, m, c);
    return hs;
}


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

ll SeqHash(const string& seq) {
    if (seq.size() != kFlankSize){
        throw std::invalid_argument("Invalid sequence to build hash!");
    }

    ll hash = 0;
    for (char nuc : seq) {
        ll val;
        switch (std::toupper(nuc)) {
            case 'A': val = 0; break;
            case 'C': val = 1; break;
            case 'G': val = 2; break;
            case 'T': val = 3; break;
            default:
                throw std::invalid_argument(std::format("Invalid nucleotide character: {}!", int(nuc)).c_str());
        }
        hash = (hash << 2) | val;
    }
    return hash;
}

string Hash2RCSeq(const ll n) {
    string seq;
    seq.reserve(kFlankSize);

    ll tmp = n;
    for (int i = 0; i < kFlankSize; ++i) {
        char nuc;
        switch (tmp & 3) {
            case 0: nuc = 'T'; break; // A -> T
            case 1: nuc = 'G'; break; // C -> G
            case 2: nuc = 'C'; break; // G -> C
            case 3: nuc = 'A'; break; // T -> A
        }
        seq.push_back(nuc);
        tmp >>= 2;
    }
    return seq;
}

string RevComp(string Seq){
    std::reverse(Seq.begin(), Seq.end());
    string RCSeq = "";
    for(auto fNuc : Seq){
        char rNuc;
        switch(std::toupper(fNuc)){
            case 'A': rNuc = 'T'; break;
            case 'T': rNuc = 'A'; break;
            case 'C': rNuc = 'G'; break;
            case 'G': rNuc = 'C'; break;
            default: throw std::invalid_argument("Invalid nucleotide character!");
        }
        RCSeq += rNuc;
    }
    return RCSeq;
}

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

template <typename Table>
void FilterInputWorker(Table &FilterTable, string fileName, int fileorder){
    std::cerr << "thread: " << fileorder << std::endl;
    SeqParser* seqfile = new SeqParser(fileName);
    while (seqfile->nextSeq()) {
        BinSeq* sequence = new BinSeq(seqfile->getSeq(1)); // fix to adapt for FASTQ real files.
        ll n, m;
        while(sequence->Yieldable()){
            std::tie(n, m) = sequence->Yield();
            HashSquare hs(n, m, EMPTY);
            
            if ((hs.discard || hs.t == 16) || (hs.t != fileorder)) continue;
            FilterTable(hs.t, hs.b) = std::min(int(FilterTable(hs.t, hs.b)), int(hs.r));
        }
    }
    std::cerr << std::format("[Task {}]", fileName) << "finished.\n";

}

template <typename FTable, typename ETable>
void CrossStatWorker(FTable &FilterTable, ETable &EsTable, string fileName, int fileorder){
    std::cerr << std::format("This is task: {}\n", fileorder);
    SeqParser* seqfile = new SeqParser(fileName);
    while (seqfile->nextSeq()) {
        BinSeq* sequence = new BinSeq(seqfile->getSeq(1));
        ll n, m;
        while(sequence->Yieldable()){
            std::tie(n, m) = sequence->Yield();
            HashSquare hs(n, m, EMPTY);
            
            if (hs.discard || hs.t == 16) continue;
            // if (hs.t >=0 && hs.t <4 && hs.t != fileorder && hs.r == FilterTable(fileorder, hs.b) && hs.r < 63){
            //     (*EsTable[fileorder][hs.t])[hs.b] = 1;
            // }
            if (hs.t >=0 && hs.t <4 && hs.t != fileorder && hs.r == FilterTable(hs.t, hs.b) && hs.r < 63){
                (*EsTable[fileorder][hs.t])[hs.b] = 1;
            }
        }
    }
    std::cerr << std::format("[task {}]: finished\n", fileorder);
}

template <typename FTable, typename ETable>
void StatDepoWorker(FTable &FilterTable, ETable &EsTable, int fileorder, ll (&bottle)[4]){
    std::cerr << std::format("This is task: {}\n", fileorder);
    for(ll b=0; b<B_SIZE; b++){
        //128 = 2*2*2*2*2*2*2 = 10000000 so 127 = 1111111 64=1000000 63=111111
        int res = ((*EsTable[0][fileorder])[b] + (*EsTable[1][fileorder])[b] + (*EsTable[2][fileorder])[b] + (*EsTable[3][fileorder])[b]);
        if (res) cerr << res << endl;
        bottle[res] += 1;
        FilterTable(fileorder, b) = (res << 6) | FilterTable(fileorder, b);
        if (b % 1000000 == 0){
            std::cerr << std::format("\r[task {}]: main loop time {}", fileorder, b);
        }
    }
    std::cerr << std::format("\n[task {}]: finished\n", fileorder);
}

int main(int argc, char** argv){
    std::cerr<<"Initializing\n";

    char* FTdata = new char[16 * (B_SIZE)];

    std::fill_n (FTdata, 16 * (B_SIZE), 63);

    auto FilterTable = [&](size_t i, size_t j) -> char& {
        return FTdata[i * (B_SIZE) + j];
    };

    if(true){
        std::cerr << "Testing reads\n";

        std::vector<std::thread> threads;
        for (int i=1; i<17; i++){
            std::cerr << "Initializing thread " << i << std::endl;
            std::string fileName = std::format("simc{}.fa", i);
            threads.emplace_back(FilterInputWorker<decltype(FilterTable)>, 
                std::ref(FilterTable), std::move(fileName), i);
        }

        for (auto& t : threads) {
            t.join();
        }
    }

    if(true){
        std::cerr << "Estimating...\n";
        std::vector<std::vector<std::unique_ptr<std::bitset<B_SIZE>>>> EsTable;
        //4*16
        std::cerr << "Initialized ES Table\n";
        EsTable.reserve(4);
        for (int i = 0; i < 4; ++i) {
            EsTable.emplace_back(4);
        }

        for (auto& innerVec : EsTable) {
            for (auto& ptr : innerVec) {
                ptr = std::make_unique<std::bitset<B_SIZE>>();
            }
        }
        
        std::cerr << "Adding working threads: 4 CrossStatWorkers\n";

        std::thread CSW1(CrossStatWorker<decltype(FilterTable), decltype(EsTable)>, std::ref(FilterTable), std::ref(EsTable), "simc1.fa", 0);
        std::thread CSW2(CrossStatWorker<decltype(FilterTable), decltype(EsTable)>, std::ref(FilterTable), std::ref(EsTable), "simc2.fa", 1);
        std::thread CSW3(CrossStatWorker<decltype(FilterTable), decltype(EsTable)>, std::ref(FilterTable), std::ref(EsTable), "simc3.fa", 2);
        std::thread CSW4(CrossStatWorker<decltype(FilterTable), decltype(EsTable)>, std::ref(FilterTable), std::ref(EsTable), "simc4.fa", 3);
        CSW1.join(); CSW2.join(); CSW3.join(); CSW4.join();

        std::cerr << "Adding working threads: 4 StatDepoWorkers\n";

        ll res1[4] = {0,0,0,0}, res2[4] = {0,0,0,0}, res3[4] = {0,0,0,0}, res4[4] = {0,0,0,0};
        std::thread SDW1(StatDepoWorker<decltype(FilterTable), decltype(EsTable)>, std::ref(FilterTable), std::ref(EsTable), 0, std::ref(res1));
        std::thread SDW2(StatDepoWorker<decltype(FilterTable), decltype(EsTable)>, std::ref(FilterTable), std::ref(EsTable), 1, std::ref(res2));
        std::thread SDW3(StatDepoWorker<decltype(FilterTable), decltype(EsTable)>, std::ref(FilterTable), std::ref(EsTable), 2, std::ref(res3));
        std::thread SDW4(StatDepoWorker<decltype(FilterTable), decltype(EsTable)>, std::ref(FilterTable), std::ref(EsTable), 3, std::ref(res4));
        SDW1.join(); SDW2.join(); SDW3.join(); SDW4.join();

        std::cerr << res1[0] << ' ' << res1[1] << ' ' << res1[2] << ' ' << res1[3] << std::endl;
        std::cerr << res2[0] << ' ' << res2[1] << ' ' << res2[2] << ' ' << res2[3] << std::endl;
        std::cerr << res3[0] << ' ' << res3[1] << ' ' << res3[2] << ' ' << res3[3] << std::endl;
        std::cerr << res4[0] << ' ' << res4[1] << ' ' << res4[2] << ' ' << res4[3] << std::endl;

    }

    std::cerr << "Finished";

    return 0;
}