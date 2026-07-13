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

using std::string;
using std::tuple;
using std::vector;
using std::array;

//kSize: odd for WASTER.
#ifdef DEBUG

    // std::vector<std::vector<unsigned char>> FilterTable(16ULL * (1ULL << 29), 64);
    

#else



#endif

const int kSize = 21;
const int kFlankSize = (int)(kSize-1) / 2;

typedef long long ll;

static_assert(kFlankSize * 2 <= 64, "Sequence too long for 64-bit hash");

struct HashSquare;
struct Filter;

ll SeqHash(const string&);
string Hash2Seq(const ll);
string RevComp(string);
template<typename T> 
T RevCompHash(T);
inline HashSquare ProcessSingleMer(string);

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

struct Filter{
    
    Filter(){
        //input by ifstream/ofstream
        //yield k-mer by coroutine
        if (true){
            string Mer = "ATCG";
            //split into n m c
        }
    }
};

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

void RandomizedFilter(){
    //ifstream
    ll Filtered[16][1<<29];

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

int main(int argc, char** argv){
    std::cerr<<"Initializing\n";

    char* FTdata = new char[16 * (1ULL<<29)];

    std::fill_n (FTdata, 16 * (1ULL<<29), 64);

    auto FilterTable = [&](size_t i, size_t j) -> char& {
        return FTdata[i * (1ULL<<29) + j];
    };
    // if(true){
    //     //TEST 1: 21-mer parsing and 
    //     vector<string> MerVec = {
    //         string("AAAAAAAAAAGTTTTTTTTTT"),
    //         string("TAAAAAAAAAGTTTTTTTTTT"),
    //         string("ATTTTTTTTTGTTTTTTTTTT")
    //     };

    //     for (auto Mer : MerVec){
    //         auto hs = ProcessSingleMer(Mer);
    //         std::cerr << "This is mer " << '"' << Mer << '"' << std::endl;
    //         std::cerr << "discard? - " << hs.discard << std::endl;
    //         std::cerr << "reverse? - " << hs.reverse << std::endl;   
            
    //         //try restore the kmer from m n c
    //         std::cerr << "kmer: " << RevComp(Hash2RCSeq(hs.n)) << ' ' << RevComp(Hash2RCSeq(hs.c)) << ' ' << Hash2RCSeq(hs.m) << std::endl;
    //     }
    // }

    if(true){     
        // TEST 2: core k-mers callback
        std::cerr << "Initialization started.\n";

        std::cerr << "Finished.\n";

        for (int i = 0; i < 16; i ++){
            std::cerr<<"\nStart to parse: sim" << i+1 << std::endl;
            auto valid = freopen(std::format("./sim{}.fa", i+1).c_str(), "r", stdin);
            if (valid == nullptr) std::cerr << "Failed to open!\n";
            std::cin.clear(); 


            std::string header;
            std::getline(std::cin, header); // discard FASTA header
            std::cerr << "Header: " << header << std::endl;
            
            std::deque<char> window;
            char ch;
            int count = 0;
            
            // FilterTable(0,0) = char(0);
            while (std::cin.get(ch)) {
                if (ch != 'a' && ch != 'c' && ch != 't' && ch != 'g') {
                    // 打印所有被跳过的字符（包括换行符等）
                    // std::cerr << "Skipping char: " << int((unsigned char)ch) << " ('" << ch << "')\n";
                    continue;
                }
                window.push_back(ch);
                if (count % (1000 * 1000 * 10) == 0) {
                    std::cerr<<'\r'<<"Processing: "<<count << '\n';
                }
                if (window.size() == kSize) {
                    std::string kmer(window.begin(), window.end());

                    
                    auto hs = ProcessSingleMer(kmer);
                    if (hs.discard || (hs.t == 16 || hs.b == 64) || (hs.t != i)) {window.pop_front();count+=1;continue;}
                    // std::cerr<<"Valid kmer: " << kmer << " - ";
                    FilterTable(hs.t, hs.b) = std::min(int(FilterTable(hs.t, hs.b)), int(hs.r));
                    // std::cerr<< FilterTable(hs.t, hs.b) << std::endl;
                    
                    window.pop_front();

                }
                count += 1;
            }

            
            fclose(stdin);
            
            std::cerr<<"\nStart to parse: sim" << 17 << std::endl;
            valid = freopen("./sim17.fa", "r", stdin);
            if (valid == nullptr) std::cerr << "Failed to open!\n";
            std::cin.clear(); 

            std::getline(std::cin, header);
            std::cerr << "Header: " << header << std::endl;

            window.clear();
            count = 0;
            while (std::cin.get(ch)) {
                if (ch != 'a' && ch != 'c' && ch != 't' && ch != 'g') {
                    // 打印所有被跳过的字符（包括换行符等）
                    // std::cerr << "Skipping char: " << int((unsigned char)ch) << " ('" << ch << "')\n";
                    continue;
                }
                window.push_back(ch);
                if (count % (1000 * 1000) == 0) std::cerr<<'\r'<<"Processing: "<<count;
                if (window.size() == kSize) {
                    std::string kmer(window.begin(), window.end());
                    
                    auto hs = ProcessSingleMer(kmer);

                    if (hs.discard || (hs.t == 16 || hs.b == 64) || (hs.t != i)) {window.pop_front();count+=1;continue;}

                    FilterTable(hs.t, hs.b) = std::min(int(FilterTable(hs.t, hs.b)), int(hs.r));
                    
                    window.pop_front();
                }
                count += 1;
            }

            fclose(stdin);
        }

        fclose(stdin);

        std::cerr <<"\nfinish parsing.\n";

        std::cerr << "Previewing partial matrix\n";

        for (int i=0; i<16; i++){
            for (int j=0; j<16; j++){
                std::cerr << (int)FilterTable(i, j) << ' ';
            }
            putchar('\n');
        }

        // std::cerr << "Start recovering k-mers\n";
        
        // for(ll i = 0; i < 16; i++){
        //     for(ll j=0; j < 1ULL<<29; j++){
        //         if (int(FilterTable(i, j)) != 64){
        //             std::cerr << "Found valid k-mer:" << i << ' ' << j;
        //             auto t=i, b=j, r=ll(FilterTable(i, j));
        //             auto v=b * 65 + r;
        //             auto k=v*17+t;
        //             auto nm = recover_n_m(t, b, r);
        //             auto n = nm.first;
        //             auto m = nm.second;
        //             std::cerr << RevComp(Hash2RCSeq(n)) << 'N' << Hash2RCSeq(m) << std::endl;
        //         }
        //     }
        // }

        std::cerr << "Start stating k-mers\n";

        auto valid = freopen("./sim17.fa", "r", stdin);
        if (valid == nullptr) std::cerr << "Failed to open!\n";
        std::cin.clear(); 

        string header;
        std::getline(std::cin, header);
        std::cerr << "Header: " << header << std::endl;

        std::deque<char> window;
        char ch;

        window.clear();
        ll count = 0;
        ll kmer_count = 0;
        ll hits = 0;
        ll cnt[16] = {};
        while (std::cin.get(ch)) {
            if (ch != 'a' && ch != 'c' && ch != 't' && ch != 'g') {
                // 打印所有被跳过的字符（包括换行符等）
                // std::cerr << "Skipping char: " << int((unsigned char)ch) << " ('" << ch << "')\n";
                continue;
            }
            window.push_back(ch);
            if (count % (1000 * 1000) == 0) std::cerr<<'\r'<<"Processing: "<<count;
            if (window.size() == kSize) {
                std::string kmer(window.begin(), window.end());
                kmer_count ++;
                
                auto hs = ProcessSingleMer(kmer);

                if (hs.discard || (hs.t == 16 || hs.b == 64)) {window.pop_front();count+=1;continue;}

                if (FilterTable(hs.t, hs.b) == hs.r && hs.r != 64) {hits++;cnt[hs.t]++;}
                
                window.pop_front();
            }
            count += 1;
        }
        std::cerr << "total kmers: " << kmer_count << '\n';
        std::cerr << "       hits: " << hits << '\n';

        for (auto i: cnt) std::cerr << i << ' ';
        std::cerr << '\n';
    } 
    return 0;
}