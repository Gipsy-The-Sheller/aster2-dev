//waster-lite for menory-save waster analysis

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
#include <cstring>
#include <cstdint>

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
constexpr unsigned long long MAX_B = (((1ULL << 39) / 17 + 1) / 65 + 1);
constexpr unsigned long long MAX_BPRIME = MAX_B / 63 + 1;
constexpr int LMERADDIFLANK = 6;
constexpr int FLANK_BITS = 2 * LMERADDIFLANK;        // bits per flank (6-mer = 12)
constexpr int FLANK_MASK = (1 << FLANK_BITS) - 1;    // 0xFFF
constexpr int LF_SHIFT   = FLANK_BITS + 2;           // lf offset in stored largeKmer (rf 12-bit + c 2-bit = 14)

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

struct KMerInfo {
    short r_prime;
    char freq;
    bool del;
    unsigned long long consensusL, consensusR;
} ;

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

    inline void YieldInit() {
        for (char i = 0 ; i < K - 1; i++){
            Yield();
        }        
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

        YieldInit();
    }
};

struct LargeBinSeq : BinSeq {
    const static int ADDIFLANK = LMERADDIFLANK;
    const static int ADDISHIFT = (ADDIFLANK + 1) * 2;
    const static int ADDIMASK  = (1 << ADDISHIFT) - 1;
    int lf = 0, rf = 0;

    inline void YieldInit() {
        for (char i = 0 ; i < K + ADDIFLANK - 1; i++){
            YieldFlanks();
        }         
    }

    bool Yieldable () {
        return count + K + ADDIFLANK > seqlength ? false : true; 
    }

    tuple<int, int, int, int, int> YieldFlanks(){
        // return n, m for core 21-mers, and n, m for additional 5+5 flanks

        n = ((n & MASK) << 2) | TwoBit(count);
        m = (m >> 2) | ((3u ^ TwoBit(count + K + 1)) << SHIFT);
        lf = ((lf & ((1 << (2*(ADDIFLANK-1))) - 1)) << 2) | TwoBit(count - K + 1);
        // rf = ((rf & ((1 << (2*(ADDIFLANK-1))) - 1)) << 2) | TwoBit(count + ADDIFLANK + 1);
        rf = ((rf & ((1 << (2*(ADDIFLANK-1))) - 1)) << 2) | ((3u ^ TwoBit(count + K + 1 + ADDIFLANK)) << (2*(ADDIFLANK-1)));
        
        count ++;
        return n > m ? std::make_tuple(n, m, lf, rf, (int)TwoBit(count+1)) : std::make_tuple(m, n, rf, lf, (int)(3-TwoBit(count+1)));
    }

    LargeBinSeq(const string& sequence): BinSeq(sequence) { }
};

struct TrieBinSeq : BinSeq {
    inline void YieldInit() { }

    bool Yieldable() {
        return count < seqlength ? true : false; 
    }

    char SingleYield() {
        count++;
        return TwoBit(count-1);
    }

    TrieBinSeq(const string& sequence): BinSeq(sequence) { }
}

struct HashSquare {
    ll n, m, c;
    ll k, v, t, b, r;
    bool reverse, discard;
    ll AdjustedFlankHash, AdjustedCoreHash;
    ll b_prime, r_prime; 

    template <typename T>
    static T RevCompHash(T x) {
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
            b_prime = b / 63;
            r_prime = (b % 63) * 65 + r;
        }

    }
};

struct Trie {
    static const int TRIEDEPTH = kFlankSize + LMERADDIFLANK;   // 16
    static constexpr unsigned long long capacity = 8ULL << 30; // 8 GB arena

    // A node is either:
    //   - internal: children[] hold real child pointers, is_tip == 0; or
    //   - a tip:    children[0] holds a KMerInfo* (back-pointer to the source kmer),
    //               is_tip == 1 (left / fullL flank) or 2 (right / fullR flank).
    // A tip stands for the whole remaining suffix of a flank WITHOUT materialising it;
    // it is expanded on demand by extend() when a second flank walks into it.
    struct Node {
        Node *parent = nullptr;
        Node *fail = nullptr;
        Node *children[4] = {nullptr, nullptr, nullptr, nullptr};
        std::bitset<2> label = 0;       // edge base (0..3) leading into this node
        std::bitset<2> is_tip = 0;      // 0 = internal; 1 = tip (left); 2 = tip (right)

        Node(Node *p = nullptr, std::bitset<2> lbl = 0b00) : parent(p), label(lbl) {
            if (parent) parent->children[label.to_ulong()] = this;
        }
    };

    static constexpr size_t NODE_SIZE = sizeof(Node);

    /*
    strategy (lazy / compacted trie over 16-mer flanks):
    - inserting a flank walks RC(low-bit-first) down the trie; where it diverges from
      existing structure it drops a *tip* (children[0] = KMerInfo*, is_tip = direction)
      instead of allocating the rest of the path.
    - when a later flank walks INTO a tip, extend() materialises one level of the mounted
      flank (read from KMerInfo.consensusL / consensusR, which hold the full 16-mer),
      pushing the tip one level deeper; repeat until the two flanks diverge or reach
      depth 16.  ("recover the mounted kmer ... extend ... judge conflicts")
    */

    char* memory;
    size_t used = 0;
    Node* root = nullptr;

    Node* allocate_node() {
        if (used + NODE_SIZE > capacity) throw std::bad_alloc();
        Node* ptr = reinterpret_cast<Node*>(memory + used);
        used += NODE_SIZE;
        return ptr;
    }

    void mountMerPtr(Node* node, KMerInfo* merptr, short direction) {
        node->is_tip = std::bitset<2>(direction); // 0: internal node; 1: left; 2: right
        node->children[0] = reinterpret_cast<Node*>(merptr);
    }

    // turn a tip node into an internal node by materialising the next base of its
    // mounted flank; the mounted flank itself moves one level deeper as a fresh tip.
    void extend(Node* node, int depth) {
        KMerInfo* mounted = reinterpret_cast<KMerInfo*>(node->children[0]);
        short dir = (short)node->is_tip.to_ulong();                 // 1 or 2
        uint32_t M = (uint32_t)(dir == 1 ? mounted->consensusL : mounted->consensusR);
        int label = 3 - (int)((M >> (2 * depth)) & 0b11);           // RC of this depth's base
        node->children[0] = nullptr;                                // clear the merptr slot
        node->is_tip = 0;                                           // node becomes internal
        Node* child = ::new (allocate_node()) Node(node, std::bitset<2>(label));
        child->is_tip = std::bitset<2>(dir);                        // mounted flank continues as a tip
        child->children[0] = reinterpret_cast<Node*>(mounted);
    }

    void add(uint32_t mer, KMerInfo* merptr, short direction) {
        Node* cur = root;
        for (int depth = 0; depth < TRIEDEPTH; depth++) {
            int label = 3 - (int)(mer & 0b11);   // RC of current low 2 bits
            mer >>= 2;
            if (cur->is_tip.any()) extend(cur, depth);      // lazily open the tip before descending
            if (cur->children[label]) {
                cur = cur->children[label];                 // shared edge: go deeper
            } else {
                // diverges here: drop a tip for the rest of this flank and stop
                Node* child = ::new (allocate_node()) Node(cur, std::bitset<2>(label));
                mountMerPtr(child, merptr, direction);
                return;
            }
        }
        // all 16 levels matched existing edges -> identical flank already indexed; (re)mount at leaf
        mountMerPtr(cur, merptr, direction);
    }

    // // query: return the indexed KMerInfo* if `mer` is present, else nullptr.
    // KMerInfo* find(uint32_t mer) const {
    //     Node* cur = root;
    //     for (int depth = 0; depth < TRIEDEPTH; depth++) {
    //         if (cur->is_tip.any()) {
    //             // tip holds a full flank; the top `depth` bases already match the path,
    //             // so only the remaining suffix needs comparing (mer has been shifted `depth` times).
    //             KMerInfo* mounted = reinterpret_cast<KMerInfo*>(cur->children[0]);
    //             short dir = (short)cur->is_tip.to_ulong();
    //             uint32_t M = (uint32_t)(dir == 1 ? mounted->consensusL : mounted->consensusR);
    //             return ((M >> (2 * depth)) == mer) ? mounted : nullptr;
    //         }
    //         int label = 3 - (int)(mer & 0b11);
    //         mer >>= 2;
    //         cur = cur->children[label];
    //         if (!cur) return nullptr;
    //     }
    //     // descended all 16 levels on existing edges -> cur is the leaf tip
    //     return cur->is_tip.any() ? reinterpret_cast<KMerInfo*>(cur->children[0]) : nullptr;
    // }

    void init_fails() {
        std::deque<Node*> queue;
        queue.push(root);

        while(!queue.empty()) {
            auto current = queue.front();
            queue.pop();

            if (current->parent == root) {
                current->fail = root;
            } else {
                current->fail = current->parent->fail->children[current->label.to_ulong()] ? current->parent->fail->children[current->label.to_ulong()] : root;
            }

            if (!current->is_tip) {
                for (auto i : current->children) {
                    if (i) queue.push(i);
                }
            }
        }
    }

    Trie(char* mem) : memory(mem), used(0) {
        root = ::new (allocate_node()) Node(nullptr, std::bitset<2>(0));
    }

    // for 16-mers, the maximum summary of nodes is:
    // 4^16 + 4^15 + ... + 4 = 4 (4^16 - 1) / (4-1) ~ 4^16 = 2^32 = ((2^10)^3) * 2^2 ~ 3GB * 4 = 12GB ?!

    // for total kmers (128MB * 16bp), maximum: 492,131,669 nodes, each with 2 pointer
};

void BuildTrie(Trie* trie, KMerInfo (*kMerInfoTable)[16][8 << 20]) {
    for (int t = 0; t < 16; t++) {
        for (ll bp = 0; bp < (ll)MAX_BPRIME; bp++) {
            KMerInfo& info = (*kMerInfoTable)[t][bp];
            if (!info.r_prime || info.del) continue;
            auto conL = info.consensusL, conR = HashSquare.RevCompHash<int>(info.consensusR);
            trie->add(conL, &info, 1);
            trie->add(conR, &info, 2);
        }
    }
}

class MerQuery {
    // mismatch-tolerate kmer query (1 multihit)
    Trie *trie;
    string filename;
    SeqParser* seqfile;

    char CallSNP () {
        //return SNP
    }    

    void MerMatch (Trie::Node* node, short lastmismatch = -1) {
        if (node->is_tip.to_ulong()) {
            // precised match
            return;
        }

        //process current node
        
        if () {
            // failed to match
            return;
        }

        for (auto i : node->children) {

        }
    }

public:
    MerQuery (Trie *_trie, string _filename, char* buffer): trie(_trie), filename(_filename) {
        seqfile = new SeqParser(filename, buffer);
        while (seqfile->nextSeq()) {
            BinSeq* sequence = new TrieBinSeq(seqfile->getSeq(1));
            ll n, m;
            while(sequence->Yieldable()){
                auto bit = sequence->SingleYield();
            }
        }
    }
}

    // mismatch-tolerate kmer query (1 multihit)

// void BuildTrie(Trie* trie, KMerInfo (*kMerInfoTable)[16][8 << 20]) {
//     ll sites = 0;
//     for (int t = 0; t < 16; t++) {
//         for (ll bp = 0; bp < (ll)MAX_BPRIME; bp++) {
//             KMerInfo& info = (*kMerInfoTable)[t][bp];
//             if (!info.r_prime || info.del) continue;
//             uint32_t fullL = (uint32_t)info.consensusL;     // [lf6 | n_rec]
//             uint32_t fullR = (uint32_t)info.consensusR;     // [m_rec | rf6]
//             uint32_t n_rec = fullL & 0xFFFFF, lf6 = (fullL >> 20) & 0xFFF;
//             uint32_t m_rec = (fullR >> 12) & 0xFFFFF, rf6 = fullR & 0xFFF;
//             uint32_t flank5p, flank3p;                       // forward biological 5'/3' 16-mers
//             if (info.reverse) {          // n_rec=SeqHash(nMer); m_rec=SeqHash(RC(mMer)); lf6=lf; rf6=RC(rf)
//                 flank5p = (lf6 << 20) | n_rec;
//                 flank3p = (rc_n(m_rec, 10) << 12) | rc_n(rf6, 6);
//             } else {                     // m_rec=SeqHash(nMer); n_rec=SeqHash(RC(mMer)); rf6=lf; lf6=RC(rf)
//                 flank5p = (rf6 << 20) | m_rec;
//                 flank3p = (rc_n(n_rec, 10) << 12) | rc_n(lf6, 6);
//             }
//             trie->add(flank3p, &info, 1);            // 3' flank forward: allele = read[i-16]
//             trie->add(rc_n(flank5p, 16), &info, 2);  // RC(5' flank):      allele = comp(read[i-16])
//             sites++;
//         }
//     }
//     std::cerr << std::format("[BuildTrie] indexed {} sites (2 patterns each), ~{} nodes used\n",
//                              sites, trie->used / Trie::NODE_SIZE);
// }

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

// recover the canonical core flank pair {n = larger, m = smaller} from the packed
// table coords. inverse of HashSquare: b_prime = b/63, r_prime = (b%63)*65 + r.
std::pair<ll, ll> recover_n_m_prime(ll t, ll b_prime, ll r_prime) {
    ll b = b_prime * 63 + r_prime / 65;   // undo b_prime = b/63
    ll r = r_prime % 65;                  // undo r_prime = (b%63)*65 + r
    return recover_n_m(t, b, r);
}

inline int seqdiff(int largeKmer, int lf, int rf) {
    int c_old = largeKmer & 0x3;
    int rf_old = (largeKmer >> 2) & FLANK_MASK;
    int lf_old = (largeKmer >> LF_SHIFT) & FLANK_MASK;
    
    int diff = 0;
    for (int i = 0; i < LMERADDIFLANK; ++i) {
        int base_lf_new = (lf >> (2 * i)) & 0x3;
        int base_lf_old = (lf_old >> (2 * i)) & 0x3;
        if (base_lf_new != base_lf_old) diff++;
        
        int base_rf_new = (rf >> (2 * i)) & 0x3;
        int base_rf_old = (rf_old >> (2 * i)) & 0x3;
        if (base_rf_new != base_rf_old) diff++;
    }
    return diff;
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
    for(ll b=0; b < MAX_B; b++){
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

void BuildKmerInfo (char **FilterTable, KMerInfo (*kmerinfo)[16][8 << 20]) {
    struct Temp {
        ll t, b, r, b_prime, r_prime;
        short sum;
        Temp (ll _t, ll _b, ll _r, short _sum): t(_t), b(_b), r(_r), sum(_sum) {
            b_prime = b / 63;
            r_prime = (b % 63) * 65 + r;
        }
    } temp (0,0,0,-1);
    // temp.sum = -1;
    short count = 0;
    for (ll b = 0; b < MAX_B; b ++) {
        for (short t = 0; t < 16; t++) {
            if (FilterTable[t][b] < 63) {
                auto sum = int((FilterTable[t][b] >> 6) & 0x03);
                if (sum > temp.sum){
                    ll r = (unsigned char)FilterTable[t][b] & 0x3F;
                    temp = Temp (t, b, r, sum);
                }
            } 
        }
        count ++;
        if (count % 63 == 0){
            count = 0;
            
            if (temp.sum >= 0){
                (*kmerinfo)[temp.t][temp.b_prime].r_prime = temp.r_prime;
            } else {
                (*kmerinfo)[temp.t][temp.b_prime].r_prime = (short)NULL;
            }
            temp.sum = -1;
        }
    }
}

void CallLargeMers(int (*largeKmer)[16][8<<20], KMerInfo (*kmerinfo)[16][8 << 20], string fileName, int fileorder, char* fileBuffer) {
    cerr << std::format("[Task {}] Initializing...\n", fileName);

    SeqParser* seqfile = new SeqParser(fileName, fileBuffer);


    while (seqfile->nextSeq()) {
        LargeBinSeq* sequence = new LargeBinSeq(seqfile->getSeq(1)); // fix to adapt for FASTQ real files.
        ll n, m, lf, rf, c;
        while(sequence->Yieldable()){
            std::tie(n, m, lf, rf, c) = sequence->YieldFlanks();
            HashSquare hs(n, m, EMPTY);

            if (hs.t >= 16) {
                continue;
            }
            if ((*kmerinfo)[hs.t][hs.b_prime].r_prime == hs.r_prime) {
                if (largeKmer[fileorder][hs.t][hs.b_prime] == -1) {
                    largeKmer[fileorder][hs.t][hs.b_prime] = ((lf << LF_SHIFT) | (rf << 2) | c);
                } else if (seqdiff(largeKmer[fileorder][hs.t][hs.b_prime], lf, rf) > 1 || (largeKmer[fileorder][hs.t][hs.b_prime] &0x3) != c) {
                    largeKmer[fileorder][hs.t][hs.b_prime] = -2;
                }
            }
        }
    }
    std::cerr << std::format("[Task {}] finished.\n", fileName);
}
// struct WorkFlowc {
//     WorkFlow {

//     }


// }

std::string full_length_int_2_mer(int val) {
    static const char base[4] = {'A', 'C', 'G', 'T'};
    std::string result;
    result.reserve(16);
    // 从最高位（bit 31-30）到最低位（bit 1-0），共16个碱基
    for (int i = 15; i >= 0; --i) {
        int idx = (val >> (2 * i)) & 0x3;
        result.push_back(base[idx]);
    }
    return result;
}

std::string encode_int32_to_10mer(int val) {
    // static const char base[4] = {'A', 'C', 'G', 'T'};
    // uint32_t uval = static_cast<uint32_t>(val) & 0xFFFFF; // 保留低20位
    // std::string result;
    // result.reserve(10);
    // for (int i = 9; i >= 0; --i) {          // 从第19-18位开始
    //     int idx = (uval >> (2 * i)) & 0x3;
    //     result.push_back(base[idx]);
    // }
    // return result;
    return full_length_int_2_mer(val);
}

int main(int argc, char** argv){
    mem = static_cast<char*>(::operator new(MEM_SIZE, std::align_val_t(ALIGNMENT)));
    // TO DELETE: ::operator delete(mem, std::align_val_t(ALIGNMENT));
    if(!mem){
        throw "Failed to allocate memory!";
        exit(-1);
    }

    
    // FilterTable         buffer
    //|----8GB----|--3GB--|-1GB-|
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

    // FilterTable EsTable buffer
    //|----8GB----|--3GB--|-1GB-|
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

        // FilterTable kMerInfo buffer
        //|----8GB----|--3GB--|-1GB-|
        KMerInfo (*kMerInfoTable)[16][8 << 20] = reinterpret_cast<KMerInfo (*)[16][8 << 20]>(EsStart);
        {
            std::cerr << "Compressing Filtertable...\n";
            
            memset(kMerInfoTable, 0, sizeof(*kMerInfoTable));

            BuildKmerInfo(FilterTable, kMerInfoTable);
            std::cerr <<"finished.\n";      
        }
        
        // int (*largeKmer)[16][16][8<<20] = (int (*)[16][16][8<<20]) memstart;
        int (*largeKmer)[16][8<<20] = new (memstart) int[16][16][8<<20];
        memset(largeKmer, -1, (16ULL * 16 * (8 << 20) * sizeof(int)));

        cerr << largeKmer << endl;

        std::cerr << "Adding working threads: 16 RecallLMWorkers\n";

        // CallLargeMers(largeKmer, kMerInfoTable, std::format("simc{}.fa", 1), 0, bufferStart);
        std::vector <std::thread> RCLthreads;
        for (int i=0; i<16; i++){
            RCLthreads.emplace_back(CallLargeMers, largeKmer, kMerInfoTable, std::format("simc{}.fa", i+1), i, bufferStart + i*BUFFER_SIZE);
        }

        for (auto &i : RCLthreads){
            i.join();
        }
        
        for (int t=0; t < 16; t++) {
            for (int b_prime = 0; b_prime < MAX_BPRIME; b_prime++) {
                if ((*kMerInfoTable)[t][b_prime].r_prime){
                    int multihits = 0;
                    std::vector<int> profile;

                    // summarize all flank info and form a profile
                    for (int f=0; f<16; f++){
                        switch(largeKmer[f][t][b_prime]){
                            case -1: continue; break;
                            case -2: multihits ++; break;
                            default: profile.push_back((largeKmer[f][t][b_prime])); (*kMerInfoTable)[t][b_prime].freq++;
                        }
                        if (multihits > 1) {(*kMerInfoTable)[t][b_prime].del = true; break;}
                    }

                    if (!(*kMerInfoTable)[t][b_prime].del) {
                        for (auto i : profile) {
                            cerr << "[profile]  " << encode_int32_to_10mer(i)  << " - raw - " << std::bitset<32>(i) << endl;
                        }

                        int left_freq[LMERADDIFLANK][4] = {{0}};  
                        int right_freq[LMERADDIFLANK][4] = {{0}}; 
                        for (auto j : profile) {
                            int lf_part = (j >> LF_SHIFT) & FLANK_MASK;  
                            int rf_part = (j >> 2) & FLANK_MASK; 
                            for (int pos = 0; pos < LMERADDIFLANK; ++pos) {
                                int base_l = (lf_part >> (2 * (LMERADDIFLANK - 1 - pos))) & 0x3;  
                                int base_r = (rf_part >> (2 * (LMERADDIFLANK - 1 - pos))) & 0x3;
                                left_freq[pos][base_l]++;
                                right_freq[pos][base_r]++;
                            }
                        }

                        // -2, 6-mer
                        int consensusL = 0;
                        for (int pos = 0; pos < LMERADDIFLANK; ++pos) {
                            int max_idx = std::max_element(left_freq[pos], left_freq[pos] + 4) - left_freq[pos];
                            consensusL = (consensusL << 2) | (max_idx & 0x3);
                        }

                        int consensusR = 0;
                        for (int pos = 0; pos < LMERADDIFLANK; ++pos) {
                            int max_idx = std::max_element(right_freq[pos], right_freq[pos] + 4) - right_freq[pos];
                            consensusR = (consensusR << 2) | (max_idx & 0x3);
                        }

                        cerr << "[consensusL] " << encode_int32_to_10mer(consensusL) << endl;
                        cerr << "[consensusR] " << encode_int32_to_10mer(consensusR) << endl;

                        // depo: store the COMPLETE 16-mer flank into consensusL/R.
                        //   consensusL/R now hold [outer 6-mer | core 10-mer] (32 bits);
                        //   the bare 6-mer consensus is still recoverable as (val >> 20).
                        // core 10-mers (n = larger, paired with L; m = smaller, paired with R)
                        // are recovered from the table coords (t, b_prime, r_prime).
                        auto [n_rec, m_rec] = recover_n_m_prime(t, b_prime, (*kMerInfoTable)[t][b_prime].r_prime);
                        // fullL = [consensusL(6) | n(10)]   forward strand, outer -> core
                        // fullR = [m(10) | consensusR(6)]   revcomp strand, core -> outer
                        (*kMerInfoTable)[t][b_prime].consensusL = ((unsigned long long)consensusL << 20) | (unsigned long long)n_rec;
                        (*kMerInfoTable)[t][b_prime].consensusR = ((unsigned long long)m_rec << 12) | (unsigned long long)consensusR;

                        cerr << "[fullL] " << full_length_int_2_mer((int)(*kMerInfoTable)[t][b_prime].consensusL) << endl;
                        cerr << "[fullR] " << full_length_int_2_mer((int)(*kMerInfoTable)[t][b_prime].consensusR) << endl;
                    }
                }
            }
        }
    }

    std::cerr << "Finished";

    return 0;
}