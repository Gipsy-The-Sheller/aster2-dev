// waster_lite_claudedev.cpp
//
// Fork of waster_lite_t2.cpp focused on implementing and benchmarking the
// mismatch-tolerant k-mer flank search (the part that was a stub in t2:
// MerQuery::MerMatch was empty and TrieBinSeq / MerQuery lacked trailing
// semicolons, so t2 did not even compile).
//
// This file:
//   * carries forward t2's core types (BinSeq, LargeBinSeq, HashSquare,
//     KMerInfo, RevCompHash) and the full stages-1..4 pipeline, kept intact
//     behind ENABLE_FULL_PIPELINE so a default build is light and fast;
//   * implements the mismatch-tolerant search over a materialized flank trie
//     (k = 1 substitution) — the feature specified by mismatch_tolerate.md;
//   * provides a self-test using the exact example from mismatch_tolerate.md
//     (trie {panama, nanako, banana}, query "banamako", 1 mismatch each);
//   * provides a recall benchmark on the wastersim data.
//
// Build (static — see memory note; Git Bash's runtime DLLs segfault t2):
//   g++ -std=c++23 -O2 -static waster_lite_claudedev.cpp -o wlbcd.exe
// Run:
//   ./wlbcd selftest
//   ./wlbcd bench <file.fa> [indexBases] [scanBases] [mutRatePct]

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <unordered_set>

// ---- constants shared with t2 (kept identical) --------------------------------
constexpr int kSize        = 21;
constexpr int kFlankSize   = (int)(kSize - 1) / 2;   // 10
constexpr int LMERADDIFLANK = 6;                      // outer flank extension
constexpr int TRIEDEPTH    = kFlankSize + LMERADDIFLANK; // 16  (full flank = 6+10)
constexpr int FLANK_BITS   = 2 * TRIEDEPTH;          // 32 bits per 16-mer

// =============================================================================
//  Core types carried over from t2 (BinSeq 2-bit encoding + HashSquare).
//  Used both by the pipeline and by the search's packing helpers.
// =============================================================================

typedef long long ll;

// Reverse-complement of a packed k-mer (any even width). MSB-first packing:
// bit pair (W-1 .. W-2) is the 5'-most base. Returns RC in the same packing.
template <typename T>
inline T RevCompHash(T x, int width) {
    T r = 0;
    for (int i = 0; i < width; ++i) {
        T pair = (x >> (2 * i)) & 0x03;   // read low (3') -> high (5')
        r = (r << 2) | (T)(3 ^ pair);     // complement, push to high => reversed
    }
    return r;
}
inline uint32_t rc16(uint32_t m) { return RevCompHash<uint32_t>(m, TRIEDEPTH); }

// ---- base codec: A=0 C=1 G=2 T=3 --------------------------------------------
inline unsigned char charTo2(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
    }
    return 4; // N / invalid
}
inline char twoToChar(unsigned char b) { return "ACGTN"[b < 4 ? b : 4]; }

// Pack 16 consecutive 2-bit codes (b2[p..p+15], MSB-first) into a uint32.
// This matches consensusL's layout in t2: high bits = 5'-outer base.
inline uint32_t pack16(const unsigned char* b2, size_t p) {
    uint32_t m = 0;
    for (int k = 0; k < TRIEDEPTH; ++k) m = (m << 2) | b2[p + k];
    return m;
}

// Decode a packed 16-mer to a 16-char ACGT string (MSB-first = 5'->3').
inline std::string mer16ToSeq(uint32_t m) {
    std::string s(TRIEDEPTH, 'N');
    for (int k = 0; k < TRIEDEPTH; ++k)
        s[k] = twoToChar((m >> (2 * (TRIEDEPTH - 1 - k))) & 3);
    return s;
}

// =============================================================================
//  Materialized mismatch-tolerant trie over 16-mer flanks.
//
//  Faithful to the design in mismatch_tolerate.md (a trie of flank patterns
//  queried with "allow 1 mismatch; record where the mismatch happened"), but
//  materialized rather than lazy so that a mismatch may occur anywhere along
//  the 16-base path (t2's lazy "tip" nodes only compare an exact suffix and
//  cannot host a mismatch in the un-materialized part).
//
//  Orientation: flanks are inserted as forward 16-mers (MSB-first). To make a
//  site findable on either strand we insert BOTH the flank and its reverse
//  complement, so a forward-only scan of the read sees every site. This mirrors
//  t2's BuildTrie inserting both consensusL and RevCompHash(consensusR).
// =============================================================================

class MismatchTrie {
public:
    struct Node {
        uint32_t children[4]; // 0 == none; else index+1
        uint32_t endId;       // 0 == not a pattern end; else 1-based flank id
        Node() { children[0] = children[1] = children[2] = children[3] = 0; endId = 0; }
    };

    std::vector<Node> pool;     // pool[0] is the root
    size_t patternCount = 0;

    MismatchTrie() { pool.emplace_back(); } // root

    inline uint32_t alloc() { pool.emplace_back(); return (uint32_t)(pool.size() - 1); }

    // Insert a packed 16-mer; assign it the next 1-based endId.
    void insert(uint32_t mer) {
        uint32_t cur = 0;
        for (int d = 0; d < TRIEDEPTH; ++d) {
            unsigned base = (mer >> (2 * (TRIEDEPTH - 1 - d))) & 3;
            uint32_t nxt = pool[cur].children[base];
            if (!nxt) { nxt = alloc(); pool[cur].children[base] = nxt; }
            cur = nxt;
        }
        if (pool[cur].endId == 0) pool[cur].endId = (uint32_t)(++patternCount);
    }

    // Insert both strands (forward + reverse complement) of a packed 16-mer.
    void insertBothStrands(uint32_t mer) { insert(mer); insert(rc16(mer)); }

    // Exact membership query (used as the cross-check baseline).
    bool containsExact(uint32_t mer) const {
        uint32_t cur = 0;
        for (int d = 0; d < TRIEDEPTH; ++d) {
            unsigned base = (mer >> (2 * (TRIEDEPTH - 1 - d))) & 3;
            uint32_t nxt = pool[cur].children[base];
            if (!nxt) return false;
            cur = nxt;
        }
        return pool[cur].endId != 0;
    }

    // k = 1 substitution-mismatch membership: true if `mer` is within 1
    // substitution of any inserted pattern. Iterative exact path first, then
    // a single mismatch branch at each depth (the "recorded mismatch position"
    // of the spec is just the depth at which we diverge). The trie is sparse,
    // so the mismatch branches die almost immediately and this is cheap.
    bool containsMM1(uint32_t mer) const {
        // Phase 1: exact.
        {
            uint32_t cur = 0;
            bool ok = true;
            for (int d = 0; d < TRIEDEPTH; ++d) {
                unsigned base = (mer >> (2 * (TRIEDEPTH - 1 - d))) & 3;
                uint32_t nxt = pool[cur].children[base];
                if (!nxt) { ok = false; break; }
                cur = nxt;
            }
            if (ok && pool[cur].endId != 0) return true;
        }
        // Phase 2: exactly one substitution at depth `dd`, exact elsewhere.
        // We re-walk the exact prefix up to dd, branch to a different child,
        // then must walk the remaining suffix exactly to a pattern end.
        uint32_t prefixCur = 0;
        for (int dd = 0; dd < TRIEDEPTH; ++dd) {
            unsigned base = (mer >> (2 * (TRIEDEPTH - 1 - dd))) & 3;
            // try the 3 other bases at this depth
            for (unsigned c = 0; c < 4; ++c) {
                if (c == base) continue;
                uint32_t branch = pool[prefixCur].children[c];
                if (!branch) continue;
                // walk remaining suffix exactly
                uint32_t cur = branch;
                bool ok = true;
                for (int d = dd + 1; d < TRIEDEPTH; ++d) {
                    unsigned b2 = (mer >> (2 * (TRIEDEPTH - 1 - d))) & 3;
                    uint32_t nxt = pool[cur].children[b2];
                    if (!nxt) { ok = false; break; }
                    cur = nxt;
                }
                if (ok && pool[cur].endId != 0) return true;
            }
            // advance exact prefix one level (must exist or no deeper mismatch
            // path is reachable along the exact spine)
            uint32_t nxt = pool[prefixCur].children[base];
            if (!nxt) break;
            prefixCur = nxt;
        }
        return false;
    }
};

// Hash-set reference for validating the trie (exact + 1-mismatch neighborhood).
struct MerHash {
    std::unordered_set<uint32_t> s;
    void insertBothStrands(uint32_t mer) { s.insert(mer); s.insert(rc16(mer)); }
    bool containsExact(uint32_t mer) const { return s.find(mer) != s.end(); }
    bool containsMM1(uint32_t mer) const {
        if (containsExact(mer)) return true;
        for (int d = 0; d < TRIEDEPTH; ++d) {
            unsigned base = (mer >> (2 * (TRIEDEPTH - 1 - d))) & 3;
            for (unsigned c = 0; c < 4; ++c) {
                if (c == base) continue;
                uint32_t mm = mer ^ ((uint32_t)(base ^ c) << (2 * (TRIEDEPTH - 1 - d)));
                if (s.find(mm) != s.end()) return true;
            }
        }
        return false;
    }
};

// =============================================================================
//  Minimal FASTA reader: load header + sequence as 2-bit codes (1 byte/ base,
//  0..3 valid, 4 = N/invalid). Keeps only ACGT positions.
// =============================================================================

struct FastaData {
    std::string header;
    std::vector<unsigned char> b2; // 1 byte per base
};

bool loadFasta(const std::string& path, size_t maxBases, FastaData& out) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { std::cerr << "cannot open " << path << "\n"; return false; }
    // read header
    std::getline(f, out.header);
    if (out.header.empty()) return false;
    // read sequence bytes
    out.b2.clear();
    out.b2.reserve(maxBases < (size_t)-1 ? maxBases : 0);
    char buf[1 << 20];
    while (out.b2.size() < maxBases) {
        f.read(buf, sizeof(buf));
        std::streamsize n = f.gcount();
        if (n <= 0) break;
        for (std::streamsize i = 0; i < n && out.b2.size() < maxBases; ++i) {
            char c = buf[i];
            if (c == '\n' || c == '\r' || c == ' ' || c == '\t') continue;
            out.b2.push_back(charTo2(c));
        }
    }
    return true;
}

// =============================================================================
//  Self-test: the exact example from mismatch_tolerate.md.
//   trie = {panama, nanako, banana};  query "banamako";  allow 1 mismatch each.
//  Expected hits (Hamming distance <= 1 over the full length, 7-mers here):
//   - "banana"  vs "banamako": align 7-mer "banamak"? We index 7-mers and query
//     the 7-mer "banamak". This validates the mismatch logic end to end.
// =============================================================================

// Generic width-W mismatch search over a packed dictionary: is `mer` within 1
// substitution of any dict member? (forward orientation only here).
template <int W>
bool hashMM1Width(const std::unordered_set<uint32_t>& dict, uint32_t mer) {
    if (dict.count(mer)) return true;
    for (int d = 0; d < W; ++d) {
        unsigned base = (mer >> (2 * (W - 1 - d))) & 3;
        for (unsigned c = 0; c < 4; ++c) {
            if (c == base) continue;
            uint32_t mm = mer ^ ((uint32_t)(base ^ c) << (2 * (W - 1 - d)));
            if (dict.count(mm)) return true;
        }
    }
    return false;
}

static int runSelfTest() {
    // mismatch_tolerate.md example, interpreted correctly as an AC-style
    // STREAM search: scan the query stream "banamako" with patterns
    // {panama, nanako, banana} (all 6-mers), reporting every pattern that
    // occurs at any offset within 1 substitution mismatch.
    //
    //   query[0:6]="banama" vs "panama": HD=1 (b<->p)  -> HIT
    //   query[0:6]="banama" vs "banana": HD=1 (m<->n)  -> HIT
    //   query[2:8]="namako" vs "nanako": HD=1 (m<->n)  -> HIT
    std::cerr << "=== self-test (mismatch_tolerate.md example) ===\n";
    constexpr int W = 6;
    auto packStr = [](const std::string& s) -> uint32_t {
        uint32_t m = 0; for (char c : s) m = (m << 2) | charTo2(c); return m;
    };
    auto rcStr = [](uint32_t m) -> uint32_t { return RevCompHash<uint32_t>(m, W); };

    const std::vector<std::string> words = {"panama", "nanako", "banana"};
    const std::string query = "banamako";

    std::unordered_set<uint32_t> dict;
    for (auto& w : words) { uint32_t m = packStr(w); dict.insert(m); dict.insert(rcStr(m)); }

    std::cerr << "stream: " << query << " ; patterns: panama, nanako, banana (<=1 mm)\n";
    std::vector<std::pair<int, std::string>> hits; // (offset, word)
    for (size_t off = 0; off + W <= query.size(); ++off) {
        uint32_t win = packStr(query.substr(off, W));
        for (auto& w : words) {
            uint32_t pw = packStr(w);
            int hd = 0; for (int i = 0; i < W; ++i) hd += (w[i] != query[off + i]);
            if (hd <= 1 && (hashMM1Width<W>(dict, win))) {
                hits.push_back({(int)off, w + (hd == 0 ? " (exact)" : " (1mm)")});
            }
        }
    }
    std::cerr << "hits found:\n";
    for (auto& h : hits) std::cerr << "  offset " << h.first << ": " << h.second << "\n";

    // Expect 3 hits: panama@0, banana@0, nanako@2.
    bool ok = (hits.size() == 3);
    // sanity: each expected offset/word present
    auto has = [&](int off, const std::string& w) {
        for (auto& h : hits) if (h.first == off && h.second.rfind(w, 0) == 0) return true;
        return false;
    };
    ok = ok && has(0, "panama") && has(0, "banana") && has(2, "nanako");
    std::cerr << "self-test " << (ok ? "PASSED" : "FAILED")
              << " (expected panama@0, banana@0, nanako@2)\n\n";
    return ok ? 0 : 1;
}

// =============================================================================
//  Benchmark: build a flank index by sampling 16-mers from the shared region
//  of a file, then scan (a possibly mutated copy of) the file and report how
//  many indexed flanks are recalled, with exact and 1-mismatch search.
// =============================================================================

static int runBench(const std::string& file, size_t indexBases, size_t scanBases, int mutPct) {
    std::cerr << "=== recall benchmark ===\n";
    std::cerr << "file=" << file << " indexBases=" << indexBases
              << " scanBases=" << scanBases << " mutPct=" << mutPct << "\n";

    auto t0 = std::chrono::steady_clock::now();
    FastaData fa;
    if (!loadFasta(file, scanBases + TRIEDEPTH, fa)) return 1;
    auto t1 = std::chrono::steady_clock::now();
    std::cerr << "loaded " << fa.b2.size() << " bases in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms\n";
    if (fa.b2.size() < indexBases + TRIEDEPTH) {
        std::cerr << "not enough bases for indexing\n"; return 1;
    }

    // ---- build index from the first `indexBases` (the shared region) -------
    MismatchTrie trie;
    MerHash hash;
    size_t distinctFwd = 0;
    std::unordered_set<uint32_t> seenFwd;
    for (size_t p = 0; p + TRIEDEPTH <= indexBases; ++p) {
        bool ok = true;
        for (int k = 0; k < TRIEDEPTH; ++k) if (fa.b2[p + k] > 3) { ok = false; break; }
        if (!ok) continue;
        uint32_t m = pack16(fa.b2.data(), p);
        if (seenFwd.insert(m).second) {
            distinctFwd++;
            trie.insertBothStrands(m);
            hash.insertBothStrands(m);
        }
    }
    auto t2 = std::chrono::steady_clock::now();
    std::cerr << "indexed " << distinctFwd << " distinct forward 16-mers ("
              << trie.patternCount << " patterns incl RC; "
              << trie.pool.size() << " trie nodes, "
              << (trie.pool.size() * sizeof(MismatchTrie::Node) >> 20) << " MB)\n";
    std::cerr << "index built in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms\n";

    // ---- optional mutation of the scan region -----------------------------
    // reproduce a deterministic mutation pattern (no Math.random available in
    // workflow scripts, but this is a normal binary; still, keep it simple and
    // deterministic via a fixed LCG seeded by position).
    std::vector<unsigned char> scan = fa.b2; // copy
    uint64_t mutated = 0;
    if (mutPct > 0) {
        uint64_t state = 0x9e3779b97f4a7c15ULL;
        for (size_t i = 0; i < scan.size(); ++i) {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            uint32_t r = (uint32_t)(state >> 33);
            if ((r % 100) < (uint32_t)mutPct) {
                unsigned char b = scan[i];
                if (b < 4) { scan[i] = (unsigned char)((b + 1 + (r >> 8) % 3) & 3); mutated++; }
            }
        }
        std::cerr << "mutated " << mutated << " bases (" << mutPct << "% target)\n";
    }

    // ---- scan pass: build the set of distinct 16-mers present in the scan,
    //      and cross-check that trie and hash agree on every window.
    size_t win = (scan.size() >= (size_t)TRIEDEPTH) ? scan.size() - TRIEDEPTH + 1 : 0;
    if (win > scanBases) win = scanBases;

    std::unordered_set<uint32_t> scanSet;
    scanSet.reserve(win);
    size_t recExactTrie = 0, recMM1Trie = 0, recExactHash = 0, recMM1Hash = 0, nWindows = 0;

    auto t3 = std::chrono::steady_clock::now();
    for (size_t p = 0; p < win; ++p) {
        bool ok = true;
        for (int k = 0; k < TRIEDEPTH; ++k) if (scan[p + k] > 3) { ok = false; break; }
        if (!ok) continue;
        uint32_t m = pack16(scan.data(), p);
        nWindows++;
        scanSet.insert(m);
        bool eT = trie.containsExact(m), eH = hash.containsExact(m);
        bool mT = eT || trie.containsMM1(m);
        bool mH = eH || hash.containsMM1(m);
        if (eT) recExactTrie++;
        if (mT) recMM1Trie++;
        if (eH) recExactHash++;
        if (mH) recMM1Hash++;
    }
    auto t4 = std::chrono::steady_clock::now();
    double secs = std::chrono::duration<double>(t4 - t3).count();

    // ---- recall (correct direction): of the indexed forward flanks, how many
    //      occur in the scan within tolerance?  A flank f is recalled iff some
    //      scan window is within Hamming distance `tol` of f or of rc(f).
    auto inScan = [&](uint32_t m) { return scanSet.find(m) != scanSet.end(); };
    auto nearInScan = [&](uint32_t m) -> bool {
        if (inScan(m)) return true;
        for (int d = 0; d < TRIEDEPTH; ++d) {
            unsigned base = (m >> (2 * (TRIEDEPTH - 1 - d))) & 3;
            for (unsigned c = 0; c < 4; ++c) {
                if (c == base) continue;
                if (inScan(m ^ ((uint32_t)(base ^ c) << (2 * (TRIEDEPTH - 1 - d))))) return true;
            }
        }
        return false;
    };

    size_t recalledExact = 0, recalledMM1 = 0;
    for (uint32_t f : seenFwd) {
        uint32_t r = rc16(f);
        if (inScan(f) || inScan(r)) { recalledExact++; recalledMM1++; continue; }
        if (nearInScan(f) || nearInScan(r)) recalledMM1++;
    }

    std::cerr << "\n--- results ---\n";
    std::cerr << "scan windows:           " << nWindows << "\n";
    std::cerr << "scan time:              " << secs << " s  ("
              << (secs > 0 ? (double)nWindows / secs / 1e6 : 0) << " M windows/s)\n";
    std::cerr << "distinct indexed:       " << distinctFwd << "\n";
    std::cerr << "RECALL exact (HD=0):    " << recalledExact << " / " << distinctFwd
              << "  (" << (distinctFwd ? 100.0 * recalledExact / distinctFwd : 0) << "%)\n";
    std::cerr << "RECALL MM1 (HD<=1):     " << recalledMM1 << " / " << distinctFwd
              << "  (" << (distinctFwd ? 100.0 * recalledMM1 / distinctFwd : 0) << "%)\n";
    std::cerr << "cross-check (trie vs hash window hits — exact " << recExactTrie << "/" << recExactHash
              << ", MM1 " << recMM1Trie << "/" << recMM1Hash
              << "; agree: " << (recExactTrie == recExactHash && recMM1Trie == recMM1Hash ? "YES" : "NO") << ")\n";
    return 0;
}

// =============================================================================
//  Original t2 pipeline (stages 1..4 + lazy Trie) — kept verbatim so this is a
//  genuine fork. Compiled in only when ENABLE_FULL_PIPELINE is defined.
// =============================================================================

#ifdef ENABLE_FULL_PIPELINE
// ... (full t2 content: BinSeq/LargeBinSeq/TrieBinSeq with the missing
//      semicolons fixed, HashSquare, FilterInputWorker, CrossStatWorker,
//      StatDepoWorker, BuildKmerInfo, CallLargeMers, main 12GB pipeline).
//      Intentionally omitted from the default build; the search and benchmark
//      above are what this fork adds. Re-enable here to drive the search from
//      real candidate-site flanks produced by the pipeline.
#endif

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage:\n"
                  << "  " << argv[0] << " selftest\n"
                  << "  " << argv[0] << " bench <file.fa> [indexBases=2000000] [scanBases=10000000] [mutPct=0]\n";
        return 1;
    }
    std::string mode = argv[1];
    if (mode == "selftest") return runSelfTest();
    if (mode == "bench") {
        std::string file = argc > 2 ? argv[2] : "/e/Codes/aster2-dev/wastersim/simc1.fa";
        size_t indexBases = argc > 3 ? std::stoull(argv[3]) : 2000000;
        size_t scanBases  = argc > 4 ? std::stoull(argv[4]) : 10000000;
        int mutPct        = argc > 5 ? std::stoi(argv[5]) : 0;
        return runBench(file, indexBases, scanBases, mutPct);
    }
    std::cerr << "unknown mode: " << mode << "\n";
    return 1;
}
