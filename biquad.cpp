#include <bits/stdc++.h>
using namespace std;

static inline string trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static inline double uni01(mt19937_64& g) {
    return (g() >> 11) * (1.0 / 9007199254740992.0);
}

static inline double nCk_fast(long long n, int k) {
    if (k < 0) return 0.0;
    if (k == 0) return 1.0;
    if (n < k) return 0.0;

    switch (k) {
        case 1: return static_cast<double>(n);
        case 2: return 0.5 * static_cast<double>(n) * static_cast<double>(n - 1);
        case 3: return static_cast<double>(n) * static_cast<double>(n - 1) *
                       static_cast<double>(n - 2) / 6.0;
        case 4: return static_cast<double>(n) * static_cast<double>(n - 1) *
                       static_cast<double>(n - 2) * static_cast<double>(n - 3) / 24.0;
        default: {
            k = min(k, static_cast<int>(n - k));
            double r = 1.0;
            for (int i = 1; i <= k; i++) {
                r *= static_cast<double>(n - (k - i));
                r /= static_cast<double>(i);
            }
            return r;
        }
    }
}

struct Reaction {
    vector<pair<int,int>> react;
    vector<pair<int,int>> delta;
    double k = 0.0;
};

struct Model {
    vector<string> species;
    unordered_map<string,int> idx;
    vector<Reaction> rxns;
};

static unordered_map<string,int> parse_pairs_side(const string& side) {
    unordered_map<string,int> cnt;
    string s = trim(side);
    if (s.empty()) return cnt;

    vector<string> tok;
    {
        istringstream iss(s);
        string t;
        while (iss >> t) tok.push_back(t);
    }

    if (tok.size() % 2 != 0) {
        throw runtime_error("Invalid side (species/coefficient pairs expected): " + side);
    }

    for (size_t i = 0; i + 1 < tok.size(); i += 2) {
        cnt[tok[i]] += stoi(tok[i + 1]);
    }
    return cnt;
}

static Model parse_lambda_r(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open " + path);

    vector<tuple<string,string,double>> raw;
    unordered_set<string> sp;
    string line;

    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vector<string> parts;
        string cur;
        for (char c : line) {
            if (c == ':') {
                parts.push_back(trim(cur));
                cur.clear();
            } else {
                cur.push_back(c);
            }
        }
        parts.push_back(trim(cur));

        if (parts.size() != 3) {
            throw runtime_error("Invalid reaction line: " + line);
        }

        string lhs = parts[0];
        string rhs = parts[1];
        double k = stod(parts[2]);

        auto L = parse_pairs_side(lhs);
        auto R = parse_pairs_side(rhs);

        for (auto& kv : L) sp.insert(kv.first);
        for (auto& kv : R) sp.insert(kv.first);

        raw.emplace_back(lhs, rhs, k);
    }

    Model M;
    M.species.assign(sp.begin(), sp.end());
    sort(M.species.begin(), M.species.end());

    for (int i = 0; i < static_cast<int>(M.species.size()); i++) {
        M.idx[M.species[i]] = i;
    }

    for (auto& t : raw) {
        string lhs, rhs;
        double k;
        tie(lhs, rhs, k) = t;

        auto L = parse_pairs_side(lhs);
        auto R = parse_pairs_side(rhs);

        Reaction rx;
        rx.k = k;

        for (auto& kv : L) {
            rx.react.push_back({M.idx[kv.first], kv.second});
        }

        unordered_map<int,int> d;
        for (auto& kv : R) d[M.idx[kv.first]] += kv.second;
        for (auto& kv : L) d[M.idx[kv.first]] -= kv.second;

        for (auto& kv : d) {
            if (kv.second != 0) rx.delta.push_back({kv.first, kv.second});
        }

        M.rxns.push_back(rx);
    }

    return M;
}

static unordered_map<string,long long> parse_lambda_in(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open " + path);

    unordered_map<string,long long> init;
    string line;

    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        istringstream iss(line);
        string name;
        long long val;
        iss >> name >> val;

        if (!iss) {
            throw runtime_error("Invalid input line: " + line);
        }

        init[name] = val;
    }

    return init;
}

static inline double propensity(const Reaction& r, const vector<long long>& x) {
    double a = r.k;
    for (auto [si, sto] : r.react) {
        double c = nCk_fast(x[si], sto);
        if (c == 0.0) return 0.0;
        a *= c;
    }
    return a;
}

static void dump_state_in_format(const string& path, const Model& M, const vector<long long>& x) {
    ofstream out(path);
    if (!out) throw runtime_error("Cannot write state file: " + path);

    for (int i = 0; i < static_cast<int>(M.species.size()); i++) {
        if (x[i] != 0) {
            out << M.species[i] << ' ' << x[i] << " N\n";
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cout << "Usage: " << argv[0]
             << " problem.r problem.in [scale] [state_out.in] [tmax] [max_events]\n";
        return 0;
    }

    const double SCALE = (argc >= 4) ? atof(argv[3]) : 1.0;
    const string state_out = (argc >= 5) ? argv[4] : "final_state.in";
    const double TMAX = (argc >= 6) ? atof(argv[5]) : 1e9;
    const long long MAX_EVENTS = (argc >= 7) ? atoll(argv[6]) : 1000000000LL;

    Model M = parse_lambda_r(argv[1]);
    auto init = parse_lambda_in(argv[2]);

    vector<long long> x(M.species.size(), 0);
    for (auto& kv : init) {
        if (M.idx.count(kv.first)) {
            x[M.idx[kv.first]] = kv.second;
        }
    }

    mt19937_64 gen(1);
    vector<double> a(M.rxns.size());

    double t = 0.0;
    long long events = 0;

    while (events < MAX_EVENTS) {
        double a0 = 0.0;

        for (int i = 0; i < static_cast<int>(M.rxns.size()); i++) {
            a[i] = propensity(M.rxns[i], x);
            a0 += a[i];
        }

        if (a0 <= 0.0) break;

        double r = uni01(gen) * a0;
        double cum = 0.0;
        int chosen = 0;

        for (int i = 0; i < static_cast<int>(a.size()); i++) {
            cum += a[i];
            if (cum >= r) {
                chosen = i;
                break;
            }
        }

        double u = uni01(gen);
        if (u <= 0.0) u = 1e-16;
        t += -log(u) / a0;

        for (auto [si, dv] : M.rxns[chosen].delta) {
            x[si] += dv;
            if (x[si] < 0) {
                cerr << "Error: negative count for species " << M.species[si] << "\n";
                return 2;
            }
        }

        ++events;
    }

    cout << fixed << setprecision(6);
    cout << "Stopped at t = " << t << ", events = " << events;
    if (events >= MAX_EVENTS) cout << " (hit max_events)";
    else if (t >= TMAX) cout << " (hit tmax)";
    else cout << " (no more reactions)";
    cout << "\n";

    if (events >= MAX_EVENTS) {
        cerr << "Warning: simulation hit max_events before settling.\n";
    }

    long long rawY  = M.idx.count("Y")  ? x[M.idx["Y"]]  : 0;
    long long rawR1 = M.idx.count("R1") ? x[M.idx["R1"]] : 0;
    long long rawR2 = M.idx.count("R2") ? x[M.idx["R2"]] : 0;

    cout << "Y  = " << rawY;
    if (SCALE != 1.0) cout << "  (scaled = " << (rawY / SCALE) << ")";
      cout << "\n";

    cout << "R1 = " << rawR1;
    if (SCALE != 1.0) cout << "  (scaled = " << (rawR1 / SCALE) << ")";
      cout << "\n";

    cout << "R2 = " << rawR2;
    if (SCALE != 1.0) cout << "  (scaled = " << (rawR2 / SCALE) << ")";
      cout << "\n";

    dump_state_in_format(state_out, M, x);
    cout << "State written to: " << state_out << "\n";
    return 0;
}