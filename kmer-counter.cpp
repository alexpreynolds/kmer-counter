#include "kmer-counter.hpp"

const std::string kmer_counter::KmerCounter::client_name = "kmer-counter";
const std::string kmer_counter::KmerCounter::client_version = "1.0";
const std::string kmer_counter::KmerCounter::client_authors = "Alex Reynolds";

int
main(int argc, char** argv)
{
    kmer_counter::KmerCounter kc;

    kc.initialize_command_line_options(argc, argv);
    kc.initialize_kmer_map();
    kc.print_kmer_map(kc.results_kmer_map_stream());
    kc.parse_input_with_emilib_hm();
    
    return EXIT_SUCCESS;
}

void
kmer_counter::KmerCounter::initialize_kmer_map(void)
{
    int i = 0;
    std::uint64_t x, y;
    int k = this->k();
    std::string mer;

    // modified from: https://www.biostars.org/p/18096/#18107
    mer.reserve(k);
    for (x = 0; x < 1ULL<<(2*k); ++x) {
        for (i = 0, y = x; i < k; ++i, y >>= 2)
            mer.push_back("ACGT"[y & 3]);
        this->add_mer_key(mer, this->offset(true));
        mer.clear();
    }
}

void
kmer_counter::KmerCounter::print_kmer_map(FILE* os)
{
    auto map = this->mer_keys();
    for (auto iter = map.begin(); iter != map.end(); ++iter) {
        std::fprintf(os, "%s\t%d\n", iter->first.c_str(), iter->second);
    }    
}

void
kmer_counter::KmerCounter::parse_input_with_emilib_hm(void)
{
    char* buf = NULL;
    size_t buf_len = 0;
    ssize_t buf_read = 0;
    char chr_str[LINE_MAX];
    char start_str[LINE_MAX];
    char stop_str[LINE_MAX];
    char id_str[LINE_MAX];
    char kv_pair[LINE_MAX];
    std::string kv_pairs("");
    emilib::HashMap<std::string, int> mer_counts;
    std::string n("N");

    while ((buf_read = getline(&buf, &buf_len, this->in_stream())) != EOF) {
        std::sscanf(buf, "%s\t%s\t%s\t%s\n", chr_str, start_str, stop_str, id_str);
        std::string seq(id_str);
        std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
        std::deque<char> window(seq.begin(), seq.begin() + this->k());
        // walk over all windows across sequence
        for (size_t i = this->k(); i <= seq.length(); ++i) {
            std::string mer_f(window.begin(), window.end());
            window.pop_front();
            window.push_back(seq[i]);
            std::size_t n_found = mer_f.find(n);
            if (n_found != std::string::npos) {
                continue;
            }
            std::string mer_r(mer_f);
            reverse_complement_string(mer_r);
            if ((mer_counts.count(mer_f) == 0) && (mer_counts.count(mer_r) == 0)) {
                mer_counts[mer_f] = 1;
                // we don't want to add a palindrome twice
                if (mer_f.compare(mer_r) == 0) {
                    continue;
                }
                mer_counts[mer_r] = 1;
            }
            else if ((mer_counts.count(mer_f) == 1) || (mer_counts.count(mer_r) == 1)) {
                mer_counts[mer_f]++;
                mer_counts[mer_r]++;
            }
        }
        // write out all the hits
        kv_pairs.clear();
        for (auto iter = mer_counts.begin(); iter != mer_counts.end(); ++iter) {
            if (iter->second != 0) {
                std::sprintf(kv_pair, "%d:%d ", this->mer_key(iter->first), iter->second);
                mer_counts[iter->first] = 0;
                kv_pairs.append(kv_pair);
            }
        }
        if (kv_pairs.length() > 0) {
            kv_pairs.pop_back();
        }
        std::fprintf(this->results_kmer_count_stream(), "%s\t%s\t%s\t%s\n", chr_str, start_str, stop_str, kv_pairs.c_str());
    }

    // cleanup
    free(buf);
}

std::string
kmer_counter::KmerCounter::client_kmer_counter_opt_string(void)
{
    static std::string _s("k:o:r:hv?");
    return _s;
}

struct option*
kmer_counter::KmerCounter::client_kmer_counter_long_options(void)
{
    static struct option _k = { "k",              required_argument,   NULL,    'k' };
    static struct option _o = { "offset",         required_argument,   NULL,    'o' };
    static struct option _r = { "results-dir",    required_argument,   NULL,    'r' };    
    static struct option _h = { "help",           no_argument,         NULL,    'h' };
    static struct option _v = { "version",        no_argument,         NULL,    'v' };
    static struct option _0 = { NULL,             no_argument,         NULL,     0  };
    static std::vector<struct option> _s;
    _s.push_back(_k);
    _s.push_back(_o);
    _s.push_back(_r);
    _s.push_back(_h);
    _s.push_back(_v);
    _s.push_back(_0);
    return &_s[0];
}

void
kmer_counter::KmerCounter::initialize_command_line_options(int argc, char** argv)
{
    int client_long_index;
    int client_opt = getopt_long(argc,
                                 argv,
                                 this->client_kmer_counter_opt_string().c_str(),
                                 this->client_kmer_counter_long_options(),
                                 &client_long_index);
    int _k = -1;
    int _offset = -1;

    opterr = 0; /* disable error reporting by GNU getopt */
    
    while (client_opt != -1) {
        switch (client_opt) {
        case 'k':
            std::sscanf(optarg, "%d", &_k);
            this->k(_k);
            break;
        case 'o':
            std::sscanf(optarg, "%d", &_offset);
            this->offset(_offset);
            break;
        case 'r':
            this->results_dir(optarg);
            break;
        case 'h':
            this->print_usage(stdout);
            std::exit(EXIT_SUCCESS);
        case 'v':
            this->print_version(stdout);
            std::exit(EXIT_SUCCESS);
        case '?':
            this->print_usage(stdout);
            std::exit(EXIT_SUCCESS);
        default:
            break;
        }
        client_opt = getopt_long(argc,
                                 argv,
                                 this->client_kmer_counter_opt_string().c_str(),
                                 this->client_kmer_counter_long_options(),
                                 &client_long_index);
    }
    
    if (optind < argc) {
        do {
            if (this->input_fn().empty()) {
                this->input_fn(argv[optind]);
                this->initialize_in_stream();
            }
            else {
                std::fprintf(stderr, "Warning: Ignoring additional input file [%s]\n", argv[optind]);
            }
        }
        while (++optind < argc);
    }

    if (this->k() == -1) {
        std::fprintf(stderr, "Error: Specify k value\n");
        this->print_usage(stderr);
        std::exit(ENODATA);
    }

    if (this->offset() == -1) {
        std::fprintf(stderr, "Error: Specify offset value\n");
        this->print_usage(stderr);
        std::exit(ENODATA);
    }

    if (this->results_dir().empty()) {
        std::fprintf(stderr, "Error: Specify results dir value\n");
        this->print_usage(stderr);
        std::exit(ENODATA);        
    }
    else {
        this->results_dir_mode(0755);
        if (this->initialize_result_dir(this->results_dir(), this->results_dir_mode())) {
            this->initialize_kmer_count_stream();
            this->initialize_kmer_map_stream();
        }
        else {
            std::fprintf(stderr, "Error: Could not create specified path [%s] with specified mode [%o]\n", this->results_dir().c_str(), this->results_dir_mode());
            this->print_usage(stderr);
            std::exit(EINVAL);
        }
    }
}

std::string
kmer_counter::KmerCounter::client_kmer_counter_name(void)
{
    static std::string _s(kmer_counter::KmerCounter::client_name);
    return _s;
}

std::string
kmer_counter::KmerCounter::client_kmer_counter_version(void)
{
    static std::string _s(kmer_counter::KmerCounter::client_version);
    return _s;
}

std::string
kmer_counter::KmerCounter::client_kmer_counter_authors(void)
{
    static std::string _s(kmer_counter::KmerCounter::client_authors);
    return _s;
}

std::string
kmer_counter::KmerCounter::client_kmer_counter_usage(void)
{
    static std::string _s("\n"                                          \
                          "  Usage:\n"                                  \
                          "\n"                                          \
                          "  $ kmer_counter [options] input\n");
    return _s;
}

std::string
kmer_counter::KmerCounter::client_kmer_counter_description(void)
{
    static std::string _s("  Count kmers in BED file and report counts for\n" \
                          "  each element, from specified offset. Write a\n" \
                          "  key-count file and key table as output.\n");
    return _s;
}

std::string
kmer_counter::KmerCounter::client_kmer_counter_io_options(void)
{
    static std::string _s("  General Options:\n\n"              \
                          "  --k=n            K-value for kmer length (integer)\n" \
                          "  --offset=n       Offset (integer)\n" \
                          "  --results-dir=s  Results directory (string)\n");
    return _s;
}

std::string
kmer_counter::KmerCounter::client_kmer_counter_general_options(void)
{
    static std::string _s("  Process Flags:\n\n"                        \
                          "  --help                  Show this usage message\n" \
                          "  --version               Show binary version\n");
    return _s;
}

void
kmer_counter::KmerCounter::print_usage(FILE* wo_stream)
{
    std::fprintf(wo_stream,
                 "%s\n"                                              \
                 "  version: %s\n"                                   \
                 "  author:  %s\n"                                   \
                 "%s\n"                                              \
                 "%s\n"                                              \
                 "%s\n"                                              \
                 "%s\n",
                 this->client_kmer_counter_name().c_str(),
                 this->client_kmer_counter_version().c_str(),
                 this->client_kmer_counter_authors().c_str(),
                 this->client_kmer_counter_usage().c_str(),
                 this->client_kmer_counter_description().c_str(),
                 this->client_kmer_counter_io_options().c_str(),
                 this->client_kmer_counter_general_options().c_str());
}

void
kmer_counter::KmerCounter::print_version(FILE* wo_stream)
{
    std::fprintf(wo_stream,
                 "%s\n"                                              \
                 "  version: %s\n"                                   \
                 "  author:  %s\n",
                 this->client_kmer_counter_name().c_str(),
                 this->client_kmer_counter_version().c_str(),
                 this->client_kmer_counter_authors().c_str());
}
