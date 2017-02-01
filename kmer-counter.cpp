#include "kmer-counter.hpp"

const std::string kmer_counter::KmerCounter::client_name = "kmer-counter";
const std::string kmer_counter::KmerCounter::client_version = "1.0";
const std::string kmer_counter::KmerCounter::client_authors = "Alex Reynolds";

int
main(int argc, char** argv)
{
    kmer_counter::KmerCounter kc;

    kc.initialize_command_line_options(argc, argv);
    kc.parse_input();
    
    return EXIT_SUCCESS;
}

void
kmer_counter::KmerCounter::parse_input(void)
{
    char* buf = NULL;
    size_t buf_len = 0;
    ssize_t buf_read = 0;
    char chr_str[LINE_MAX];
    char start_str[LINE_MAX];
    char stop_str[LINE_MAX];
    char id_str[LINE_MAX];
    uint64_t start_val = 0;
    uint64_t stop_val = 0;
    emilib::HashMap<std::string, int> mer_counts;
    std::string n("N");
    
    while ((buf_read = getline(&buf, &buf_len, this->get_in_stream())) != EOF) {
        std::sscanf(buf, "%s\t%s\t%s\t%s\n", chr_str, start_str, stop_str, id_str);
        std::sscanf(start_str, "%" SCNu64, &start_val);
        std::sscanf(stop_str, "%" SCNu64, &stop_val);
        std::string seq(id_str);
        std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
        std::deque<char> window(seq.begin(), seq.begin() + this->get_k());
        for (size_t i = this->get_k(); i <= seq.length(); ++i) {
            std::string mer_f(window.begin(), window.end());
            window.pop_front();
            window.push_back(seq[i]);
            std::size_t n_found = mer_f.find(n);
            if (n_found != std::string::npos) {
                continue;
            }
            std::string mer_r(mer_f);
            reverse_complement_string(mer_r);
            if ((!mer_counts[mer_f]) || (!mer_counts[mer_r])) {
                mer_counts[mer_f] = 1;
                mer_counts[mer_r] = 1;
            }
            else if ((mer_counts[mer_f]) || (mer_counts[mer_r])) {
                mer_counts[mer_f]++;
                mer_counts[mer_r]++;
            }
        }
        for (auto iter = mer_counts.begin(); iter != mer_counts.end(); ++iter) {
            std::fprintf(stdout, "[%s:%d]\n", iter->first.c_str(), iter->second);
        }
    }
    
    free(buf);
}

std::string
kmer_counter::KmerCounter::get_client_kmer_counter_opt_string(void)
{
    static std::string _s("k:o:hv?");
    return _s;
}

struct option*
kmer_counter::KmerCounter::get_client_kmer_counter_long_options(void)
{
    static struct option _k = { "k",              required_argument,   NULL,    'k' };
    static struct option _o = { "offset",         required_argument,   NULL,    'o' };
    static struct option _h = { "help",           no_argument,         NULL,    'h' };
    static struct option _v = { "version",        no_argument,         NULL,    'v' };
    static struct option _0 = { NULL,             no_argument,         NULL,     0  };
    static std::vector<struct option> _s;
    _s.push_back(_k);
    _s.push_back(_o);
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
                                 this->get_client_kmer_counter_opt_string().c_str(),
                                 this->get_client_kmer_counter_long_options(),
                                 &client_long_index);
    int k = -1;
    int offset = -1;

    opterr = 0; /* disable error reporting by GNU getopt */
    
    while (client_opt != -1) {
        switch (client_opt) {
        case 'k':
            std::sscanf(optarg, "%d", &k);
            this->set_k(k);
            break;
        case 'o':
            std::sscanf(optarg, "%d", &offset);
            this->set_offset(offset);
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
                                 this->get_client_kmer_counter_opt_string().c_str(),
                                 this->get_client_kmer_counter_long_options(),
                                 &client_long_index);
    }
    
    if (optind < argc) {
        do {
            if (this->get_input_fn().empty()) {
                this->set_input_fn(argv[optind]);
                this->initialize_in_stream();
            }
            else {
                std::fprintf(stderr, "Warning: Ignoring additional input file [%s]\n", argv[optind]);
            }
        }
        while (++optind < argc);
    }

    if (this->get_k() == -1) {
        std::fprintf(stderr, "Error: Specify k value\n");
        this->print_usage(stderr);
        std::exit(ENODATA);
    }

    if (this->get_offset() == -1) {
        std::fprintf(stderr, "Error: Specify offset value\n");
        this->print_usage(stderr);
        std::exit(ENODATA);
    }
}

std::string
kmer_counter::KmerCounter::get_client_kmer_counter_name(void)
{
    static std::string _s(kmer_counter::KmerCounter::client_name);
    return _s;
}

std::string
kmer_counter::KmerCounter::get_client_kmer_counter_version(void)
{
    static std::string _s(kmer_counter::KmerCounter::client_version);
    return _s;
}

std::string
kmer_counter::KmerCounter::get_client_kmer_counter_authors(void)
{
    static std::string _s(kmer_counter::KmerCounter::client_authors);
    return _s;
}

std::string
kmer_counter::KmerCounter::get_client_kmer_counter_usage(void)
{
    static std::string _s("\n"                                          \
                          "  Usage:\n"                                  \
                          "\n"                                          \
                          "  $ kmer_counter [options] input\n");
    return _s;
}

std::string
kmer_counter::KmerCounter::get_client_kmer_counter_description(void)
{
    static std::string _s("  Count kmers in BED file and report counts for\n" \
                          "  each element, from specified offset. Write a\n" \
                          "  key-count file and key table as output.\n");
    return _s;
}

std::string
kmer_counter::KmerCounter::get_client_kmer_counter_io_options(void)
{
    static std::string _s("  General Options:\n\n"              \
                          "  --k=n        K-value for kmer length (integer)\n" \
                          "  --offset=n   Offset (integer)\n");
    return _s;
}

std::string
kmer_counter::KmerCounter::get_client_kmer_counter_general_options(void)
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
                 this->get_client_kmer_counter_name().c_str(),
                 this->get_client_kmer_counter_version().c_str(),
                 this->get_client_kmer_counter_authors().c_str(),
                 this->get_client_kmer_counter_usage().c_str(),
                 this->get_client_kmer_counter_description().c_str(),
                 this->get_client_kmer_counter_io_options().c_str(),
                 this->get_client_kmer_counter_general_options().c_str());
}

void
kmer_counter::KmerCounter::print_version(FILE* wo_stream)
{
    std::fprintf(wo_stream,
                 "%s\n"                                              \
                 "  version: %s\n"                                   \
                 "  author:  %s\n",
                 this->get_client_kmer_counter_name().c_str(),
                 this->get_client_kmer_counter_version().c_str(),
                 this->get_client_kmer_counter_authors().c_str());
}
