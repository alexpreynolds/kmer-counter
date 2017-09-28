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
    if (kc.map_keys)
        kc.print_kmer_map(kc.results_kmer_map_stream());
    switch (kc.input_type) {
        case kmer_counter::KmerCounter::bedInput:
            kc.parse_bed_input_to_counts();
            break;
        case kmer_counter::KmerCounter::fastaInput:
            kc.parse_fasta_input_to_counts();
            break;
        default:
            std::fprintf(stderr, "Undefined input type!\n");
            exit(EXIT_FAILURE);
    }
    kc.close_output_streams();
    
    return EXIT_SUCCESS;
}

void
kmer_counter::KmerCounter::close_output_streams(void)
{
    this->close_in_stream();
    this->close_kmer_count_stream();
    this->close_kmer_map_stream();
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
        this->set_mer_key(mer, this->offset(true));
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
kmer_counter::KmerCounter::parse_bed_input_to_counts(void)
{
    char* buf = NULL;
    size_t buf_len = 0;
    ssize_t buf_read = 0;
    char chr_str[LINE_MAX];
    char start_str[LINE_MAX];
    char stop_str[LINE_MAX];
    char* id_str = NULL;
    std::string n("N");

    id_str = (char*) malloc(KMER_COUNTER_LINE_MAX);
    if (!id_str) {
        std::fprintf(stderr, "Error: Could not allocate memory for ID buffer\n");
        std::exit(ENOMEM);
    }

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
            if ((mer_count(mer_f) == 0) && (mer_count(mer_r) == 0)) {
                set_mer_count(mer_f, 1);
                // we don't want to add a palindrome twice
                if (mer_f.compare(mer_r) == 0) {
                    continue;
                }
                set_mer_count(mer_r, 1);
            }
            else if ((mer_count(mer_f) == 1) || (mer_count(mer_r) == 1)) {
                increment_mer_count(mer_f);
                increment_mer_count(mer_r);
            }
        }
        this->print_kmer_count(this->results_kmer_count_stream(), chr_str, start_str, stop_str);
    }

    // cleanup
    free(buf);
    free(id_str);
}

void
kmer_counter::KmerCounter::parse_fasta_input_to_counts(void)
{
    char* buf = NULL;
    size_t buf_len = 0;
    ssize_t buf_read = 0;
    char header_str[LINE_MAX] = {0};
    char* sequence_str = NULL;
    ssize_t sequence_read = 0;
    char* sequence_intermediate_buf = NULL;
    ssize_t line_count = 0;

    sequence_str = (char*) malloc(KMER_COUNTER_LINE_MAX);
    if (!sequence_str) {
        std::fprintf(stderr, "Error: Could not allocate memory for sequence buffer\n");
        std::exit(ENOMEM);
    }

    sequence_intermediate_buf = (char*) malloc(KMER_COUNTER_LINE_MAX);
    if (!sequence_intermediate_buf) {
        std::fprintf(stderr, "Error: Could not allocate memory for sequence intermediate buffer\n");
        std::exit(ENOMEM);
    }

    while ((buf_read = getline(&buf, &buf_len, this->in_stream())) != EOF) {
        if (buf[0] == '>') {
            if ((strlen(header_str) > 0) && (strlen(sequence_str) > 0)) {
                this->process_fasta_record(header_str, sequence_str);
                sequence_read = 0;
            }
            // read in next header
            std::sscanf(buf, ">%s\n", header_str);
        }
        else {
            std::sscanf(buf, "%s\n", sequence_intermediate_buf);
            std::memcpy(sequence_str + sequence_read, sequence_intermediate_buf, buf_read);
            sequence_read += (buf_read - 1);
        }
        line_count += 1;
    }

    // process final record
    if ((strlen(header_str) > 0) && (strlen(sequence_str) > 0)) {
        this->process_fasta_record(header_str, sequence_str);
        sequence_read = 0;
    }

    // cleanup
    free(buf);
    free(sequence_str);
    free(sequence_intermediate_buf);
}

void
kmer_counter::KmerCounter::process_fasta_record(char* header, char* sequence)
{
    std::string seq(sequence);
    std::string n("N");

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
        if ((mer_count(mer_f) == 0) && (mer_count(mer_r) == 0)) {
            set_mer_count(mer_f, 1);
            // we don't want to add a palindrome twice
            if (mer_f.compare(mer_r) == 0) {
                continue;
            }
            set_mer_count(mer_r, 1);
        }
        else if ((mer_count(mer_f) == 1) || (mer_count(mer_r) == 1)) {
            increment_mer_count(mer_f);
            increment_mer_count(mer_r);
        }
    }

    this->print_kmer_count(this->results_kmer_count_stream(), header);
}

void
kmer_counter::KmerCounter::print_kmer_count(FILE* os, char header[])
{
    char kv_pair[LINE_MAX];
    std::string kv_pairs;

    if (!os)
        os = stdout;

    // filter, if we do not want to print reverse complement hits
    std::vector<std::string> mer_keys;
    std::vector<std::string> mer_keys_to_remove;
    for (auto iter = this->mer_counts().begin(); iter != this->mer_counts().end(); ++iter) {
        std::string test_key = iter->first;
        std::string rc_mer_key(test_key);
        reverse_complement_string(rc_mer_key);
        if ((!this->write_reverse_complement) && (test_key.compare(rc_mer_key) != 0)) {
            auto found = std::find(mer_keys.begin(), mer_keys.end(), rc_mer_key);
            if (found == mer_keys.end()) {
                mer_keys_to_remove.push_back(rc_mer_key);
            }
        }
        mer_keys.push_back(test_key);
    }
    for (auto iter = mer_keys_to_remove.begin(); iter != mer_keys_to_remove.end(); ++iter) {
        std::string erase_key = *iter;
        this->erase_mer_count(erase_key);
    }

    // write the hits
    kv_pairs.clear();
    for (auto iter = mer_keys.begin(); iter != mer_keys.end(); ++iter) {
        std::string mer_key = *iter;
        auto mer_key_lookup = this->mer_counts().find(mer_key);
        int mer_key_count = (mer_key_lookup == this->mer_counts().end()) ? 0 : mer_key_lookup->second;
        if (mer_key_count != 0) {
            if (this->map_keys)
                std::sprintf(kv_pair, "%d:%d ", this->mer_key(mer_key), mer_key_count);
            else 
                std::sprintf(kv_pair, "%s:%d ", mer_key.c_str(), mer_key_count);
            set_mer_count(mer_key, 0);
            kv_pairs.append(kv_pair);
        }
    }
    if (kv_pairs.length() > 0) {
        kv_pairs.pop_back();
    }

    std::fprintf(os, ">%s\t%s\n", header, kv_pairs.c_str());    
}

void
kmer_counter::KmerCounter::print_kmer_count(FILE* os, char chr[], char start[], char stop[])
{
    char kv_pair[LINE_MAX];
    std::string kv_pairs;

    if (!os)
        os = stdout;

    // filter, if we do not want to print reverse complement hits
    std::vector<std::string> mer_keys;
    std::vector<std::string> mer_keys_to_remove;
    for (auto iter = this->mer_counts().begin(); iter != this->mer_counts().end(); ++iter) {
        std::string test_key = iter->first;
        std::string rc_mer_key(test_key);
        reverse_complement_string(rc_mer_key);
        if ((!this->write_reverse_complement) && (test_key.compare(rc_mer_key) != 0)) {
            auto found = std::find(mer_keys.begin(), mer_keys.end(), rc_mer_key);
            if (found == mer_keys.end()) {
                mer_keys_to_remove.push_back(rc_mer_key);
            }
        }
        mer_keys.push_back(test_key);
    }
    for (auto iter = mer_keys_to_remove.begin(); iter != mer_keys_to_remove.end(); ++iter) {
        std::string erase_key = *iter;
        this->erase_mer_count(erase_key);
        std::fprintf(stderr, "erased [%s]\n", erase_key.c_str());
    }

    // write the hits
    kv_pairs.clear();
    for (auto iter = mer_keys.begin(); iter != mer_keys.end(); ++iter) {
        std::string mer_key = *iter;
        auto mer_key_lookup = this->mer_counts().find(mer_key);
        int mer_key_count = (mer_key_lookup == this->mer_counts().end()) ? 0 : mer_key_lookup->second;
        if (mer_key_count != 0) {
            if (this->map_keys)
                std::sprintf(kv_pair, "%d:%d ", this->mer_key(mer_key), mer_key_count);
            else
                std::sprintf(kv_pair, "%s:%d ", mer_key.c_str(), mer_key_count);
            set_mer_count(mer_key, 0);
            kv_pairs.append(kv_pair);
        }
    }
    if (kv_pairs.length() > 0) {
        kv_pairs.pop_back();
    }

    std::fprintf(os, "%s\t%s\t%s\t%s\n", chr, start, stop, kv_pairs.c_str());    
}

std::string
kmer_counter::KmerCounter::client_kmer_counter_opt_string(void)
{
    static std::string _s("k:o:r:bfchv?");
    return _s;
}

struct option*
kmer_counter::KmerCounter::client_kmer_counter_long_options(void)
{
    static struct option _k = { "k",              required_argument,   NULL,    'k' };
    static struct option _o = { "offset",         required_argument,   NULL,    'o' };
    static struct option _r = { "results-dir",    required_argument,   NULL,    'r' };
    static struct option _b = { "bed",            no_argument,         NULL,    'b' };
    static struct option _f = { "fasta",          no_argument,         NULL,    'f' };
    static struct option _c = { "no-rc",          no_argument,         NULL,    'c' };
    static struct option _h = { "help",           no_argument,         NULL,    'h' };
    static struct option _v = { "version",        no_argument,         NULL,    'v' };
    static struct option _0 = { NULL,             no_argument,         NULL,     0  };
    static std::vector<struct option> _s;
    _s.push_back(_k);
    _s.push_back(_o);
    _s.push_back(_r);
    _s.push_back(_b);
    _s.push_back(_f);
    _s.push_back(_c);
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

    this->input_type = KmerCounter::undefinedInput;
    this->write_results_to_stdout = false;
    this->write_reverse_complement = true;

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
        case 'b':
            this->input_type = KmerCounter::bedInput;
            break;
        case 'f':
            this->input_type = KmerCounter::fastaInput;
            break;
        case 'c':
            this->write_reverse_complement = false;
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

    this->map_keys = true;
    if (this->offset() == -1) {
        this->map_keys = false;
    }

    if (this->input_type == KmerCounter::undefinedInput) {
        std::fprintf(stderr, "Error: Specify input type value (BED or FASTA)\n");
        this->print_usage(stderr);
        std::exit(ENODATA);
    }

    if (this->results_dir().empty()) {
        this->write_results_to_stdout = true;    
    }
    else {
        this->results_dir_mode(0755);
        if (this->initialize_result_dir(this->results_dir(), this->results_dir_mode())) {
            switch (this->input_type) {
                case kmer_counter::KmerCounter::bedInput:
                    if (!this->write_results_to_stdout)
                        this->initialize_kmer_count_stream("count.bed");
                    break;
                case kmer_counter::KmerCounter::fastaInput:
                    if (!this->write_results_to_stdout)
                        this->initialize_kmer_count_stream("count.txt");
                    break;
                default:
                    std::fprintf(stderr, "Undefined input type!\n");
                    exit(EXIT_FAILURE);
            }
            if (this->map_keys)
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
                          "  --no-rc           Disable writing of non-palindrome reverse complement counts\n" \
                          "  --k=n             K-value for kmer length (integer)\n" \
                          "  --offset=n        Offset (integer)\n" \
                          "  --results-dir=s   Results directory (string)\n");
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
