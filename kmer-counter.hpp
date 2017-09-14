#ifndef KMER_COUNTER_H_
#define KMER_COUNTER_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif /* getline() support */

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#include <string>
#include <vector>
#include <deque>
#include <algorithm>
#include <iostream>
#include <exception>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cinttypes>
#include <cstring>
#include <getopt.h>
#include <pthread.h>
#include <sys/stat.h>
#include "hash_map.hpp"

#define KMER_COUNTER_LINE_MAX 524280

namespace kmer_counter
{
    class KmerCounter
    {
        
    private:
        int _k;
        int _offset;
        std::string _input_fn;
        FILE* _in_stream;
        std::string _results_dir;
        std::string _results_kmer_count_fn;
        FILE* _results_kmer_count_stream = NULL;
        std::string _results_kmer_map_fn;
        FILE* _results_kmer_map_stream = NULL;
        mode_t _results_dir_mode;
        emilib::HashMap<std::string, int> _mer_keys;
        emilib::HashMap<std::string, int> _mer_counts;        
        
    public:
        enum KmerCounterInput {
            undefinedInput = 0,
            bedInput,
            fastaInput
        };

        void parse_bed_input_to_counts(void);
        void parse_fasta_input_to_counts(void);
        void initialize_command_line_options(int argc, char** argv);
        void initialize_kmer_map(void);
        void print_kmer_map(FILE* wo_stream);
        void print_kmer_count(FILE* os, char header[]);
        void print_kmer_count(FILE* wo_stream, char chr[], char start[], char stop[]);
        void close_output_streams(void);
        // --
        static const std::string client_name;
        static const std::string client_version;
        static const std::string client_authors;
        KmerCounterInput input_type;
        bool map_keys;
        bool write_results_to_stdout;
        std::string client_kmer_counter_opt_string(void);
        struct option* client_kmer_counter_long_options(void);
        std::string client_kmer_counter_name(void);
        std::string client_kmer_counter_version(void);
        std::string client_kmer_counter_authors(void);
        std::string client_kmer_counter_usage(void);
        std::string client_kmer_counter_description(void);
        std::string client_kmer_counter_io_options(void);
        std::string client_kmer_counter_general_options(void);
        void print_usage(FILE* wo_stream);
        void print_version(FILE* wo_stream);
        const std::string& results_dir(void);
        void results_dir(const std::string& s);
        const std::string& results_kmer_count_fn(void);
        void results_kmer_count_fn(const std::string& s);
        FILE* results_kmer_count_stream(void);
        void results_kmer_count_stream(FILE **wo_stream_ptr);
        const std::string& results_kmer_map_fn(void);
        void results_kmer_map_fn(const std::string& s);
        FILE* results_kmer_map_stream(void);
        void results_kmer_map_stream(FILE **wo_stream_ptr);
        const mode_t& results_dir_mode(void);
        void results_dir_mode(const mode_t& m);
        bool initialize_result_dir(const std::string& s, const mode_t m);
        void initialize_kmer_count_stream(const std::string& fn);
        void close_kmer_count_stream(void);
        void initialize_kmer_map_stream(void);
        void close_kmer_map_stream(void);
        FILE* in_stream(void);
        void in_stream(FILE** ri_stream_ptr);
        void initialize_in_stream(void);
        void close_in_stream(void);
        const std::string& input_fn(void);
        void input_fn(const std::string& s);
        const int& k(void);
        void k(const int& k);
        const int& offset(void);
        int offset(const bool& increment);
        void offset(const int& o);
        void increment_offset(void);
        const emilib::HashMap<std::string, int>& mer_counts(void);
        void mer_counts(const emilib::HashMap<std::string, int>& mc);
        auto mer_count(const std::string& k);
        void set_mer_count(const std::string& k, const int& v);
        void increment_mer_count(const std::string& k);        
        const emilib::HashMap<std::string, int>& mer_keys(void);
        void mer_keys(const emilib::HashMap<std::string, int>& mk);
        void set_mer_key(const std::string& k, const int& v);
        auto mer_key(const std::string& k);
        
        static void reverse_complement_string(std::string &s) {
            std::reverse(s.begin(), s.end());
            for (auto i = s.begin(); i != s.end(); ++i) {
                switch(*i) {
                case 'A':
                    *i = 'T';
                    break;
                case 'G':
                    *i = 'C';
                    break;
                case 'C':
                    *i = 'G';
                    break;
                case 'T':
                    *i = 'A';
                    break;
                }
            }
        }

        static int do_mkdir(const char *path, mode_t mode) {
            struct stat st;
            int status = 0;
            
            if (stat(path, &st) != 0) {
                /* Directory does not exist. EEXIST for race condition */
                if (mkdir(path, mode) != 0 && errno != EEXIST)
                    status = -1;
            }
            else if (!S_ISDIR(st.st_mode)) {
                errno = ENOTDIR;
                status = -1;
            }
            
            return status;
        }
        
        static int mkpath(const char *path, mode_t mode) {
            char* pp = NULL;
            char* sp = NULL;
            int status = 0;
            char* copypath = NULL;

            size_t path_len = strlen(path) + 1;
            copypath = (char*) malloc(path_len);
            if (!copypath) {
                std::fprintf(stderr, "Error: Could not allocate memory for copypath\n");
                std::exit(ENOMEM);
            }
            std::memcpy(copypath, path, path_len);
            
            status = 0;
            pp = copypath;
            while (status == 0 && (sp = strchr(pp, '/')) != 0) {
                if (sp != pp) {
                    /* Neither root nor double slash in path */
                    *sp = '\0';
                    status = do_mkdir(copypath, mode);
                    *sp = '/';
                }
                pp = sp + 1;
            }
            if (status == 0) {
                status = do_mkdir(path, mode);
            }
            if (status == -1) {
                char msg_buf[LINE_MAX] = {0};
                sprintf(msg_buf, "Error: Could not create specified path [%s] with specified mode [%o]\n", path, mode);
                throw std::runtime_error(msg_buf);
            }
            free(copypath);
            
            return status;
        }

        KmerCounter();
        ~KmerCounter();
    };

    const emilib::HashMap<std::string, int>& KmerCounter::mer_counts(void) { return _mer_counts; }
    auto KmerCounter::mer_count(const std::string& k) { return _mer_counts.count(k); }
    void KmerCounter::mer_counts(const emilib::HashMap<std::string, int>& mc) { _mer_counts = mc; }
    void KmerCounter::set_mer_count(const std::string& k, const int& v) { _mer_counts[k] = v; }
    void KmerCounter::increment_mer_count(const std::string& k) { _mer_counts[k]++; }

    const emilib::HashMap<std::string, int>& KmerCounter::mer_keys(void) { return _mer_keys; }
    void KmerCounter::mer_keys(const emilib::HashMap<std::string, int>& mk) { _mer_keys = mk; }
    void KmerCounter::set_mer_key(const std::string& k, const int& v) { _mer_keys[k] = v; }
    auto KmerCounter::mer_key(const std::string& k) { return _mer_keys[k]; }
    
    const std::string& KmerCounter::results_dir(void) { return _results_dir; }
    void KmerCounter::results_dir(const std::string& s) { _results_dir = s; }

    const std::string& KmerCounter::results_kmer_count_fn(void) { return _results_kmer_count_fn; }
    void KmerCounter::results_kmer_count_fn(const std::string& s) { _results_kmer_count_fn = s; }
    FILE* KmerCounter::results_kmer_count_stream(void) { return _results_kmer_count_stream; }
    void KmerCounter::results_kmer_count_stream(FILE **wsp) { _results_kmer_count_stream = *wsp; }
    
    const std::string& KmerCounter::results_kmer_map_fn(void) { return _results_kmer_map_fn; }
    void KmerCounter::results_kmer_map_fn(const std::string& s) { _results_kmer_map_fn = s; }
    FILE* KmerCounter::results_kmer_map_stream(void) { return _results_kmer_map_stream; }
    void KmerCounter::results_kmer_map_stream(FILE **wsp) { _results_kmer_map_stream = *wsp; }    
    
    const mode_t& KmerCounter::results_dir_mode(void) { return _results_dir_mode; }
    void KmerCounter::results_dir_mode(const mode_t& m) { _results_dir_mode = m; }
    
    bool KmerCounter::initialize_result_dir(const std::string& s, const mode_t m) { return (mkpath(s.c_str(), m) != 0) ? false : true; }
    
    void KmerCounter::initialize_kmer_count_stream(const std::string& fn) {
        FILE* out_fp = NULL;
        std::string _kmer_count_fn(this->results_dir() + "/" + fn);
        this->results_kmer_count_fn(_kmer_count_fn);
        out_fp = this->results_kmer_count_fn().empty() ? NULL : std::fopen(this->results_kmer_count_fn().c_str(), "w");
        if (!out_fp) {
            std::fprintf(stderr, "Error: Output file handle to count stream could not be created\n");
            std::exit(ENODATA); /* No message is available on the STREAM head read queue (POSIX.1) */
        }
        this->results_kmer_count_stream(&out_fp);
    }
    void KmerCounter::close_kmer_count_stream(void) { if (_results_kmer_count_stream) { std::fclose(_results_kmer_count_stream); } }

    void KmerCounter::initialize_kmer_map_stream(void) {
        FILE* out_fp = NULL;
        std::string _kmer_map_fn(this->results_dir() + "/" + "map.txt");
        this->results_kmer_map_fn(_kmer_map_fn);
        out_fp = this->results_kmer_map_fn().empty() ? NULL : std::fopen(this->results_kmer_map_fn().c_str(), "w");
        if (!out_fp) {
            std::fprintf(stderr, "Error: Output file handle to map stream could not be created\n");
            std::exit(ENODATA); /* No message is available on the STREAM head read queue (POSIX.1) */
        }
        this->results_kmer_map_stream(&out_fp);
    }
    void KmerCounter::close_kmer_map_stream(void) { if (map_keys && _results_kmer_map_stream) { std::fclose(_results_kmer_map_stream); } }    
    
    const int& KmerCounter::k(void) { return _k; }
    void KmerCounter::k(const int& k) { _k = k; }
    
    const int& KmerCounter::offset(void) { return _offset; }
    void KmerCounter::offset(const int& o) { _offset = o; }
    int KmerCounter::offset(const bool& increment) { int _o = _offset; if (increment) { _offset++; } return _o; }
    void KmerCounter::increment_offset(void) { ++_offset; }

    FILE* KmerCounter::in_stream(void) { return _in_stream; }
    void KmerCounter::in_stream(FILE** isp) { _in_stream = *isp; }
    void KmerCounter::initialize_in_stream(void) {
        FILE* in_fp = NULL;
        in_fp = this->input_fn().empty() ? stdin : std::fopen(this->input_fn().c_str(), "r");
        if (!in_fp) {
            std::fprintf(stderr, "Error: Input file handle could not be created\n");
            std::exit(ENODATA); /* No message is available on the STREAM head read queue (POSIX.1) */
        }
        this->in_stream(&in_fp);
    }
    void KmerCounter::close_in_stream(void) {
        std::fclose(this->in_stream());
    }

    const std::string& KmerCounter::input_fn(void) { return _input_fn; }
    void KmerCounter::input_fn(const std::string& s) {
        struct stat s_stat;
        if (stat(s.c_str(), &s_stat) == 0) {
            _input_fn = s;
        }
        else {
            std::fprintf(stderr, "Error: Input file does not exist (%s)\n", s.c_str());
            std::exit(ENODATA); /* No message is available on the STREAM head read queue (POSIX.1) */
        }
    }

    KmerCounter::KmerCounter() {
        k(-1);
        offset(-1);
    }
    
    KmerCounter::~KmerCounter() {
    }
}

#endif // KMER_COUNTER_H_
