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
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cinttypes>
#include <cstring>
#include <getopt.h>
#include <pthread.h>
#include <sys/stat.h>
#include "hash_map.hpp"

namespace kmer_counter
{
    class KmerCounter
    {
        
    private:
        std::string _input_fn;
        FILE* _in_stream;
        int _k;
        int _offset;
        
    public:
        void parse_input(void);
        void initialize_command_line_options(int argc, char** argv);
        static const std::string client_name;
        static const std::string client_version;
        static const std::string client_authors;
        std::string get_client_kmer_counter_opt_string(void);
        struct option* get_client_kmer_counter_long_options(void);
        std::string get_client_kmer_counter_name(void);
        std::string get_client_kmer_counter_version(void);
        std::string get_client_kmer_counter_authors(void);
        std::string get_client_kmer_counter_usage(void);
        std::string get_client_kmer_counter_description(void);
        std::string get_client_kmer_counter_io_options(void);
        std::string get_client_kmer_counter_general_options(void);
        void print_usage(FILE* wo_stream);
        void print_version(FILE* wo_stream);
        FILE* get_in_stream(void);
        FILE** get_in_stream_ptr(void);
        void set_in_stream(FILE* ri_stream);
        void initialize_in_stream(void);
        void close_in_stream(void);
        std::string get_input_fn(void);
        void set_input_fn(std::string s);
        int get_k(void);
        void set_k(int k);
        int get_offset(void);
        void set_offset(int o);

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

        KmerCounter();
        ~KmerCounter();
    };

    int KmerCounter::get_k(void) { return _k; }
    void KmerCounter::set_k(int k) { _k = k; }
    
    int KmerCounter::get_offset(void) { return _offset; }
    void KmerCounter::set_offset(int o) { _offset = o; }

    FILE* KmerCounter::get_in_stream(void) { return _in_stream; }
    FILE** KmerCounter::get_in_stream_ptr(void) { return &_in_stream; }    
    void KmerCounter::set_in_stream(FILE* is) { _in_stream = is; }
    void KmerCounter::initialize_in_stream(void) {
        FILE* in_fp = NULL;
        in_fp = this->get_input_fn().empty() ? stdin : std::fopen(this->get_input_fn().c_str(), "r");
        if (!in_fp) {
            std::fprintf(stderr, "Error: Input file handle could not be created\n");
            std::exit(ENODATA); /* No message is available on the STREAM head read queue (POSIX.1) */
        }
        this->set_in_stream(in_fp);
    }
    void KmerCounter::close_in_stream(void) {
        std::fclose(this->get_in_stream());
    }

    std::string KmerCounter::get_input_fn(void) { return _input_fn; }
    void KmerCounter::set_input_fn(std::string s) {
        struct stat buf;
        if (stat(s.c_str(), &buf) == 0) {
            _input_fn = s;
        }
        else {
            std::fprintf(stderr, "Error: Input file does not exist (%s)\n", s.c_str());
            std::exit(ENODATA); /* No message is available on the STREAM head read queue (POSIX.1) */
        }
    }
    
    KmerCounter::KmerCounter() {
        set_k(-1);
        set_offset(-1);
    }
    
    KmerCounter::~KmerCounter() {
    }
    
}

#endif // KMER_COUNTER_H_
