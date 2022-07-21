//
// Created by henrique on 02/04/19.
//

#ifndef GRAMMARINDCTION_GRAMMAR_H
#define GRAMMARINDCTION_GRAMMAR_H


#include "Rule.h"
#include <random>
#include <unordered_map>
#include <vector>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/access.hpp>
#include <stack>
#include <set>
#include <map>
namespace Grammar {
class Grammar
{

public:
    //Grammar();
    enum training_method {n_gram_count, pfa_baum_welch, pfa_collapsed_gibbs_sample, pfa_alergia, pcfg_inside_outside, pcfg_metropolis_hastings, pcfg_gibbs_sampling, pcsg_metropolis_hastings, pcsg_gibbs_sampling, pcfg_pumping_inference};
    enum grammar_type {n_gram, pfa, pcfg, pcsg};
    const double ALFA = 0.1;
    std::vector<Symbol::Symbol> terminals;
    std::vector<Symbol::Symbol> non_terminals;
    std::vector<Rule::Rule> rules;
    std::vector<std::pair<std::string, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>>> parse_trees;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>>> parse_trees_vec;
    Symbol::Symbol start;
    size_t n_terminals{};
    size_t n_non_terminals{};
    std::pair<size_t, size_t> context_amount;
    std::pair<unsigned long, unsigned long> context_size;
    std::vector<Symbol::Symbol> actual_production;
    std::vector<std::vector<Symbol::Symbol>> words;
    std::random_device rd;
    grammar_type g_tp;

    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive & ar, const unsigned int version)
            {
                //ar & ALFA;
                ar & terminals;
                ar & non_terminals;
                ar & rules;
                ar & parse_trees;
                ar & parse_trees_vec;
                ar & start;
                ar & n_terminals;
                ar & n_terminals;
                ar & context_amount;
                ar & context_size;
                ar & actual_production;
                ar & words;
                //ar & rd;
                ar & g_tp;
            }

    Grammar(const std::vector<Symbol::Symbol> &terminals, int nNonTerminals, std::vector<std::vector<Symbol::Symbol>> words, enum grammar_type g_tp, std::pair<int,int> contextSize);
    Grammar(){};
    void print_grammar();
    [[maybe_unused]] std::string grammar_to_str();
    void print_rules();

    [[maybe_unused]] std::string rules_to_string();
    void train(training_method algorithm, int iterations, double alpha_alergia, double p_ratio, double time_limite);
    std::pair<double,double> perplexity(const std::vector<std::vector<Symbol::Symbol>>& test_data);
    std::pair<double,double> perplexity_kl(const std::vector<std::vector<Symbol::Symbol>>& test_data);
    void prob_sequitur();
    void convert_to_cnf();
    void remove_unused_rules();
    void remove_unused_rule_zero_righties(std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> &right);
    void group_equal_rhs(std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> & right);
    int verify_rule_existence_cnf(std::vector<Symbol::Symbol> rhs_cnf, std::vector<Rule::Rule> rules_nt_cnf);
    bool equal_rhs(std::vector<Symbol::Symbol> rh1, std::vector<Symbol::Symbol> rh2);
    void count_pumping_str(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, int sub_amount, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word);
    void count_pumping_str_by_slice(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<int>> &map_pump_to_word, int w_index, double reduction_share);
    void count_pumping_size1(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word);
    void count_pumping_size2(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word);
    void count_pumping_size3(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word);
    void count_pumping_size4(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word);
    void count_pumping_size5(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word);
    Symbol::Symbol verify_pumping_times(Symbol::Symbol s, int pumping_size);
    std::vector<std::string> convert_symbol_to_vector_string(Symbol::Symbol s);
    std::string convert_vector_to_string(std::vector<Symbol::Symbol> vs);
    void pumping_inference (std::unordered_map<std::string, int> & map, std::unordered_map<std::string, std::vector<int>> &map_pump_to_word);
    static bool equal_word(std::vector<Symbol::Symbol> w1, std::vector<Symbol::Symbol> w2);
    void group_equal_initial_rules();
    void create_new_nt_rule_substring(std::string sb, std::unordered_map<std::string, Symbol::Symbol> & non_terminal_string_map);
    double calculate_prob_word(std::vector<Symbol::Symbol> word);
    double calculate_parse_tree_prob_top_down(std::vector<Symbol::Symbol> word);
    void fpta_pumping_inference(std::unordered_map<std::string, int> & map, std::unordered_map<std::string, std::vector<int>> &map_pump_to_word, int max_ntv, int max_ntw, int max_ntx, int max_ntz);
    void build_node_to_pump_map(std::unordered_map<std::string, std::vector<int>> &map_pump_to_word, std::unordered_map<int, std::vector<std::string>> & node_pump_map, int max_ntv, int max_ntw, int max_ntx, int max_ntz);
    std::set<std::string> build_nt_pumping_set(Symbol::Symbol nt, std::unordered_map<std::string, std::vector<int>> &map_pump_to_word, std::string u, int max_ntv, int max_ntw, int max_ntx, int max_ntz);
    std::unordered_map<std::string, std::set<std::string>> build_nt_pumping_sets(std::unordered_map<std::string, int> &map, int max_ntv, int max_ntw, int max_ntx, int max_ntz);
    bool fpta_pumping_compatible(std::set<std::string> snt1, std::set<std::string> snt2, std::string u1, std::string z1, std::string u2, double pump_percentage_tolerance);
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> build_pumping_righties(std::set<std::string> nt1_pumping_set, size_t index_nt2, size_t index_nt1);
    void find_best_pumping_coverage(std::vector<Symbol::Symbol> word);
    bool pumping_merge_nt(Symbol::Symbol nt1, Symbol::Symbol nt2, Grammar &g, std::vector<Symbol::Symbol> &pumped_nts);
    bool pumping_merge_nt_2(Symbol::Symbol nt1, Symbol::Symbol nt2, Grammar &g, std::vector<Symbol::Symbol> &pumped_nts);
    bool pumping_merge_nt_3(Symbol::Symbol nt1, Symbol::Symbol nt2, Grammar &g, std::vector<Symbol::Symbol> &pumped_nts);
    bool verify_nt_in_subtree_of_nt2(Symbol::Symbol nt1, Symbol::Symbol nt2);
    bool verify_nt_in_subtree_if_any_pumped(Symbol::Symbol nt1, std::vector<Symbol::Symbol> pumped_nts);
    bool verify_nt_in_subtree_if_any_pumped_not_equal(Symbol::Symbol nt1, std::vector<Symbol::Symbol> pumped_nts);
    void add_remaining_subtrees_to_grammar(Symbol::Symbol nt, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> blocked_rhs, Grammar &g, std::vector<Symbol::Symbol> &pumped_nts);
    size_t add_subtrees_of_nt1_to_grammar(Symbol::Symbol nt1, Symbol::Symbol nt2, std::set<std::string> &nt1_pumping_set, Grammar &g, std::vector<Symbol::Symbol> &pumped_nts, size_t index_nt1);
    bool fpta_pumping_compatible_tree(Symbol::Symbol nt1, Symbol::Symbol nt2, double tolerance, std::vector<Rule::Rule> &vector_rules, std::vector<Symbol::Symbol> &vector_symbol);
    bool fpta_pumping_compatible_tree2(Symbol::Symbol nt1, Symbol::Symbol nt2, double tolerance, std::vector<Rule::Rule> &vector_rules, std::vector<Symbol::Symbol> &vector_symbol);
    bool fpta_single_pumping_compatible_tree(Symbol::Symbol nt1, Symbol::Symbol nt2, double tolerance, std::vector<Rule::Rule> &vector_rules, std::vector<Symbol::Symbol> &vector_symbol);
    void fold_fpta_subtree(Symbol::Symbol nt1, Symbol::Symbol nt2, std::vector<Symbol::Symbol> &pumped_nts, std::vector<Rule::Rule> &vector_rules, std::vector<Symbol::Symbol> &vector_symbol);
    std::map<std::pair<std::string,std::string>,bool>  build_compatible_matrix();
    std::map<std::vector<std::pair<int,int>>, std::map<int, double>> find_next_prefix_probabilities(std::vector<Symbol::Symbol> prefix);
    std::vector<std::pair<int, double>>  find_prefix_ranking_probabilities(std::vector<Symbol::Symbol> prefix);
    double find_word_probabilities(std::vector<Symbol::Symbol> word);
    double find_word_probabilities_from_pcfg_inside_table(std::vector<Symbol::Symbol> word);
    bool compatible_alergia(const Symbol::Symbol &a, const Symbol::Symbol &b, double alpha, std::vector<Rule::Rule> & vector_rules );
    void pumping_alergia(double alpha, std::vector<Symbol::Symbol> & pumped_nts);
    void eliminate_covered_pumpings(std::unordered_map<std::string, std::vector<int>> &map_pump_to_word, std::unordered_map<std::string, int> &map);
    void find_pumping_rule(Symbol::Symbol nt1, Symbol::Symbol nt2);
    std::vector<Rule::Rule> find_pumping_rule_by_auto_similarity(Symbol::Symbol nt, std::set<int> &not_search_nts, std::map<int, std::pair<std::vector<std::pair<Symbol::Symbol, int>>, int>> compatible_lists, double p_ratio);
    bool check_auto_similarity(std::vector<Symbol::Symbol> path_to_accept, Symbol::Symbol start);
    bool exist_empty_rule(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right);
    void build_rhs_until_empty_rule(Symbol::Symbol nt1, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> &right, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs_pumping);
    void build_pumping_fronts();
    std::pair<int,int> find_first_compatible_pump();
    void find_v_w_x_z_from_path(std::vector<Symbol::Symbol> path, std::vector<Symbol::Symbol> &v, std::vector<Symbol::Symbol> &w, std::vector<Symbol::Symbol> &x, int v_size, int pumping_size, std::vector<Symbol::Symbol> z_path, std::vector<Symbol::Symbol> &z);
    bool check_derivation_from_nt_with_v_w_x_z(Symbol::Symbol nt, std::vector<Symbol::Symbol> &v, std::vector<Symbol::Symbol> &w, std::vector<Symbol::Symbol> &x, std::vector<Symbol::Symbol> &z, int pumping_times);
    bool check_reachable_node_from_nt_with_v_w_x_z(Symbol::Symbol nt, std::vector<Symbol::Symbol> &v, std::vector<Symbol::Symbol> &w, std::vector<Symbol::Symbol> &x, std::vector<Symbol::Symbol> &z, int pumping_times, int already_pumped_v);
    bool check_v_pumping_use(Symbol::Symbol nt, Symbol::Symbol nt_target, std::vector<Symbol::Symbol> &v, int pumping_times);
    std::vector<Rule::Rule> mount_pumping_rules(std::map<int, std::vector<std::tuple<std::vector<Symbol::Symbol>, std::vector<Symbol::Symbol>, std::vector<Symbol::Symbol>, std::vector<Symbol::Symbol>, std::set<int>, int>>> nts_pumpings, std::set<int> &not_search_nts, Symbol::Symbol nt, int max_height);
    void super_duper_pumping_inference(double alpha, double p_ratio, double time_limite);
    void pump_and_reduce_w(Rule::Rule &r);
    std::map<int, std::pair<std::vector<std::pair<Symbol::Symbol, int>>, int>> build_compatible_lists(double alpha);
    void decrease_rules_from_pumping(Symbol::Symbol nt, std::vector<Rule::Rule> & rs, int max_height, std::vector<Symbol::Symbol> & z);
    void recursive_build_parse_tree_indexes(std::vector<Rule::Rule> & rs, int & max_height, std::vector<Symbol::Symbol> & z, std::vector<std::pair<int,int>> parse_tree_indexes, std::vector<std::vector<std::pair<int,int>>> & parse_trees);
    bool check_reacheable_multiple_pumpings(Symbol::Symbol nt, std::vector<std::pair<int, int>> parse_tree_indexes, std::vector<Rule::Rule> &rs);




private:
    bool check_id_nt_rule();
    void baum_welch(int iteration);
    void inside_outside(int iterations);
    void gen_fpta();
    void alergia(double alpha);
    void collapsed_gibbs_sample_pfa(int iterations);
    void train_n_gram();
    void generate_non_termnals();
    void generate_permutation(std::vector<std::vector<Symbol::Symbol>> & permutations, std::vector<Symbol::Symbol> symbols, size_t size, std::vector<Symbol::Symbol> word, bool context);
    void generate_rules_cnf();
    void generate_rules_regular();
    void generate_n_gram_rules();
    void generate_n_gram_non_terminals();
    void metropolis_hastings_pcfg(int iterations);
    void gibbs_sampling_pcfg(int iterations);
    void gibbs_sampling_pcsg(int iterations);
    void metropolis_hastings_pcsg(int iterations, double time_limit);
    void check_digram_position_integrity(std::unordered_map<std::string, std::tuple<int, int, int>> & digram_position);
    void check_digram_position_integrity_by_rules(std::unordered_map<std::string, std::tuple<int, int, int>> &digram_position);

    bool enforce_digram_uniqueness(std::unordered_map<std::string, int> & digram_map, std::unordered_map<std::string, std::tuple<int, int, int>> & digram_position, int i);
    void enforce_rule_utility(std::unordered_map<std::string, std::tuple<int, int, int>> & digram_position, int ii);
    bool verify_duplicate_digram(std::unordered_map<std::string, std::tuple<int, int, int>> digram_map, std::string digram);
    double **** cyk_prob_kl_vec(std::vector<Symbol::Symbol> w);
    /*double **** CYKProbKLVecOpt(std::vector<Symbol::Symbol> w);*/
    double **** cyk_prob_kln_vec(const std::vector<Symbol::Symbol>& w);
    void print_inside_table(double ***p, int wSize) const;
    /*void printInsideTableKL(double ****p, int wSize);*/
    static void print_outside_table(const std::vector<std::vector<std::vector<double>>>& p);
    void sample_parse_tree(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, Rule::Rule r, const std::string& w, double ***inside_table, size_t i, size_t k);
    void sample_parse_tree_kl(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, Rule::Rule r, const std::string& w, double ****inside_table, size_t i, size_t k);
    void sample_parse_tree_kl_vec(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, Rule::Rule r, std::vector<Symbol::Symbol> w, double ****inside_table, size_t i, size_t k);
    void sample_parse_tree_kl_vec_opt(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, Rule::Rule r, std::vector<Symbol::Symbol> w, double ****inside_table, size_t i, size_t k);
    void calculate_new_theta_vec_opt(int i);
    int calculate_producton_counts (const std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>& production, const std::string& w);
    std::vector<std::pair<double, int>> calculate_rule_frequence(Rule::Rule r, const std::string& w);
    static double c_constant(std::vector<std::pair<double, int>> ruleFrequence);
    static double prob_tree (std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> tree);
    static bool equal_productions(std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> prod_a,
                                  std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> prod_b);
    void update_parse_tress_theta();
    void apply_production(std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> prod_a);
    Rule::Rule find_rule_by_lhs(std::vector<Symbol::Symbol> lhs);
    Symbol::Symbol find_terminal_by_name(std::string name);
    Symbol::Symbol find_symbol_in_vector_by_name(std::string name, std::vector<Symbol::Symbol> & vec);

    int convert_context_to_id(int side, std::vector<Symbol::Symbol> context);
    void get_actual_context(std::vector<Symbol::Symbol> & leftContext, std::vector<Symbol::Symbol> & rightContext);
    /*static void printProduction(
            std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>> tree);
    static void
    printOneProduction(std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> prd);*/
    void sample_parse_tree_vec(
            std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>> & vr,
            Rule::Rule r, std::vector<Symbol::Symbol> w, double ***p_double, size_t i, size_t k);
    double ***cyk_prob_vec(std::vector<Symbol::Symbol> w);
    std::vector<std::vector<std::vector<double>>> outside_table(const std::vector<Symbol::Symbol>& w, double ***inside_table);
    double ***cyk_prob_n_vec(const std::vector<Symbol::Symbol>& w);
    void calculate_new_theta_vec(const std::vector<Symbol::Symbol>& w);
    int calculate_producton_counts_vec(
            const std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>& pair,
            const std::vector<Symbol::Symbol>& w);

    double p_ti_ti_minus_1_vec(const std::vector<Symbol::Symbol>& w);
    double p_ti_ti_minus_1_vec_opt(int i) ;
    std::vector<std::pair<double, int>> calculate_rule_frequence_vec(Rule::Rule &rule, const std::vector<Symbol::Symbol>& w);
    void p_ti_minus_1_frequence(int i);
    void p_ti_minus_1_plus_frequence(int i);
    void free_inside_table(double ***p, size_t wSize) const;
    void free_inside_table_kl(double ****p, size_t wSize) const;
    void baum_welch_expectation(std::vector<double> &count_nd, std::vector<double> &count_nl, std::vector<std::vector<std::vector<double>>> &count_nt);
    void baum_welch_maximization(std::vector<double> &count_nd_r, std::vector<double> &count_nl_r, std::vector<std::vector<std::vector<double>>> &count_ntr);
    void calculate_wv(std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &w, std::vector<std::vector<std::vector<double>>> &v, const std::vector<Symbol::Symbol>& word, double***inside_table, std::vector<std::vector<std::vector<double>>>o_t);

    static bool test_alergia(double probFa, double freqa, double probFb, double freqb, double alpha);
    void stochastic_merge(const Symbol::Symbol& a, const Symbol::Symbol& b);
    void stochastic_pumping_merge(const Symbol::Symbol& a, const Symbol::Symbol& b, std::vector<Symbol::Symbol> & pumped_nts);
    void stochastic_fold(const Symbol::Symbol& a, const Symbol::Symbol& b);
    double stochastic_pumping_fold(const Symbol::Symbol& a, const Symbol::Symbol& b, std::vector<Symbol::Symbol> & pumped_nts);
    void remove_unused_nt();
    void recursive_insert_unused(std::vector<Symbol::Symbol> &unused, const Symbol::Symbol& nt);
    void normalize_probs();
    void sample_regular_rules(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, const std::vector<Symbol::Symbol>& w);
    void recursive_add_terminal_to_nt(Symbol::Symbol nt, size_t n, int &nIds);
    void add_n_gram_rule_frequency(const Symbol::Symbol& lhs, const Symbol::Symbol& next_symbol);
    Symbol::Symbol is_handle(std::vector<Symbol::Symbol> sub_word);
    std::vector<Symbol::Symbol> remove_empty_substring(std::vector<Symbol::Symbol> word);
    void nullify_unreacheble_rules();
    bool build_words_from_rule(Symbol::Symbol nt, std::vector<std::pair<int, int>> parse_tree_indexes, std::vector<Rule::Rule> &rs);


public:
    virtual ~Grammar();

    bool compare_pumping_use(std::vector<Symbol::Symbol> vector, std::vector<int> & pumping_string);
    void find_pumping_rules(int &v, std::string uvxyz, std::unordered_map<std::string, Symbol::Symbol> & non_terminal_string_map);
    void calculate_new_rule_from_starting_symbol(std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> &new_start_right, int v, std::string uvxyz, std::unordered_map<std::string, Symbol::Symbol> &non_terminal_string_map);
    std::vector<Symbol::Symbol> yield_string(std::vector<Symbol::Symbol> vector);
    int load_map_pump_to_word(std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> map_pump_to_word);
    std::vector<std::vector<Symbol::Symbol>> generate_max_size_words_from_rules (int max_t);
    void generate_nt_for_t ();
};
}

#endif //GRAMMARINDCTION_GRAMMAR_H
