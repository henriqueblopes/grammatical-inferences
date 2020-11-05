//
// Created by henrique on 02/04/19.
//

#ifndef GRAMMARINDCTION_GRAMMAR_H
#define GRAMMARINDCTION_GRAMMAR_H



#include "Rule.h"
#include <vector>
#include <random>
namespace Grammar {
class Grammar
{

public:
    //Grammar();
    enum training_method {n_gram_count, pfa_baum_welch, pfa_collapsed_gibbs_sample, pfa_alergia, pcfg_inside_outside, pcfg_metropolis_hastings, pcfg_gibbs_sampling, pcsg_metropolis_hastings, pcsg_gibbs_sampling};
    enum grammar_type {n_gram, pfa, pcfg, pcsg};
    const double ALFA = 0.1;
    std::vector<Symbol::Symbol> terminals;
    std::vector<Symbol::Symbol> non_terminals;
    std::vector<Rule::Rule> rules;
    std::vector<std::pair<std::string, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>>> parse_trees;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>>> parse_trees_vec;
    Symbol::Symbol start;
    int n_terminals{};
    int n_non_terminals{};
    std::pair<int, int> context_amount;
    std::pair<unsigned long, unsigned long> context_size;
    std::vector<Symbol::Symbol> actual_production;
    std::vector<std::vector<Symbol::Symbol>> words;
    std::random_device rd;
    grammar_type g_tp;

    Grammar(const std::vector<Symbol::Symbol> &terminals, int nNonTerminals, std::vector<std::vector<Symbol::Symbol>> words, enum grammar_type g_tp, std::pair<int,int> contextSize);

    void print_grammar();
    [[maybe_unused]] std::string grammar_to_str();
    void print_rules();

    [[maybe_unused]] std::string rules_to_string();
    void train(training_method algorithm, int iterations);
    std::pair<double,double> perplexity(const std::vector<std::vector<Symbol::Symbol>>& test_data);
    std::pair<double,double> perplexity_kl(const std::vector<std::vector<Symbol::Symbol>>& test_data);

private:
    void baum_welch(int iteration);
    void inside_outside(int iterations);
    void gen_fpta();
    void alergia(double alpha);
    void collapsed_gibbs_sample_pfa(int iterations);
    void train_n_gram();
    void generate_non_termnals();
    void generate_permutation(std::vector<std::vector<Symbol::Symbol>> & permutations, std::vector<Symbol::Symbol> symbols, unsigned long size, std::vector<Symbol::Symbol> word, bool context);
    void generate_rules_cnf();
    void generate_rules_regular();
    void generate_n_gram_rules();
    void generate_n_gram_non_terminals();
    void metropolis_hastings_pcfg(int iterations);
    void gibbs_sampling_pcfg(int iterations);
    void gibbs_sampling_pcsg(int iterations);
    void metropolis_hastings_pcsg(int iterations);
    double **** cyk_prob_kl_vec(std::vector<Symbol::Symbol> w);
    /*double **** CYKProbKLVecOpt(std::vector<Symbol::Symbol> w);*/
    double **** cyk_prob_kln_vec(const std::vector<Symbol::Symbol>& w);
    void print_inside_table(double ***p, int wSize) const;
    /*void printInsideTableKL(double ****p, int wSize);*/
    static void print_outside_table(const std::vector<std::vector<std::vector<double>>>& p);
    void sample_parse_tree(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, Rule::Rule r, const std::string& w, double ***inside_table, unsigned int i, unsigned int k);
    void sample_parse_tree_kl(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, Rule::Rule r, const std::string& w, double ****inside_table, unsigned int i, unsigned int k);
    void sample_parse_tree_kl_vec(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, Rule::Rule r, std::vector<Symbol::Symbol> w, double ****inside_table, unsigned int i, unsigned int k);
    void sample_parse_tree_kl_vec_opt(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, Rule::Rule r, std::vector<Symbol::Symbol> w, double ****inside_table, unsigned int i, unsigned int k);
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
    [[nodiscard]] int convert_context_to_id(int side, std::vector<Symbol::Symbol> context) const;
    void get_actual_context(std::vector<Symbol::Symbol> & leftContext, std::vector<Symbol::Symbol> & rightContext);
    /*static void printProduction(
            std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>> tree);
    static void
    printOneProduction(std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> prd);*/
    void sample_parse_tree_vec(
            std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>> & vr,
            Rule::Rule r, std::vector<Symbol::Symbol> w, double ***p_double, unsigned int i, unsigned int k);
    double ***cyk_prob_vec(std::vector<Symbol::Symbol> w);
    std::vector<std::vector<std::vector<double>>> outside_table(const std::vector<Symbol::Symbol>& w, double ***inside_table);
    double ***cyk_prob_n_vec(const std::vector<Symbol::Symbol>& w);
    void calculate_new_theta_vec(const std::vector<Symbol::Symbol>& w);
    int calculate_producton_counts_vec(
            const std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>& pair,
            const std::vector<Symbol::Symbol>& w);
    static bool equal_word(std::vector<Symbol::Symbol> w1, std::vector<Symbol::Symbol> w2);
    double p_ti_ti_minus_1_vec(const std::vector<Symbol::Symbol>& w);
    double p_ti_ti_minus_1_vec_opt(int i) ;
    std::vector<std::pair<double, int>> calculate_rule_frequence_vec(Rule::Rule &rule, const std::vector<Symbol::Symbol>& w);
    void p_ti_minus_1_frequence(int i);
    void p_ti_minus_1_plus_frequence(int i);
    void free_inside_table(double ***p, int wSize) const;
    void free_inside_table_kl(double ****p, int wSize) const;
    void baum_welch_expectation(std::vector<double> &count_nd, std::vector<double> &count_nl, std::vector<std::vector<std::vector<double>>> &count_nt);
    void baum_welch_maximization(std::vector<double> &count_nd_r, std::vector<double> &count_nl_r, std::vector<std::vector<std::vector<double>>> &count_ntr);
    void calculate_wv(std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &w, std::vector<std::vector<std::vector<double>>> &v, const std::vector<Symbol::Symbol>& word, double***inside_table, std::vector<std::vector<std::vector<double>>>o_t);
    bool compatible_alergia(const Symbol::Symbol& a, const Symbol::Symbol& b, double alpha);
    static bool test_alergia(double probFa, double freqa, double probFb, double freqb, double alpha);
    void stochastic_merge(const Symbol::Symbol& a, const Symbol::Symbol& b);
    void stochastic_fold(const Symbol::Symbol& a, const Symbol::Symbol& b);
    void remove_unused_nt();
    void recursive_insert_unused(std::vector<Symbol::Symbol> &unused, const Symbol::Symbol& nt);
    void normalize_probs();
    void sample_regular_rules(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, const std::vector<Symbol::Symbol>& w);
    void recursive_add_terminal_to_nt(Symbol::Symbol nt, int n, int &nIds);
    void add_n_gram_rule_frequency(const Symbol::Symbol& lhs, const Symbol::Symbol& next_symbol);

public:
    virtual ~Grammar();
};
}

#endif //GRAMMARINDCTION_GRAMMAR_H
