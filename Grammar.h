//
// Created by henrique on 02/04/19.
//

#ifndef GRAMMARINDCTION_GRAMMAR_H
#define GRAMMARINDCTION_GRAMMAR_H



#include "Rule.h"
#include <vector>
class Grammar
{

public:
    const double ALFA = 0.1;
    std::vector<Symbol> terminals;
    std::vector<Symbol> nonterminals;
    std::vector<Rule> rules;
    std::vector<std::pair<std::string, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>>> parseTrees;
    std::vector<std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>>> parseTreesVec;
    Symbol start;
    int maxRightHandSize;
    int maxProductionRules;
    int nTerminals;
    int nNonTerminals;
    std::pair<int, int> contextAmount;
    std::pair<int,int> contextSize;
    std::vector<Symbol> actualProduction;
    std::vector<std::vector<Symbol>> words;

    Grammar(const std::vector<Symbol> &terminals, const std::vector<Symbol> &nonterminals,
            const std::vector<Rule> &rules, const Symbol &start);

    Grammar(const std::vector<Symbol> &terminals, std::pair<int, int> contextSize,
            int nNonTerminals, std::vector<std::vector<Symbol>> words, int trainingM);

    Grammar(const std::vector<Symbol> &terminals, const std::vector<Symbol> &nonterminals,
            const std::vector<Rule> &rules, const Symbol &start, int maxRightHandSize, int maxProductionRules,
            int nTerminals, int nNonTerminals);

    void printGrammar();
    std::string grammarToStr();
    void printRules();
    std::string rulesToString();
    void train(int algorithm, int iterations);
    std::pair<double,double> perplexity(std::vector<std::vector<Symbol>> testData, bool normalized);
    std::pair<double,double> perplexityKL(std::vector<std::vector<Symbol>> testData, bool normalized);

private:
    void generateNonTermnals();
    void generatePermutation(std::vector<std::vector<Symbol>> & permutations, std::vector<Symbol> symbols, int size, std::vector<Symbol> word, bool context);
    void generateRulesCNF();
    void inductGFG();
    void metropolisHastingsPCFG(int iterations);
    void gibbsSamplingPCFG(int iterations);
    void gibbsSamplingPCSG(int iterations);
    void metropolisHastingsPCSG(int iterations);
    double *** CYKProb(std::string w);
    double **** CYKProbKL(std::string w);
    double **** CYKProbKLVec(std::vector<Symbol> w);
    double **** CYKProbKLVecOpt(std::vector<Symbol> w);
    double **** CYKProbKLNVec(std::vector<Symbol> w);
    static inline void free(Grammar& G);
    void printInsideTable(double ***p, int wSize);
    void printInsideTableKL(double ****p, int wSize);
    void printInsideTableKLVec(double ****p, int wSize);
    void sampleParseTree(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, std::string w, double ***insideTable, int i, int k);
    void sampleParseTreeKL(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, std::string w, double ****insideTable, int i, int k);
    void sampleParseTreeKLVec(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, std::vector<Symbol> w, double ****insideTable, int i, int k);
    void sampleParseTreeKLVecOpt(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, std::vector<Symbol> w, double ****insideTable, int i, int k);

    void calculateNewTheta(std::string w);
    void calculateNewThetaVecOpt(std::vector<Symbol> w, int i);
    int calculateProductonCounts (std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production, std::string w);
    std::vector<std::pair<double, int>> calculateRuleFrequence(Rule r, std::string w);
    double cConstant(std::vector<std::pair<double, int>> ruleFrequence);
    double pTiTiMinus1 (std:: string w);
    double probTree (std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> tree);
    bool equalProductions(std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prodA,
                     std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prodB);

    void updateParseTressTheta();
    void applyProduction(std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prodA);
    Rule findRuleByLHS(std::vector<Symbol> lhs);
    void printProduction(
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> tree);

    int convertContextToID(int side, std::vector<Symbol> context);
    void getActualContext(std::vector<Symbol> & leftContext, std::vector<Symbol> & rightContext);

    void
    printOneProduction(std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prd);

    void sampleParseTreeVec(
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> & vr,
            Rule r, std::vector<Symbol> w, double ***pDouble, int i, unsigned long i1);

    double ***CYKProbVec(std::vector<Symbol> w);
    double ***CYKProbNVec(std::vector<Symbol> w);

    void calculateNewThetaVec(std::vector<Symbol> w);

    int calculateProductonCountsVec(
            std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> pair,
            std::vector<Symbol> w);

    bool equalWord(std::vector<Symbol> w1, std::vector<Symbol> w2);

    double pTiTiMinus1Vec(std::vector<Symbol> w);
    double pTiTiMinus1VecOpt(std::vector<Symbol> w, int i) ;

    std::vector<std::pair<double, int>> calculateRuleFrequenceVec(Rule &rule, std::vector<Symbol> w);
    void pTiMinus1Frequence(std::vector<Symbol> w, int i);
    void pTiMinus1PlusFrequence(std::vector<Symbol> w, int i);
    void freeInsideTable(double ***p, int wSize);
    void freeInsideTableKL(double ****p, int wSize);

public:
    virtual ~Grammar();
};


#endif //GRAMMARINDCTION_GRAMMAR_H
