//
// Created by henrique on 02/04/19.
//

#ifndef GRAMMARINDCTION_GRAMMAR_H
#define GRAMMARINDCTION_GRAMMAR_H



#include "Rule.h"
#include <vector>
#include <random>

class Grammar
{

public:
    //Grammar();
    const double ALFA = 0.1;
    std::vector<Symbol> terminals;
    std::vector<Symbol> nonterminals;
    std::vector<Rule> rules;
    std::vector<std::pair<std::string, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>>> parseTrees;
    std::vector<std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>>> parseTreesVec;
    Symbol start;
    int nTerminals{};
    int nNonTerminals{};
    std::pair<int, int> contextAmount;
    std::pair<unsigned long, unsigned long> contextSize;
    std::vector<Symbol> actualProduction;
    std::vector<std::vector<Symbol>> words;
    int type{};
    std::random_device rd;

    Grammar(const std::vector<Symbol> &terminals, int nNonTerminals, std::vector<std::vector<Symbol>> words, int type, std::pair<int,int> contextSize);

    void printGrammar();
    std::string grammarToStr();
    void printRules();
    std::string rulesToString();
    void train(int algorithm, int iterations);
    std::pair<double,double> perplexity(const std::vector<std::vector<Symbol>>& testData);
    std::pair<double,double> perplexityKL(const std::vector<std::vector<Symbol>>& testData);

private:
    void baumWelch(int iteration);
    void insideOutside(int iterations);
    void genFPTA();
    void ALERGIA(double alpha);
    void collapsedGibbsSamplePFA(int iterations);
    void trainNGram();
    void generateNonTermnals();
    void generatePermutation(std::vector<std::vector<Symbol>> & permutations, std::vector<Symbol> symbols, unsigned long size, std::vector<Symbol> word, bool context);
    void generateRulesCNF();
    void generateRulesRegular();
    void generateNGramRules();
    void generateNGramNonTerminals();
    void metropolisHastingsPCFG(int iterations);
    void gibbsSamplingPCFG(int iterations);
    void gibbsSamplingPCSG(int iterations);
    void metropolisHastingsPCSG(int iterations);
    double **** CYKProbKLVec(std::vector<Symbol> w);
    /*double **** CYKProbKLVecOpt(std::vector<Symbol> w);*/
    double **** CYKProbKLNVec(const std::vector<Symbol>& w);
    void printInsideTable(double ***p, int wSize) const;
    /*void printInsideTableKL(double ****p, int wSize);*/
    static void printOutsideTable(const std::vector<std::vector<std::vector<double>>>& p);
    void sampleParseTree(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, const std::string& w, double ***insideTable, unsigned int i, unsigned int k);
    void sampleParseTreeKL(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, const std::string& w, double ****insideTable, unsigned int i, unsigned int k);
    void sampleParseTreeKLVec(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, std::vector<Symbol> w, double ****insideTable, unsigned int i, unsigned int k);
    void sampleParseTreeKLVecOpt(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, std::vector<Symbol> w, double ****insideTable, unsigned int i, unsigned int k);
    void calculateNewThetaVecOpt(int i);
    int calculateProductonCounts (const std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>& production, const std::string& w);
    std::vector<std::pair<double, int>> calculateRuleFrequence(Rule r, const std::string& w);
    static double cConstant(std::vector<std::pair<double, int>> ruleFrequence);
    static double probTree (std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> tree);
    static bool equalProductions(std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prodA,
                     std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prodB);
    void updateParseTressTheta();
    void applyProduction(std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prodA);
    Rule findRuleByLHS(std::vector<Symbol> lhs);
    [[nodiscard]] int convertContextToID(int side, std::vector<Symbol> context) const;
    void getActualContext(std::vector<Symbol> & leftContext, std::vector<Symbol> & rightContext);
    /*static void printProduction(
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> tree);
    static void
    printOneProduction(std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prd);*/
    void sampleParseTreeVec(
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> & vr,
            Rule r, std::vector<Symbol> w, double ***pDouble, unsigned int i, unsigned int k);
    double ***CYKProbVec(std::vector<Symbol> w);
    std::vector<std::vector<std::vector<double>>> outsideTable(const std::vector<Symbol>& w, double ***insideTable);
    double ***CYKProbNVec(const std::vector<Symbol>& w);
    void calculateNewThetaVec(const std::vector<Symbol>& w);
    int calculateProductonCountsVec(
            const std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>& pair,
            const std::vector<Symbol>& w);
    static bool equalWord(std::vector<Symbol> w1, std::vector<Symbol> w2);
    double pTiTiMinus1Vec(const std::vector<Symbol>& w);
    double pTiTiMinus1VecOpt(int i) ;
    std::vector<std::pair<double, int>> calculateRuleFrequenceVec(Rule &rule, const std::vector<Symbol>& w);
    void pTiMinus1Frequence(int i);
    void pTiMinus1PlusFrequence(int i);
    void freeInsideTable(double ***p, int wSize) const;
    void freeInsideTableKL(double ****p, int wSize) const;
    void baumWelchExpectation(std::vector<double> &countNd, std::vector<double> &countNl, std::vector<std::vector<std::vector<double>>> &countNT);
    void baumWelchMaximization(std::vector<double> &countNdR, std::vector<double> &countNlR, std::vector<std::vector<std::vector<double>>> &countNTR);
    void calculateWV(std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &W, std::vector<std::vector<std::vector<double>>> &V, const std::vector<Symbol>& w, double***insideTable, std::vector<std::vector<std::vector<double>>>oT);
    bool compatibleALERGIA(const Symbol& a, const Symbol& b, double alpha);
    static bool testALERGIA(double probFa, double freqa, double probFb, double freqb, double alpha);
    void stochasticMerge(const Symbol& a, const Symbol& b);
    void stochasticFold(const Symbol& a, const Symbol& b);
    void removeUnusedNT();
    void recursiveInsertUnused(std::vector<Symbol> &unused, const Symbol& nt);
    void normalizeProbs();
    void sampleRegularRules(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, const std::vector<Symbol>& w);
    void recursiveAddTerminalToNT(Symbol nt, int n, int &nIds);
    void addNGramRuleFrequency(const Symbol& lhs, const Symbol& nextSymbol);

public:
    virtual ~Grammar();
};


#endif //GRAMMARINDCTION_GRAMMAR_H
