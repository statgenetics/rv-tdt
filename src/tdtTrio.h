#ifndef TDTTRIO_H
#define TDTTRIO_H
#include "gw_utilities.h"

/** Others **/
inline void trimNonSenseSites(vector2F& genotype, const vectorL& criteria)
{
    for(UINT i=0; i<genotype.size(); i++){
        for(UINT j=0; j<genotype[i].size(); j++){
            if(!criteria[j]) genotype[i][j]=-9.0;
        }
    }
    return;
}

struct shuffleOutput
{
    shuffleOutput() {}
    shuffleOutput(const vector2F& table, const vector2F& untrans): tdtTable(table), untransmitted(untrans) {}
    vector2F tdtTable, untransmitted;
};

vector2F phasedToUnphased(const vector2F& phased);

/** Base Trio **/
class  trio
{
    public:
        trio() {}
        virtual void load() {}
        virtual ~trio() {}

        bool checkInformTrio() const {return __isInformative;}
        vector2F getTdtTable() const {return __tdtTable;}
        vector2F getNuTransPatChrs() const {return __unTransParentalChrs;}
        vectorL getAnalyzed() const {return __tobeAnalyzed;}
        vectorL getDenovo() const {return __denovoSite;}
        vector2F getGenotype() const {return __genotype;}

        virtual void tdtTableCount() {}
        virtual shuffleOutput shuffle(std::string GenoHapo, gsl_rng* gslr);
        virtual shuffleOutput shuffleWithGenos(std::string GenoHapo, gsl_rng* gslr, vector2F& genos);

    protected:
        vector2F __genotype;
        vectorF __popMafs;
        vectorL __tobeAnalyzed;
        bool __skipMiss;
        vector2F __unTransParentalChrs;
        // four vectors for phased trio, first two vector is for father: v[0] is B count; v[1] is C count; Two vector for unphased trio
        vector2F __tdtTable;
        vectorL __denovoSite;
        bool __PesudoGenoReady;
        bool __isInformative; // check if this trio informative
        vector2F __pesudoGeno; // genotype for shuffle

        vector2F m_GenoShuffle(gsl_rng* gslr);
        vector2F m_HapoShuffle(gsl_rng* gslr);
        virtual void m_setUpPesudoGeno(gsl_rng* gslr) {} // For shuffle: convert unphased to phased, guess missing hapotype
};

/** Phased Trio **/
class phasedTrio: public trio
{
    public:
        phasedTrio() {}
        //~phasedTrio() {}
        void load(const vector2F& genotype, const vectorF& popMafs, const vectorL& tobeAnalyzed, bool skipMiss);
        void tdtTableCount();
        shuffleOutput shuffle(std::string GenoHapo, gsl_rng* gslr);

    protected:
        vectorUI __chrOrigin;
        bool m_ChromOrigin(bool allowMissing, bool allowDenovo);
        void m_phasedTdtCount(vectorF& transParent, vectorF& unTransParent, vectorF& Child, UINT base);
        void m_setUpPesudoGeno(gsl_rng* gslr);
};

bool ChromEqual(const vectorF& ChromOne, const vectorF& ChromTwo, bool Missing, bool Denovo);

/** unPhased Trio **/
class unPhasedTrio: public trio
{
    public:
        unPhasedTrio() {}
        void load(const vector2F& genotype, const vectorF& popMafs, const vectorL& tobeAnalyzed, bool skipMiss);
        void tdtTableCount();
        shuffleOutput shuffle(std::string GenoHapo, gsl_rng* gslr);
        shuffleOutput shuffleWithGenos(std::string GenoHapo, gsl_rng* gslr, vector2F& genos);

    protected:
        void m_unPhasedTdtCount(double fat, double mot, double kid, UINT index);
        void m_setUpPesudoGeno(gsl_rng* gslr);
        double guessMissing(double knownParent, double kid, gsl_rng* gslr);
};

#endif
