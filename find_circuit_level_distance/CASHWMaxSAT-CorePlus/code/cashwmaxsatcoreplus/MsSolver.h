/**************************************************************************************[MsSolver.h ]
  Copyright (c) 2018-2019, Marek Piotrów

  Based on PbSolver.h ( Copyright (c) 2005-2010, Niklas Een, Niklas Sorensson)

  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
  associated documentation files (the "Software"), to deal in the Software without restriction,
  including without limitation the rights to use, copy, modify, merge, publish, distribute,
  sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all copies or
  substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
  NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
  OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
  **************************************************************************************************/

#ifndef MsSolver_h
#define MsSolver_h

#include "PbSolver.h"
#include <vector>
#include <algorithm>
#include <string.h>
#include "Sort.h"

Int evalGoal(const vec<Pair<weight_t, Minisat::vec<Lit>* > >& soft_cls, vec<bool>& model, Minisat::vec<Lit>& soft_unsat);

static inline int hleft (int i)  { return i * 2; }
static inline int hright(int i)  { return i * 2 + 1; }
static inline int hparent(int i) { return i / 2; }

#define __STR_FUNCTION__ _CutParenthesesNTail(std::string(__PRETTY_FUNCTION__))
#define __Print__ true

std::string _CutParenthesesNTail(std::string&& prettyFuncon); 

class IntLitQueue {
  private:
    vec<Pair<Int, Lit> > heap;
    std::map<Lit, int> core_size;

    bool cmp(int x, int y) { 
        if(weighted_instance)
            return heap[x].fst > heap[y].fst;
        else
        {
            if(heap[x].fst != heap[y].fst)  return heap[x].fst > heap[y].fst;
            else                            return core_size[heap[x].snd] < core_size[heap[y].snd];
        }
    }

  public:
    bool weighted_instance = false;

    IntLitQueue() { heap.push(Pair_new(1, lit_Undef)); }

    bool empty() { return heap.size() <= 1; }

    int size() { return heap.size(); }

    const vec<Pair<Int, Lit> >& getHeap() const { return heap; }

    void clear() { heap.shrink(heap.size() - 1); }

    const Pair<Int, Lit>& top() { return heap[1]; }

    void push(Pair<Int, Lit> p, int core_size = 0) { 
        this->core_size[p.snd] = core_size;
        heap.push();
        int i = heap.size() - 1;
        heap[0] = std::move(p);
        while (hparent(i) != 0 && cmp(0, hparent(i))) { // percolate up
            heap[i] = std::move(heap[hparent(i)]);
            i       = hparent(i);
        }
        heap[i] = std::move(heap[0]);
    }

    void pop(void) {
        heap[1] = std::move(heap.last());
        heap.pop();
        if (heap.size() > 1) { // percolate down
            int i = 1;
            heap[0] = std::move(heap[1]);
            while (hleft(i) < heap.size()){
                int child = hright(i) < heap.size() && cmp(hright(i), hleft(i)) ? hright(i) : hleft(i);
                if (!cmp(child, 0)) break;
                heap[i] = std::move(heap[child]);
                i       = child;
            }
            heap[i] = std::move(heap[0]);
        }
    }

} ;

template<class Comp>
class Strat_Array {

    Comp     lt;
    vec<int> array;       
    vec<int> indices;
    
    int strat_indices;
    
    public:
        Strat_Array(const Comp& c) : lt(c), strat_indices(0) { }

        int  size      ()          const { return array.size(); }
        bool empty     ()          const { return array.size() == 0; }
        int  get_strat ()          const { return strat_indices; }
        int  operator[](int index) const { assert(index < array.size()); return array[index]; }

        void update(int v) {
            assert(v < indices.size());
            if(indices[v] < 0) exit(0);
            if(lt(v) && indices[v] >= strat_indices) swap_array(array, indices, v, strat_indices++);
            else if(!lt(v) && indices[v] < strat_indices) swap_array(array, indices, v, --strat_indices);
            return;
        }

        void insert(int n) {
            indices.growTo(n + 1, -1);
            if(indices[n] > 0) exit(0);
            indices[n] = array.size();
            array.push(n);
            update(n); 
            return;
        }

        void clear(bool dealloc = false) {
            strat_indices = 0;
            for (int i = 0; i < array.size(); i++) indices[array[i]] = -1;
            array.clear(dealloc); 
            return;
        }

        void swap(int x, int pos) {
            swap_array(array, indices, x, pos);
            return;
        }

    private:
        void swap_array(vec<int>& array, vec<int>& indices, int v, int pos) {
            array[indices[v]] = array[pos];
            indices[array[pos]] = indices[v];
            indices[v] = pos;
            array[pos] = v;
            return;
        }

        void remove_array(vec<int>& array, vec<int>& indices, int v) {
            array[indices[v]] = array[array.size() - 1];
            indices[array.size() - 1] = indices[v];
            array.pop();
            indices[v] = -1;
            return;
        }        

};

#ifdef USE_SCIP
//#include <vector>
//#include <algorithm>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#endif

class MsSolver : public PbSolver {
  public:
    MsSolver(bool use_preprocessing = false) : 
          PbSolver(use_preprocessing)
        , harden_goalval(0)
        , fixed_goalval(0)
        , goal_gcd(1) {}

    Int                 harden_goalval,  //  Harden goalval used in the MaxSAT preprocessing 
                        fixed_goalval;   // The sum of weights of soft clauses that must be false
    vec<Pair<weight_t, Minisat::vec<Lit>* > > orig_soft_cls; // Soft clauses before preprocessing by MaxPre; empty if MaxPre is not used
    vec<Pair<weight_t, Minisat::vec<Lit>* > > soft_cls; // Relaxed non-unit soft clauses with weights; a relaxing var is the last one in a vector. 
    std::vector<std::vector<Lit>>    ls_hard_cls;
    std::vector<std::vector<Lit>>    ls_soft_cls;
    std::vector<weight_t>    ls_hard_weight;
    std::vector<weight_t>    ls_soft_weight;
    weight_t            goal_gcd; // gcd of soft_cls weights
    int                 top_for_strat, top_for_hard; // Top indices to soft_cls for stratification and hardening operations.
    Map<Lit, Int>       harden_lits;    // The weights of literals included into "At most 1" clauses (MaxSAT preprocessing of soft clauese).
    vec<Pair<Lit,Int> > am1_rels;       // The weights of relaxing vars in "At most 1" clauses

    void    storeSoftClause(const vec<Lit>& ps, weight_t weight) {
                Minisat::vec<Lit> *ps_copy = new Minisat::vec<Lit>; 
                for (int i = 0; i < ps.size(); i++) ps_copy->push(ps[i]); 
                soft_cls.push(Pair_new(weight, ps_copy)); }


    void    harden_soft_cls(Minisat::vec<Lit>& assump_ps, vec<Int>& assump_Cs, vec<weight_t>& sorted_assump_Cs, IntLitQueue& delayed_assump, Int& delayed_assump_sum);
    void    optimize_last_constraint(vec<Linear*>& constrs, Minisat::vec<Lit>& assump_ps, Minisat::vec<Lit>& new_assump);

#ifdef USE_SCIP
    bool scip_solve(const Minisat::vec<Lit> *assump_ps, const vec<Int> *assump_Cs, const IntLitQueue *delayed_assump,
            bool weighted_instance, int sat_orig_vars, int sat_orig_cls);
#endif    

    bool    merge(Minisat::vec<Lit>& ps_1, Minisat::vec<Lit>& ps_2);

    void    maxsat_solve(solve_Command cmd = sc_Minimize); 
    void    preprocess_soft_cls(Minisat::vec<Lit>& assump_ps, vec<Int>& assump_Cs, const Lit max_assump, const Int& max_assump_Cs, 
                                           IntLitQueue& delayed_assump, Int& delayed_assump_sum);
} ;

#endif
