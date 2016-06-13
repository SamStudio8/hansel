from math import log,log10

import numpy as np
#TODO Abstract symbols

#NOTE
# Although it is possible to 'fully' subclass numpy's ndarray, the
# accepted practice appears to be instructing users to create an ndarray first
# and then pass it to an object constructor responsible for 'converting' that
# ndarray to one of your own type.
# The benefit appears to be a simpler constructor (it need not emulate the
# actual constructor of ndarray) and I guess we are somewhat less coupled
# to the ndarray interface, too.
# You can read more about subclassing ndarray via the official documentation:
# http://docs.scipy.org/doc/numpy-1.10.1/user/basics.subclassing.html
class Hansel(np.ndarray):

    def __new__(cls, input_arr):
        # Force our class on the input_arr
        obj = np.asarray(input_arr).view(cls)
        obj.is_weighted = False
        obj.symbols = ['A', 'C', 'G', 'T', 'N']
        obj.n_slices= 0
        obj.n_crumbs = 0
        return obj

    #NOTE Provides support for construction mechanisms of numpy
    def __array_finalize__(self, obj):
        if obj is None:
            # We're in a __new__ call (so we're covered)
            return

        # View casting or from-template, add the missing info from the
        # existing object before finalizing
        self.is_weighted = getattr(obj, 'is_weighted', False)
        self.symbols = getattr(obj, 'symbols', ['A', 'C', 'G', 'T', 'N'])
        self.n_slices = getattr(obj, 'n_slices', 0)
        self.n_crumbs = getattr(obj, 'n_crumbs', 0)

    def __orient_positions(self, i, j):
        if i < j:
            return i, j
        else:
            return j, i

    @property
    def sources(self):
        return self.n_slices

    @property
    def observations(self):
        return self.n_crumbs

    def add_observation(self, snp_from, snp_to, pos_from, pos_to):
        self.n_crumbs += 1
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        self[self.__base_num(snp_from)][self.__base_num(snp_to)][pos_from][pos_to] += 1

    def get_observation(self, snp_from, snp_to, pos_from, pos_to):
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        return self[self.__base_num(snp_from)][self.__base_num(snp_to)][pos_from][pos_to]

    def reweight_observation(self, snp_from, snp_to, pos_from, pos_to, ratio):
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        old_v = self[self.__base_num(snp_from)][self.__base_num(snp_to)][pos_from][pos_to]
        new_v = old_v - (ratio*old_v)

        if old_v != 0:
            if new_v < 1:
                self[self.__base_num(snp_from)][self.__base_num(snp_to)][pos_from][pos_to] = 0
            else:
                print "Reducing support between %d(%s) -> %d(%s) by %.2f (%.2f -> %.2f)" % (pos_from, snp_from, pos_to, snp_to, ratio, old_v, new_v)
                self[self.__base_num(snp_from)][self.__base_num(snp_to)][pos_from][pos_to] = new_v
        self.is_weighted = True

    def __base_num(self, base):
        #TODO Catch potential IndexError
        return self.symbols.index(base.upper())

    def __base_unnum(self, num):
        #TODO Catch potential IndexError
        return self.symbols[num]

    def get_marginal_at(self, at_snp):
        marg = {"total": 0}
        for snp_a in self.symbols:
            obs = 0
            for snp_b in self.symbols:
                obs += self.get_observation(snp_a, snp_b, at_snp, at_snp+1)

            if obs > 0:
                marg[snp_a] = obs
                marg["total"] += obs

        #print "\t[MARG] %d %s" % (snp_pos, str(marg))
        return marg

    def get_marginal_of_at(self, of_snp, at_snp):
        marginal = self.get_marginal_at(at_snp)
        return marginal[of_snp] / marginal["total"]

    def get_edge_weights_at(self, snp_pos, current_path, L=5):
        # ...work out each probability, for each branch
        curr_branches = {}
        curr_branches_tot = 0.0
        for mallele in self.get_marginal_at(snp_pos):
            if mallele == "total":
                continue

            curr_branches[mallele] = 0.0

            # ...with conditionals on each part of LAST L ENTRIES in the current path
            # (where the case for there being less than L entries in the path being
            #  accounted for)
            LI = L
            if len(current_path) <= L:
                LI = len(current_path)-1

            for l in range(0, LI):
                curr_i = len(current_path) - LI + l
                curr_branches[mallele] += log10(self.get_conditional_of_at(current_path[curr_i], mallele, curr_i, snp_pos))

            curr_branches[mallele] += log10(self.get_marginal_of_at(mallele, snp_pos))
            curr_branches_tot += curr_branches[mallele]

        return curr_branches, curr_branches_tot

    #TODO This doesn't belong here, the API should provide an interface to the
    #     underlying evidence only, we need not concern ourselves with the end use.
    def select_next_edge_at(self, snp_pos, current_path, L=5):
        curr_branches = self.get_edge_weights_at(snp_pos, current_path, L=L)[0]
        print "\t[TREE] %s" % curr_branches
        # Return the mallele and probability of the next base to add to the
        # current path based on the best marginal
        next_v = 0.0
        next_m = None

        """
        if ENTROPY.any() and random.random() > ENTROPY[i] and path_i > 0 and snp == path_i:
            for mallele in curr_branches:
                if max_v is None:
                    max_v = curr_branches[mallele]
                    max_m = mallele
                elif curr_branches[mallele] < max_v:
                    max_v = curr_branches[mallele]
                    max_m = mallele
        """
        for mallele in curr_branches:
            if mallele == "total":
                continue
            if next_m is None:
                next_v = curr_branches[mallele]
                next_m = mallele
            elif curr_branches[mallele] > next_v:
                next_v = curr_branches[mallele]
                next_m = mallele

        return next_m, next_v

    def get_conditional_of_at(self, snp_from, snp_to, pos_from, pos_to):
        marg_from = self.get_marginal_at(pos_from) #TODO pos_from - 1?
        obs = self.get_observation(snp_from, snp_to, pos_from, pos_to)
        total = self.get_spanning_support(snp_to, pos_from, pos_to)
        return self.__estimate_conditional(len(marg_from)-1, obs, total)

    def get_spanning_support(self, snp_to, pos_from, pos_to):
        total = 0
        for symbol in self.symbols:
            total += self.get_observation(symbol, snp_to, pos_from, pos_to)
        return total

    def __estimate_conditional(self, av, obs, total):
        return (1 + obs)/float(av + total)

    #TODO
    def validate_path(self):
        pass

    #TODO
    def probability_path(self):
        pass
