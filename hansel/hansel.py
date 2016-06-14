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

    def add_observation(self, symbol_from, symbol_to, pos_from, pos_to):
        self.n_crumbs += 1
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        self[self.__symbol_num(symbol_from)][self.__symbol_num(symbol_to)][pos_from][pos_to] += 1

    def get_observation(self, symbol_from, symbol_to, pos_from, pos_to):
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        return self[self.__symbol_num(symbol_from)][self.__symbol_num(symbol_to)][pos_from][pos_to]

    def reweight_observation(self, symbol_from, symbol_to, pos_from, pos_to, ratio):
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        old_v = self[self.__symbol_num(symbol_from)][self.__symbol_num(symbol_to)][pos_from][pos_to]
        new_v = old_v - (ratio*old_v)

        if old_v != 0:
            if new_v < 1:
                self[self.__symbol_num(symbol_from)][self.__symbol_num(symbol_to)][pos_from][pos_to] = 0
            else:
                print "Reducing support between %d(%s) -> %d(%s) by %.2f (%.2f -> %.2f)" % (pos_from, symbol_from, pos_to, symbol_to, ratio, old_v, new_v)
                self[self.__symbol_num(symbol_from)][self.__symbol_num(symbol_to)][pos_from][pos_to] = new_v
        self.is_weighted = True

    def __symbol_num(self, symbol):
        #TODO Catch potential IndexError
        #TODO Generic mechanism for casing (considering non-alphabetically named states, too...)
        return self.symbols.index(symbol)

    def __symbol_unnum(self, num):
        #TODO Catch potential IndexError
        return self.symbols[num]

    def get_marginal_at(self, at_symbol):
        marg = {"total": 0}
        for symbol_a in self.symbols:
            obs = 0
            for symbol_b in self.symbols:
                obs += self.get_observation(symbol_a, symbol_b, at_symbol, at_symbol+1)

            if obs > 0:
                marg[symbol_a] = obs
                marg["total"] += obs

        #print "\t[MARG] %d %s" % (symbol_pos, str(marg))
        return marg

    def get_marginal_of_at(self, of_symbol, at_symbol):
        marginal = self.get_marginal_at(at_symbol)
        return marginal[of_symbol] / marginal["total"]

    def get_edge_weights_at(self, symbol_pos, current_path, L=5):
        # ...work out each probability, for each branch
        curr_branches = {}
        curr_branches_tot = 0.0
        for symbol in self.get_marginal_at(symbol_pos):
            if symbol == "total":
                continue

            curr_branches[symbol] = 0.0

            # ...with conditionals on each part of LAST L ENTRIES in the current path
            # (where the case for there being less than L entries in the path being
            #  accounted for)
            LI = L
            if len(current_path) <= L:
                LI = len(current_path)-1

            for l in range(0, LI):
                curr_i = len(current_path) - LI + l
                curr_branches[symbol] += log10(self.get_conditional_of_at(current_path[curr_i], symbol, curr_i, symbol_pos))

            curr_branches[symbol] += log10(self.get_marginal_of_at(symbol, symbol_pos))
            curr_branches_tot += curr_branches[symbol]

        return curr_branches, curr_branches_tot

    #TODO This doesn't belong here, the API should provide an interface to the
    #     underlying evidence only, we need not concern ourselves with the end use.
    def select_next_edge_at(self, symbol_pos, current_path, L=5):
        curr_branches = self.get_edge_weights_at(symbol_pos, current_path, L=L)[0]
        print "\t[TREE] %s" % curr_branches
        # Return the symbol and probability of the next base to add to the
        # current path based on the best marginal
        next_v = 0.0
        next_m = None

        """
        if ENTROPY.any() and random.random() > ENTROPY[i] and path_i > 0 and symbol == path_i:
            for symbol in curr_branches:
                if max_v is None:
                    max_v = curr_branches[symbol]
                    max_m = symbol
                elif curr_branches[symbol] < max_v:
                    max_v = curr_branches[symbol]
                    max_m = symbol
        """
        for symbol in curr_branches:
            if symbol == "total":
                continue
            if next_m is None:
                next_v = curr_branches[symbol]
                next_m = symbol
            elif curr_branches[symbol] > next_v:
                next_v = curr_branches[symbol]
                next_m = symbol

        return next_m, next_v

    def get_conditional_of_at(self, symbol_from, symbol_to, pos_from, pos_to):
        marg_from = self.get_marginal_at(pos_from) #TODO pos_from - 1?
        obs = self.get_observation(symbol_from, symbol_to, pos_from, pos_to)
        total = self.get_spanning_support(symbol_to, pos_from, pos_to)
        return self.__estimate_conditional(len(marg_from)-1, obs, total)

    def get_spanning_support(self, symbol_to, pos_from, pos_to):
        total = 0
        for symbol in self.symbols:
            total += self.get_observation(symbol, symbol_to, pos_from, pos_to)
        return total

    def __estimate_conditional(self, av, obs, total):
        return (1 + obs)/float(av + total)

    #TODO
    def validate_path(self):
        pass

    #TODO
    def probability_path(self):
        pass
