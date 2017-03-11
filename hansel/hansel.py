from math import log10

import numpy as np

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
    """API to a numpy-backed graph-inspired data structure for determining
    likely sequences of symbols from breadcrumbs of evidence.

    Given a numpy array and a list of permitted symbols, Hansel presents a
    user friendly API to store and operate on a graph-like data structure
    that can be used to investigate likely sequences of symbols based on
    observations of pairwise co-occurrence of those symbols in a dimension
    such as space or time.

    Parameters
    ----------

    input_arr : 4D numpy array
        A numpy array (typically initialised with zeros) of size (A, A, B+2, B+2),
        where `A` is the number of `symbols` + `unsymbols` and `B` are the points in time or
        space on which pairwise observations between symbols can be observed.

    symbols : list{str}
        A list of permitted states or symbols as strings.

    unsymbols : list{str}
        A list of permitted states that may represent the known absence of a symbol.

    Attributes
    ----------

    n_slices(sources) : int
        The number of distinct sources from which observations were received.

    n_crumbs(observations) : int
        The number of pairwise observations that were observed.

    symbols : list{str}
        The list of permitted states or symbols.

    unsymbols : list{str}
        A list of states that represent an absence of a symbol.

    is_weighted : boolean
        Whether or not the underlying numpy array has been modified by the
        `reweight_observation` function at least once.

    L : int, optional(default=1)
        The number of positions back from the current position (inclusive)
        to consider when calculating conditional probabilities.

    """
    def __new__(cls, input_arr, symbols, unsymbols, L=1):
        # Force our class on the input_arr
        #TODO Is there an overhead in casting a view here (are we copying the
        #     big matrix to a new object? :(
        obj = np.asarray(input_arr).view(cls)

        obj.is_weighted = False
        obj.symbols = symbols
        obj.symbols_d = {symbol: i for i, symbol in enumerate(symbols)}
        obj.unsymbols = unsymbols
        obj.n_slices= 0
        obj.n_crumbs = 0
        obj.L = L

        obj.valid_symbols = set(symbols) - set(unsymbols)
        return obj

    #NOTE Provides support for construction mechanisms of numpy
    def __array_finalize__(self, obj):
        #NOTE This is probably quite gross.
        #     I am quite sure that `obj` will never be `None` here (despite docs)
        #     as __new__ now invokes a view cast and thus always sends us to
        #     __array_finalize__ *before* it has had a chance to set variables.
        #     I try and catch us inside a "__new__ view cast" by checking whether
        #     obj is a plain ndarray (which it *should* be if this is new).
        #     Elsewise type(obj) should be hansel.Hansel for a view cast or
        #     from-template on an existing Hansel structure.
        if obj is None or type(obj) == np.ndarray:
            return

        #NOTE Defaults are set here, not in __new__, as it is this method
        #     that actually oversees creation of all objects (__new__) sends
        #     us here anyway... Unfortunately this leaves us no choice but to
        #     provide an impractical default for the symbols list...
        # View casting or from-template, add the missing info from the
        # existing object before finalizing
        self.is_weighted = getattr(obj, 'is_weighted', False)
        self.n_slices = getattr(obj, 'n_slices', 0)
        self.n_crumbs = getattr(obj, 'n_crumbs', 0)
        self.symbols = getattr(obj, 'symbols', [])
        self.symbols_d = getattr(obj, 'symbols_d', {})
        self.unsymbols = getattr(obj, 'unsymbols', [])
        self.valid_symbols = getattr(obj, 'valid_symbols', set())
        self.L = getattr(obj, 'L', 1)

        #TODO Safer warning?
        if self.symbols is None or len(self.symbols) == 0:
            import sys
            sys.stderr.write("[FAIL] Attempted to allocate Hansel structure without symbols.\n")
            sys.exit(1)

    def __orient_positions(self, i, j):
        if i < j:
            return i, j
        else:
            return j, i

    @property
    def sources(self):
        """An alias for `n_slices`"""
        return self.n_slices

    @property
    def observations(self):
        """An alias for `n_crumbs`"""
        return self.n_crumbs

    def add_observation(self, symbol_from, symbol_to, pos_from, pos_to):
        """Add a pairwise observation to the data structure.

        Parameters
        ----------

        symbol_from : str
            The first observed symbol of the pair (in space or time).

        symbol_to : str
            The second observed symbol of the pair (in space or time).

        pos_from : int
            The "position" at which `symbol_from` was observed.

        pos_to : int
            The "position" at which `symbol_to` was observed.

        """
        self.n_crumbs += 1
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        self[self.__symbol_num(symbol_from), self.__symbol_num(symbol_to), pos_from, pos_to] += 1

    def get_observation(self, symbol_from, symbol_to, pos_from, pos_to):
        """Get the number of co-occurrences of a pair of positioned symbols.

        Parameters
        ----------

        symbol_from : str
            The first observed symbol of the pair (in space or time).

        symbol_to : str
            The second observed symbol of the pair (in space or time).

        pos_from : int
            The "position" at which `symbol_from` was observed.

        pos_to : int
            The "position" at which `symbol_to` was observed.

        Returns
        -------
        Number of observations : float
            The number of times `symbol_from` was observed at `pos_from` with
            `symbol_to` at `pos_to` across all sources (slices).

            .. note::
                It is possible for the number of observations returned to be
                a `float` if :attr:`hansel.hansel.Hansel.is_weighted`
                is `True`.
        """
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        return self[self.__symbol_num(symbol_from), self.__symbol_num(symbol_to), pos_from, pos_to]

    #TODO Describe the conditions under which reweighting halts in the documentation
    def reweight_observation(self, symbol_from, symbol_to, pos_from, pos_to, ratio):
        """Alter the number of co-occurrences between a pair of positioned symbols by some ratio.

        .. note:: This function will set :attr:`hansel.hansel.Hansel.is_weighted` to `True`.

        Parameters
        ----------

        symbol_from : str
            The first observed symbol of the pair to be reweighted (in space or time).

        symbol_to : str
            The second observed symbol of the pair to be reweighted (in space or time).

        pos_from : int
            The "position" at which `symbol_from` was observed.

        pos_to : int
            The "position" at which `symbol_to` was observed.

        ratio : float
            The ratio by which to subtract the current number of observations.
            That is, `new_value = old_value - (ratio * old_value)`.
        """
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        old_v = self[self.__symbol_num(symbol_from), self.__symbol_num(symbol_to), pos_from, pos_to]
        new_v = old_v - (ratio*old_v)

        if old_v != 0:
            if new_v < 1:
                self[self.__symbol_num(symbol_from), self.__symbol_num(symbol_to), pos_from, pos_to] = 0
                return old_v
            else:
                self[self.__symbol_num(symbol_from), self.__symbol_num(symbol_to), pos_from, pos_to] = new_v
                return old_v - new_v
        return 0.0

        #TODO This is a bit gross as we should maybe handle it with Gretel instead
        """
        if self.get_observation(symbol_from, "_", pos_from, pos_from+1) > 0:
            #print "Reducing support between %d(%s) -> %d(%s) by %.2f (%.2f -> %.2f)" % (pos_from, symbol_from, pos_to, symbol_to, ratio, old_v, new_v)
            old_v = self[self.__symbol_num(symbol_from)][self.__symbol_num("_")][pos_from][pos_from+1]

            # Provisional testing seems to indicate this works best on small sets...
            #new_v = old_v - ((ratio*old_v) / (len( list(set(self.symbols) - set(self.unsymbols)) ) ))
            new_v = old_v - (ratio*old_v)

            try:
                marg_from = self.get_counts_at(pos_from+1)
            except IndexError:
                pass
            else:

                valid_symbols_seen = 0
                for s in list(set(self.symbols) - set(self.unsymbols)):
                    if s in marg_from:
                        valid_symbols_seen += 1

                new_v = old_v - ((ratio*old_v)/ valid_symbols_seen)

                if old_v != 0:
                    if new_v < 1:
                        self[self.__symbol_num(symbol_from)][self.__symbol_num("_")][pos_from][pos_from+1] = 0
                    else:
                        self[self.__symbol_num(symbol_from)][self.__symbol_num("_")][pos_from][pos_from+1] = new_v
        """
        self.is_weighted = True

    def __symbol_num(self, symbol):
        #TODO Catch potential KeyError
        #TODO Generic mechanism for casing (considering non-alphabetically named states, too...)
        return self.symbols_d[symbol]

    def __symbol_unnum(self, num):
        #TODO Catch potential IndexError
        return self.symbols[num]

    #TODO What about non-whole genes (need a special symbol)
    def get_counts_at(self, at_pos):
        """Get the counts for each symbol that appears at a given position.

        Parameters
        ----------

        at_pos : int
            The "position" for which to return the number of occurrences of each symbol.

        Returns
        -------
        Symbol counts : dict{str, float}
            A dictionary whose keys are each of the symbols that were observed
            at position `at_pos` and a special "total" key. The values are the
            number of observations of that symbol at `at_pos`. The "total" is
            the sum of all observation counts.
        """
        marg = {"total": 0}
        if(at_pos == 0):
            permitted_a_symbols = self.symbols_d.keys()
        else:
            permitted_a_symbols = self.valid_symbols

        for symbol_a in permitted_a_symbols:
            obs = 0
            for symbol_b in self.symbols_d:
                obs += self.get_observation(symbol_a, symbol_b, at_pos, at_pos+1)

            if obs > 0:
                marg[symbol_a] = obs
                marg["total"] += obs

        #print "\t[MARG] %d %s" % (symbol_pos, str(marg))
        return marg

    def get_marginal_of_at(self, of_symbol, at_symbol):
        """Get the marginal distribution of a symbol appearing at a position.

        Parameters
        ----------

        of_symbol : str
            The symbol for which to calculate the marginal distribution.

        at_symbol : int
            The position at which to calculate the marginal distribution.

        Returns
        -------

        Marginal probability : float
            The probability a random symbol drawn from all observations at
            `at_symbol` being equal to `of_symbol`. Alternatively, the
            proportion of all symbols observed at `at_symbol` being equal
            to `of_symbol`.
        """
        marginal = self.get_counts_at(at_symbol)
        return marginal[of_symbol] / marginal["total"]

    def get_edge_weights_at(self, symbol_pos, current_path):
        """Get the outgoing weighted edges at some position, given a path to that position.

        Parameters
        ----------

        symbol_pos : int
            The index of the current position.

        current_path : list{str}
            A list of symbols representing the path of selected symbols that led
            to the current position, `symbol_pos`.

        Returns
        -------
        Conditional distribution : dict{str, float}
            A dictionary whose keys are each of the possible symbols that may
            be reached from the current position, given the observed path.
            The values are `log10` conditional probabilities of the next symbol
            in the path (or sequence) being that of the key.
        """
        # ...work out each probability, for each branch
        curr_branches = {}
        for symbol in self.get_counts_at(symbol_pos):
            if symbol in self.unsymbols or symbol == "total":
                continue

            curr_branches[symbol] = 0.0

            if symbol_pos > 1:

                # If the length of the path, without the sentinel is less than L,
                # we can only inspect the available members of L so far...
                if len(current_path)-1 < self.L:
                    l_limit = len(current_path)-1
                else:
                    l_limit = self.L

                for l in range(0, l_limit):
                    curr_i = (len(current_path)-1) - l
                    curr_branches[symbol] += log10(self.get_conditional_of_at(current_path[curr_i], symbol, curr_i, symbol_pos))

            # Append the marginal of symbol at desired position
            curr_branches[symbol] += log10(self.get_marginal_of_at(symbol, symbol_pos))

        return curr_branches

    #TODO Given/predicted is a bit misleading as they turn out to be the "wrong way around"
    def get_conditional_of_at(self, symbol_from, symbol_to, pos_from, pos_to):
        """Given a symbol and position, calculate the conditional for co-occurrence with another positioned symbol.

        Parameters
        ----------

        symbol_from : str
            The first (given) symbol of the pair on which to condition.

        symbol_to : str
            The second (predicted) symbol of the pair on which to condition.

        pos_from : int
            The "position" at which the given `symbol_from` was observed.

        pos_to : int
            The "position" at which the predicted `symbol_to` was observed.


        Returns
        -------
        Conditional probability : float
            The conditional probability of `symbol_from` occurring at `pos_from`
            given observation of a predicted `symbol_to` at `pos_to`.
        """
        marg_from = self.get_counts_at(pos_from)
        obs = self.get_observation(symbol_from, symbol_to, pos_from, pos_to)
        total = self.get_spanning_support(symbol_to, pos_from, pos_to)

        valid_symbols_seen = 0
        for s in self.valid_symbols:
            if s in marg_from:
                valid_symbols_seen += 1
        return self.__estimate_conditional(valid_symbols_seen, obs, total)

    #TODO Should this be "number of sources", rather than "number of observations"
    def get_spanning_support(self, symbol_to, pos_from, pos_to):
        """Get the number of observations that span over two positions of interest, that also feature some symbol.

        Parameters
        ----------

        symbol_to : str
            The symbol that should appear at `pos_to`.

        pos_from : int
            A position that appears "before" `pos_to` in space or time that
            must be overlapped by a source to be counted. The symbol at
            `pos_from` is not relevant.

        pos_to : int
            The second position a source must overlap (but not necessarily
            terminate at) that must be the symbol `symbol_to`.


        Returns
        -------
        Number of observations : float
            The number of observations yielded from sources that overlap both
            `pos_from` and `pos_to`, that also feature `symbol_to` at `pos_to`.

            .. note::
                It is possible for the number of observations returned to be
                a `float` if :attr:`hansel.hansel.Hansel.is_weighted`
                is `True`.
        """
        total = 0
        for symbol_from in self.valid_symbols:
            total += self.get_observation(symbol_from, symbol_to, pos_from, pos_to)
        return total

    def __estimate_conditional(self, av, obs, total):
        return (1 + obs)/float(av + total)

