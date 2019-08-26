from math import log10

import numpy as np

class HanselSymbol(int):
    """Shim integer that allows an int to masquerade as a str.
    HanselSymbol essentially acts as an int, but when coerced to a str, it will
    use an internal dictionary to map the value of itself to the corresponding
    Hansel symbol."""

    def __new__(cls, symbol_lookup, value):
        cls.symbol_lookup = symbol_lookup
        return int.__new__(cls, value)

    def __str__(self):
        return self.symbol_lookup.get(self, None)


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

    symbols_d : dict{str: int}
        A mapping of symbols to their integer index in the matrix

    symbols_i : dict{int: str}
        A mapping of integer indices to their corresponding symbol

    valid_symbols_i : dict{int: str}
        A mapping of integer indices to their corresponding symbol, iff the symbol
        does not appear in the given list of unsymbols. `valid_symbols_i` is therefore
        a subset of `symbols_i`.

    unsymbols : list{str}
        A list of states that represent an absence of a symbol.

    is_weighted : boolean
        Whether or not the underlying numpy array has been modified by the
        `reweight_observation` function at least once.

    L : int, optional(default=1)
        The number of positions back from the current position (inclusive)
        to consider when calculating conditional probabilities.

    """
    def __new__(cls, input_arr, symbols, unsymbols, L=0):
        # Force our class on the input_arr
        #TODO Is there an overhead in casting a view here (are we copying the
        #     big matrix to a new object? :(
        obj = np.asarray(input_arr).view(cls)

        obj.n_slices= 0
        obj.n_crumbs = 0

        obj.symbols = symbols
        obj.symbols_d = {symbol: i for i, symbol in enumerate(symbols)}
        obj.symbols_i = {i: symbol for i, symbol in enumerate(symbols)}

        obj.is_weighted = False

        obj.unsymbols = unsymbols
        obj.L = L

        obj.valid_symbols_i = {i: symbol for i, symbol in enumerate(symbols) if symbol not in unsymbols}
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
        self.L = getattr(obj, 'L', 0)

        self.symbols_i = getattr(obj, 'symbols_i', {})
        self.valid_symbols_i = getattr(obj, 'valid_symbols_i', {})

        #TODO Safer warning?
        if self.symbols is None or len(self.symbols) == 0:
            import sys
            sys.stderr.write("[FAIL] Attempted to allocate Hansel structure without symbols.\n")
            sys.exit(1)

    @staticmethod
    def init_matrix(symbols, unsymbols, n_positions):
        from multiprocessing import Array
        import ctypes

        n_symbols = len(symbols)
        n_positions += 2 # Add a position for the start source and end sink
        hanselx = np.frombuffer(Array(ctypes.c_float, n_symbols * n_symbols * n_positions * n_positions, lock=False), dtype=ctypes.c_float)
        hanselx = hanselx.reshape(n_symbols, n_symbols, n_positions, n_positions)
        hanselx.fill(0.0) # Initialise as empty
        return Hansel(hanselx, symbols, unsymbols)

    def __orient_positions(self, i, j):
        if i < j:
            return i, j
        else:
            return j, i

    #def __orient_symbols(self, symbol_a, symbol_b, pos_from, pos_to, mirror=False):
    #    if not mirror:
    #        if pos_from < pos_to:
    #            return symbol_a, symbol_b, pos_from, pos_to
    #        else:
    #            return symbol_b, symbol_a, pos_to, pos_from
    #    else:
    #        if pos_from > pos_to:
    #            return symbol_a, symbol_b, pos_from, pos_to
    #        else:
    #            return symbol_b, symbol_a, pos_to, pos_from

    @property
    def sources(self):
        """An alias for `n_slices`"""
        return self.n_slices

    @property
    def observations(self):
        """An alias for `n_crumbs`"""
        return self.n_crumbs

    def add_observation(self, symbol_from, symbol_to, pos_from, pos_to, value=1):
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

        value : float, optional(default=1)
            Magnitude of the observation (defaults to 1).

        """
        #self.n_crumbs += 1 # This doesn't work when updating the matrix in parallel, users should set it manually instead
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        self[self.__symbol_num(symbol_from), self.__symbol_num(symbol_to), pos_from, pos_to] += value

    def __get_observation(self, symbol_from, symbol_to, pos_from, pos_to):
        pos_from, pos_to = self.__orient_positions(pos_from, pos_to)
        v = self[symbol_from, symbol_to, pos_from, pos_to]

        # reweight_matrix destroys evidence without checking for <1 cells for speed
        # so let's check for it here instead
        if v >= 1.0:
            return v
        else:
            return 0

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
        return self.__get_observation(self.__symbol_num(symbol_from), self.__symbol_num(symbol_to), pos_from, pos_to)

    def reweight_observation(self, symbol_from, symbol_to, pos_from, pos_to, ratio):
        """Alter the number of co-occurrences between a pair of positioned symbols by some ratio.

        .. note:: This function will set :attr:`hansel.hansel.Hansel.is_weighted` to `True`.

        Parameters
        ----------

        symbol_from : :py:class:`hansel.hansel.HanselSymbol`
            The first observed symbol of the pair to be reweighted (in space or time).

        symbol_to : :py:class:`hansel.hansel.HanselSymbol`
            The second observed symbol of the pair to be reweighted (in space or time).

        pos_from : int
            The "position" at which `symbol_from` was observed.

        pos_to : int
            The "position" at which `symbol_to` was observed.

        ratio : float
            The ratio by which to subtract the current number of observations.
            That is, `new_value = old_value - (ratio * old_value)`.
        """

        s_a = symbol_from
        s_b = symbol_to
        pos_i, pos_j = self.__orient_positions(pos_from, pos_to)

        old_v = self[s_a, s_b, pos_i, pos_j]
        new_v = old_v - (ratio*old_v)

        if old_v != 0:
            if new_v < 1:
                # Once the last whole crumb of evidence has been used, set it to 0
                # otherwise we will ~infinitely take smaller and smaller decimal crumbs away
                self[s_a, s_b, pos_i, pos_j] = 0
                return old_v
            else:
                self[s_a, s_b, pos_i, pos_j] = new_v
                return old_v - new_v
        self.is_weighted = True
        return 0.0

    def reweight_matrix(self, ratio):
        total = self.sum()
        self = self - (self*ratio)
        self.is_weighted = True
        return total - self.sum()

    def __symbol_num(self, symbol):
        #TODO Catch potential KeyError
        #TODO Generic mechanism for casing (considering non-alphabetically named states, too...)
        return self.symbols_d[symbol]

    def __symbol_unnum(self, num):
        #TODO Catch potential IndexError
        return self.symbols[num]

    def write_path_support_matrix(self, fname, path_a, path_b=None):
        # from    to  count_on    count_off   marginal_on marginal_off
        fh = open(fname, 'w')
        fh.write("\t".join([
            "symbol_a",
            "symbol_b",
            "i",
            "j",
            "obs",
            "total",
        ])+'\n')

        def count_and_write(s_a, s_b, pos_i, pos_j, flip=False):
            #s_a, s_b, pos_i, pos_j = self.__orient_symbols(symbol_a, symbol_b, i, j, mirror=mirror)
            obs = self.__get_observation(self.__symbol_num(s_a), self.__symbol_num(s_b), pos_i, pos_j)
            total = self.get_counts_at(pos_j)["total"]

            i, j = pos_i, pos_j
            if flip:
                i, j = pos_j, pos_i

            fh.write("\t".join([str(x) for x in [
                symbol_a,
                symbol_b,
                i,
                j,
                obs,
                total,
            ]]) + '\n')

        for i, symbol_a in enumerate(path_a):
            for j, symbol_b in enumerate(path_a):
                if path_b:
                    if i > j:
                        # easier to do this here whatev
                        continue
                count_and_write(symbol_a, symbol_b, i, j)
        if path_b:
            for i, symbol_a in enumerate(path_b):
                for j, symbol_b in enumerate(path_b):
                    if i < j:
                        continue
                    count_and_write(symbol_b, symbol_a, j, i, flip=True)


        fh.close()

    def load_hansel_dump(self, prefix):
        import glob
        dump_names = glob.glob(prefix+"*hansel*txt")
        for dump_fn in dump_names:
            symbol_a, symbol_b = dump_fn.split('.')[-2].split("~")
            self[self.symbols_d[symbol_a], self.symbols_d[symbol_b]] = np.loadtxt(dump_fn, delimiter=',')
        return(dump_names)

    def save_hansel_dump(self, prefix):
        for symbol_a in self.symbols_d:
            for symbol_b in self.symbols_d:
                dump_fn = prefix + ".hansel.%s~%s.txt" % (symbol_a, symbol_b)
                np.savetxt(dump_fn, self[self.symbols_d[symbol_a], self.symbols_d[symbol_b]], delimiter=',')

    def get_counts_at(self, at_pos):
        """Get the counts for each symbol that appears at a given position.

        Parameters
        ----------

        at_pos : int
            The "position" for which to return the number of occurrences of each symbol.

        Returns
        -------
        Symbol counts : dict{str: float}
            A dictionary whose keys are each of the symbols that were observed
            at position `at_pos` and a special "total" key. The values are the
            number of observations of that symbol at `at_pos`. The "total" is
            the sum of all observation counts.
        """
        marg = {"total": 0}
        if(at_pos == 0):
            permitted_a_symbols = self.symbols_i
        else:
            permitted_a_symbols = self.valid_symbols_i

        for symbol_a in permitted_a_symbols:
            obs = 0
            #TODO Should this be the other way around?
            for symbol_b in self.symbols_i:
                obs += self.__get_observation(symbol_a, symbol_b, at_pos, at_pos+1)

            if obs > 0:
                marg[HanselSymbol(self.symbols_i, symbol_a)] = obs
                marg["total"] += obs

        return marg

    def get_marginal_of_at(self, of_symbol, at_symbol):
        """Get the marginal distribution of a symbol appearing at a position.

        Parameters
        ----------

        of_symbol : :py:class:`hansel.hansel.HanselSymbol`
            The symbol for which to calculate the marginal distribution.

        at_symbol : :py:class:`hansel.hansel.HanselSymbol`
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

    def get_edge_weights_at(self, symbol_pos, current_path, debug=False):
        """Get the outgoing weighted edges at some position, given a path to that position.

        Parameters
        ----------

        symbol_pos : int
            The index of the current position.

        current_path : list{HanselSymbol}
            A list of symbols representing the path of selected symbols that led
            to the current position, `symbol_pos`.

        Returns
        -------
        Conditional distribution : dict{str: float}
            A dictionary whose keys are each of the possible symbols that may
            be reached from the current position, given the observed path.
            The values are `log10` conditional probabilities of the next symbol
            in the path (or sequence) being that of the key.
        """
        # ...work out each probability, for each branch
        curr_branches = {}
        for symbol in self.get_counts_at(symbol_pos):
            if str(symbol) in self.unsymbols or str(symbol) == "total":
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
                    if debug:
                        print("%.15f" % log10(self.get_conditional_of_at(current_path[curr_i], symbol, curr_i, symbol_pos)), current_path[curr_i], symbol, curr_i, symbol_pos)

            # Append the marginal of symbol at desired position
            curr_branches[symbol] += log10(self.get_marginal_of_at(symbol, symbol_pos))

        return curr_branches

    #TODO Given/predicted is a bit misleading as they turn out to be the "wrong way around"
    def get_conditional_of_at(self, symbol_from, symbol_to, pos_from, pos_to):
        """Given a symbol and position, calculate the conditional for co-occurrence with another positioned symbol.

        Parameters
        ----------

        symbol_from : :py:class:`hansel.hansel.HanselSymbol`
            The first (given) symbol of the pair on which to condition.

        symbol_to : :py:class:`hansel.hansel.HanselSymbol`
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
        obs = self.__get_observation(symbol_from, symbol_to, pos_from, pos_to)
        total = self.get_spanning_support(symbol_to, pos_from, pos_to)
        total_from = self.get_counts_at(pos_from)["total"]
        marg_sym_from = self.get_marginal_of_at(symbol_from, pos_from)

        valid_symbols_seen = 0
        for s in self.valid_symbols_i:
            if s in marg_from:
                valid_symbols_seen += 1
        return self.__estimate_conditional(valid_symbols_seen, obs, total)
        #return self.__estimate_conditional_wmarginal(valid_symbols_seen, obs, total_from, total, marg_sym_from)

    #TODO Should this be "number of sources", rather than "number of observations"
    def get_spanning_support(self, symbol_to, pos_from, pos_to):
        """Get the number of observations that span over two positions of interest, that also feature some symbol.

        Parameters
        ----------

        symbol_to : :py:class:`hansel.hansel.HanselSymbol`
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
        for symbol_from in self.valid_symbols_i:
            total += self.__get_observation(symbol_from, symbol_to, pos_from, pos_to)
        return total

    def __estimate_conditional(self, av, obs, total):
        return (1 + obs)/float(av + total)
        #if obs < 1:
        #    return 0.0000001
        #else:
        #    #return (1 + obs)/float(av + total)
        #    return (obs)/float(total)

    def __estimate_conditional_wmarginal(self, av, obs, total_from, span_j, marg_from_symbol):
        # this is a cool idea but not sure if statistically sound because
        # a) does it destroy the calculation that this is a conditional, not a marginal
        # b) it considers the delta between total_from and span_j, which becomes ~invalid after the first reweight...
        delta = total_from - span_j
        if delta < 0:
            delta = 0
        c = ((delta * marg_from_symbol) + obs + 1) / float(av + span_j + delta)
        return c

