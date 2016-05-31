__author__ = 'eifu'

from numpy.random import sample
import numpy as np


class Histone(object):
    def __init__(self, position=0,
                 kp=0.12,
                 kp2=0.12,
                 km=0.117,
                 a_bool=False,
                 ka=0.12,
                 nextnode=None,
                 prenode=None,
                 copy=False,
                 copy_histone=None):
        """
        Initialize the histone. The position must be
        specified, but K_PLUS, K_MINUS, nextNode, preNode,
        and percentage can be specified. If they are not,
        default values variables will be assigned to those
        instance variables.

        position  ... a position of histone, to differentiate
                  from other histone instances.

        A_bool ... a boolean value, activator of transcription factors
                   if A is true, then the K_ACE is the given paramater.
                   else, K_ACE is 0, which means no histones can be turned to be acetilated.
        preNode ... a neighbor of the instances, judged by position - 1
        nextNode ... a neighbor of the instances, judged by position + 1
        """

        self.status = ""

        # copy the info of histone object
        if copy:
            self.position = copy_histone.position

            self.prenode = copy_histone.prenode  # doubly-linked
            if copy_histone.prenode is not None:
                copy_histone.prenode.nextnode = self

            self.nextnode = copy_histone.nextnode  # doubly-linked
            if copy_histone.nextnode is not None:
                copy_histone.nextnode.prenode = self

            self.K_ACE = copy_histone.K_ACE
            self.K_PLUS = copy_histone.K_PLUS
            self.K_PLUS2 = copy_histone.K_PLUS2
            self.K_MINUS = copy_histone.K_MINUS

        # create new histone object based on the arguments
        else:
            self.position = position
            self.prenode = prenode
            self.nextnode = nextnode
            # if A_bool is false, K_ACE is set to be 0.
            if a_bool:
                self.K_ACE = ka
            else:
                self.K_ACE = 0
            self.K_PLUS = kp
            self.K_PLUS2 = kp2
            self.K_MINUS = km

    def set_adjhistone(self, nextnode):
        """
        Set the input to self.nextNode, one of its instance variables.
        """
        self.nextnode = nextnode
        nextnode.prenode = self

    def set_ka(self, a_bool):
        # set K_ace default value if A_bool is True
        if a_bool:
                self.K_ACE = 0.12
        else:
            self.K_ACE = 0

    def dna_methylation(self):
        return self

    def k_plus(self):
        """
        the unmethylated histone will get methylated by K_PLUS probability
        return methylated object if the histone will get methylated
        """
        return self

    def k_minus(self):
        """
        the methylated histone will be get unmethylated by K_MINUS probability
        or
        the acetilated histone will be get unacetylated by K_MINUS probability
        return unmethylated histone object if the histone will get unmethylated
        """

        if sample() < self.K_MINUS:
            return UHistone(copy=True, copy_histone=self)
        else:
            return self

    def k_ace(self):
        """
        the unmethylated hitone will be get acetilated by K_ACE probability
        return acetilated histone if the histone will get acetilated
        """
        return self

    def __lt__(self, other):
        """
        return true when the position of self is smaller than that of other
        """
        return self.position < other.position

    def __str__(self):
        """
        return all info about the instance of this object.
        """
        if self.status == 'm':
            st = "methylated"
        elif self.status == 'a':
            st = "acetylated"
        else:
            st = "unmethylated"
        sentence = "pos: {}  status: {}\t".format(self.position, st)
        return sentence

    def display(self):
        sentence = "address: {}\nposition: {}\nstatus: {}\n".format(self, self.position, self.status)

        try:
            sentence += "Prev Histone:{}\n".format(self.prenode)
        except AttributeError:
            sentence += "Prev Histone: None\n"
        try:
            sentence += "Next Histone: {}".format(self.nextnode)
        except AttributeError:
            sentence += "Next Histone: None"
        finally:
            print(sentence)
            print("###\n")


class MHistone(Histone):
    def __init__(self, **kwarg):
        super().__init__(**kwarg)
        self.status = "m"


class UHistone(Histone):
    def __init__(self, **kwarg):
        super().__init__(**kwarg)
        self.status = "u"

    def k_plus(self):
        if (self.prenode is not None and self.prenode.status == "m"
            and sample() < self.K_PLUS):
            # preNode is methylated and sample() gets smaller than K_PLUS
            return MHistone(copy=True, copy_histone=self)
        if (self.nextnode is not None and self.nextnode.status == "m"
            and sample() < self.K_PLUS):
            # nextNode is methylated and sample() gets smaller than K_PLUS
            return MHistone(copy=True, copy_histone=self)
        return self

    def k_minus(self):
        return self

    def k_ace(self):
        if sample() < self.K_ACE:
            return AHistone(copy=True, copy_histone=self)
        else:
            return self


class AHistone(Histone):
    def __init__(self, **kwarg):
        super().__init__(**kwarg)
        self.status = "a"

    def k_plus(self):
        if (self.prenode is not None and self.prenode.status == "m"
            and sample() < self.K_PLUS2):
            # preNode is methylated and sample() gets smaller than KPLUS2
            return UHistone(copy=True, copy_histone=self)
        if (self.nextnode is not None and self.nextnode.status == "m"
            and sample() < self.K_PLUS2):
            # nextNode is methylated and sample() gets smaller than KPLUS2
            return UHistone(copy=True, copy_histone=self)
        return self


class HistoneOct4(Histone):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # copy histone info from copy_histone
        if kwargs.pop('copy', None):
            self.CpGislandlist = kwargs.pop('copy_histone', None).CpGislandlist

        # create new histone object based on the argument
        else:
            position = kwargs.pop('position', 0)
            if position == -2:
                self.CpGislandlist = [0, 0, 0, 0]
            elif position == -1:
                self.CpGislandlist = [0, 0, 0, 0]
            elif position == 0:
                self.CpGislandlist = [0, 0]
            else:
                self.CpGislandlist = []

    def k_minus(self):
        if sample() < self.K_MINUS:
            return UHistoneOct4(copy=True, copy_histone=self)
        else:
            return self

    def dna_methylation(self):
        if sample() < sum(self.CpGislandlist) * self.K_PLUS:
            # the chance gets higher, depends on the num of CpG island sites
            return MHistoneOct4(copy=True, copy_histone=self)
        else:
            return self


class MHistoneOct4(HistoneOct4):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.status = "m"


class UHistoneOct4(HistoneOct4):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.status = "u"

    def k_plus(self):
        if (self.prenode is not None and self.prenode.status == "m"
            and sample() < self.K_PLUS):
            # preNode is methylated and sample() gets smaller than K_PLUS
            return MHistoneOct4(copy=True, copy_histone=self)
        if (self.nextnode is not None and self.nextnode.status == "m"
            and sample() < self.K_PLUS):
            # nextNode is methylated and sample() gets smaller than KPLUS
            return MHistoneOct4(copy=True, copy_histone=self)
        return self

    def k_minus(self):
        return self

    def k_ace(self):
        # if histone has CpG island on it
        if 1 in self.CpGislandlist:
            # histone cannnot get acetylated
            return self

        else:
            if sample() < self.K_ACE:
                return AHistoneOct4(copy=True, copy_histone=self)
            else:
                return self


class AHistoneOct4(HistoneOct4):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.status = "a"

    def k_plus(self):
        if (self.prenode is not None and self.prenode.status == "m"
            and sample() < self.K_PLUS2):
            # preNode is methylated and sample() gets smaller than PLUS2
            return UHistoneOct4(copy=True, copy_histone=self)

        if (self.nextnode is not None and self.nextnode.status == "m"
            and sample() < self.K_PLUS2):
            # nextNode is methylated and sample() gets smaller  than PLUS2
            return UHistoneOct4(copy=True, copy_histone=self)

        return self


def init_genome(percentage=50,
                a_bool=1,
                hst_n=81,
                kp=0.12,
                kp2=0.12,
                km=0.117,
                ka=0.12):
    """
    percentage ... the probability of having methylated hitone.
    this method returns a list of histone randomly generated with respect to
    the inputs.
    """
    before_promoter = hst_n // 2
    hst_list = []  # hst_list stores histones
    ratio = percentage / 100  ## ratio should be float number between 0 and 1

    for i in range(hst_n):

        if sample() < ratio:
            hst_list.append(MHistone(position=i - before_promoter,
                                     kp=kp,
                                     kp2=kp2,
                                     km=km,
                                     ka=ka,
                                     a_bool=a_bool
                                     )
                            )
        else:
            hst_list.append(AHistone(position=i - before_promoter,
                                     kp=kp,
                                     kp2=kp2,
                                     km=km,
                                     ka=ka,
                                     a_bool=a_bool
                                     )
                            )

        hst_list[i - 1].set_adjhistone(hst_list[i])

    hst_list[0].prenode = None  # disjoint the edge histone to itself.
    return hst_list


def init_genome_oct4(percentage=50,
                                 a_bool=1,
                                 hst_n=81,
                                 kp=0.176,
                                 kp2=0.17,
                                 km=0.117,
                                 ka=0.12):
    """
    this method is a modified version of createRandomHistoneList
    used for Oct4 histones.
    """
    before_promoter = hst_n // 2
    hst_list = []  # hst_list stores histones
    ratio = percentage / 100  # ratio should be float number between 0 and 1

    for i in range(hst_n):

        if sample() < ratio:
            hst_list.append(MHistoneOct4(position=i - before_promoter,
                                         kp=kp,
                                         kp2=kp2,
                                         km=km,
                                         ka=ka,
                                         a_bool=a_bool
                                         )
                            )
        else:
            hst_list.append(AHistoneOct4(position=i - before_promoter,
                                         kp=kp,
                                         kp2=kp2,
                                         km=km,
                                         ka=ka,
                                         a_bool=a_bool
                                         )
                            )

        hst_list[i - 1].set_adjHistone(hst_list[i])

    hst_list[0].prenode = None  # disjoint the edge histone to itself.
    return hst_list


def next_genome(hst_list, a_bool, r_bool, window):
    """
    this method takes histone list and returns the next generation out of them.
    """
    nexthst_list = []
    ahst_n = 0
    mhst_n = 0

    for hst in hst_list:

        if -window // 2 <= hst.position <= window // 2:
            if hst.status == "a":
                ahst_n += 1
            elif hst.status == "m":
                mhst_n += 1

        hst = hst.k_minus()
        hst = hst.k_ace()
        hst = hst.k_plus()

        nexthst_list.append(hst)

    t_bool = a_bool and (ahst_n > 5)
    """
    WINDOW is size 10(11 histones note that there is E0 between E(-1) and E(1)), 
    so acetylated histones will be dominant if non-acetylated histones are less than 5.
    """
    if r_bool == 1:
        eext_bool = not t_bool
    else:
        eext_bool = mhst_n > 2

    if eext_bool:
        center = len(hst_list) // 2
        nexthst_list[center] = MHistone(copy=True, copy_histone=nexthst_list[center])

    return {"hstL": nexthst_list, "T": t_bool, "Eext": eext_bool}


def next_genome_oct4(hst_list, a_bool, r_bool, window, p_off):
    """
    this method takes histone list and returns the next generation of them.
    """
    nexthst_list = []
    ahst_n = 0
    mhst_n = 0

    for hst in hst_list:
        if -window // 2 <= hst.position <= window // 2:
            if hst.status == "a":
                ahst_n += 1
            elif hst.status == "m":
                mhst_n += 1

            for index in range(len(hst.CpGislandlist)):
                # p_on probability
                if hst.status == 'm' and sample() < 0.001133:
                    hst.CpGislandlist[index] = 1

                # p_off probability   
                if sample() < p_off:
                    hst.CpGislandlist[index] = 0

        hst = hst.k_ace()
        hst = hst.k_plus()
        hst = hst.k_minus()


        hst = hst.DNAmethylation()

        nexthst_list.append(hst)

    t_bool = a_bool and (ahst_n > window // 2)
    """
    WINDOW is size 10(11 histones note that there is E0 between E(-1) and E(1)), 
    so acetylated histones will be dominant if non-acetylated histones are less than 5.
    """
    if r_bool == 1:
        eext_bool = (not t_bool)
    else:
        eext_bool = mhst_n > 2

    if eext_bool:
        center = len(hst_list) // 2
        nexthst_list[center] = MHistoneOct4(copy=True, copy_histone=nexthst_list[center])

    return {"hstL": nexthst_list, "T": t_bool, "Eext": eext_bool}


def vectorize(hst_list):
    """
    this method takes a list of histone objects, and returns three dimention numpy array
    in which three bit vectors are stored. 
    """
    v_mlist = np.array([1 if h.status == "m" else 0 for h in hst_list])
    v_ulist = np.array([1 if h.status == "u" else 0 for h in hst_list])
    v_alist = np.array([1 if h.status == "a" else 0 for h in hst_list])

    return np.array([v_mlist,
                     v_ulist,
                     v_alist],
                    np.int32)


def vectorize_oct4(hst_list):
    v_mlist = np.array([1 if h.status == "m" else 0 for h in hst_list])
    v_ulist = np.array([1 if h.status == "u" else 0 for h in hst_list])
    v_alist = np.array([1 if h.status == "a" else 0 for h in hst_list])
    v_cpg = [sum(h.CpGislandlist) for h in hst_list]

    return np.array([v_mlist,
                     v_ulist,
                     v_alist,
                     v_cpg],
                    np.int32)


def track_epigenetic_process(hst_list,  # initial histone list
                             time,  # time for tracking (in hour)
                             a_bool,  # activator bool
                             r_bool,  # repressor bool
                             t_bool,  # transcription bool
                             window=10  # default is 10
                             ):
    for i in range(len(hst_list)):
        hst_list[i].set_ka(a_bool)
    vectorizedgene_list = []  # array of compressed data of vectors
    t_list = []  # one dimension array
    for _ in range(time):
        vectorizedgene_list.append(vectorize(hst_list))
        t_list.append(t_bool)
        dict_hst = next_genome(hst_list, a_bool, r_bool, window)
        hst_list = dict_hst["hstL"]
        t_bool = dict_hst["T"]

    return {"vectorize": np.array(vectorizedgene_list), "hstL": hst_list, "TList": t_list}


def track_epigenetic_process_oct4(hst_list,  # initial histone list
                                  time,  # time for tracking (in hour)
                                  a_bool,  # activator bool
                                  r_bool,  # repressor bool
                                  t_bool,  # transcription bool
                                  p_off,  # prob of CpG island site gets OFF
                                  window=10  # default is 10
                                  ):
    for i in range(len(hst_list)):
        hst_list[i].set_K_ACE(a_bool)
    vectorizedgene_list = []  # array of compressed data of vectors
    t_list = []  # one dimension array
    for _ in range(time):
        vectorizedgene_list.append(vectorize_oct4(hst_list))
        t_list.append(t_bool)
        dict_hst = next_genome_oct4(hst_list, a_bool, r_bool, window, p_off)
        hst_list = dict_hst["hstL"]
        t_bool = dict_hst["T"]

    return {"vectorize": np.array(vectorizedgene_list), "hstL": hst_list, "TList": t_list}


def gettime_mhst_dieout(vectorizedgene_list):
    for time in range(len(vectorizedgene_list)):
        if vectorizedgene_list[time][0][40] == 0:
            for else_in_window in range(35, 45):
                if vectorizedgene_list[time][0][else_in_window] == 1:
                    break
                return time
    return -1
