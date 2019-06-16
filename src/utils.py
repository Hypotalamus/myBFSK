""" File with some auxiliary classes and functions
    - IQdata -  Complex data type for myHDL blocks, class
    - RTLblocks - Set of rtl-blocks writing in myHDL, class
"""

from myhdl import *

class IQdata(object):
    """ Complex data type for myHDL blocks """

    def __init__(self, bWidth):
        self.I = Signal(intbv(0, min=-2 ** (bWidth - 1), max=2 ** (bWidth - 1)))
        self.Q = Signal(intbv(0, min=-2 ** (bWidth - 1), max=2 ** (bWidth - 1)))


class RTLblocks(object):
    """ Set of rtl-blocks writing in myHDL
        List of methods:
        - SingleAsRdRAM - Single-Port RAM with Asynchronous Read;
        - MAC - Multiply-accumulate block;
        - CMult - Complex multiplier
    """

    @staticmethod
    @block
    def SingleAsRdRAM(i_clk, i_d, i_wren, i_addr, o_d):
        """ Single-Port RAM with Asynchronous Read
            ~~~~ input ports ~~~~
            - i_clk - clock;
            - i_d - data to write;
            - i_wren - write enable;
            - i_addr - address for read/write;
            ~~~~ output ports ~~~~
            - o_d - read data
        """
        mem = [Signal(intbv(0, min=i_d.min, max=i_d.max)) for _ in
               range(i_addr.max)]

        @always(i_clk.posedge)
        def write():
            if i_wren:
                mem[i_addr].next = i_d

        @always_comb
        def read():
            o_d.next = mem[i_addr]

        return instances()

    @staticmethod
    @block
    def MAC(i_clk, i_rst, i_d1, i_d2, i_dv, i_accClr, o_d, o_acc, o_dv, wAcc=None):
        """ Multiply-accumulate block
            ~~~~ input ports ~~~~
            - i_clk - clock;
            - i_rst - synchronous reset;
            - i_d1 - input data 1st port;
            - i_d2 - input data 2nd port;
            - i_dv - input data valid flag;
            - i_accClr - clear accumulator;
            ~~~~ output ports ~~~~
            - o_d - output data (before accumulator's register);
            - o_acc - accumulator value, can be set None
            - o_dv - output data valid (== i_dv), can be set None
            ~~~~ parameters ~~~~
            - wAcc - bitwidth of accumulator. If None then wAcc = bitwidth of multipliers;
        """
        multMin = -1 * i_d1.min * i_d2.min
        multMax = i_d1.max * i_d2.max
        mult = Signal(intbv(0, min=multMin, max=multMax))
        if wAcc is None:
            accIn = Signal(mult.val)
        else:
            accIn = Signal(intbv(0, min=-2 ** (wAcc - 1), max=2 ** (wAcc - 1)))
        acc = Signal(accIn.val)

        @always(i_clk.posedge)
        def regs():
            if i_rst:
                mult.next = 0
                acc.next = 0
            else:
                if i_dv:
                    mult.next = i_d1 * i_d2
                    if i_accClr:
                        acc.next = 0
                    else:
                        acc.next = accIn

        @always_comb
        def intComb():
            accIn.next = acc + mult

        @always_comb
        def output():
            o_d.next = accIn

        if o_acc is not None:
            @always_comb
            def out_opt1():
                o_acc.next = acc

        if o_dv is not None:
            @always_comb
            def out_opt2():
                o_dv.next = i_dv

        return instances()

    @staticmethod
    @block
    def CMult(i_clk, i_rst, i_d1, i_d2, i_dv, o_d, o_dv, shOut):
        """ Complex multiplication block using Karatsuba algorithm (consume three multipliers and 5 adders)
        ~~~~ input ports ~~~~
        :param i_clk: clock;
        :param i_rst: synchronous reset;
        :param i_d1: input data 1, IQdata();
        :param i_d2: input data 2, IQdata();
        :param i_dv: input data valid flag;
        ~~~~ output ports ~~~~
        :param o_d: output data, IQdata();
        :param o_dv: output data valid flag, can be set None;
        ~~~~ parameters ~~~~
        :param shOut: right-shift output - scaling by 1/2**shOut, integer
        :return:
        """
        intMax = i_d1.I.max * i_d2.I.max * 8
        intMin =  -1 * i_d1.I.min * i_d2.I.min * 8
        d1I, d1I_reg = [Signal(i_d1.I.val) for _ in range(2)]
        d1Q, d1Q_reg = [Signal(i_d1.Q.val) for _ in range(2)]
        d2I, d2I_reg = [Signal(i_d2.I.val) for _ in range(2)]
        d2Q, d2Q_reg = [Signal(i_d2.Q.val) for _ in range(2)]
        s1, s2, p1, p2, p3, p3_reg, s3, s3_reg, s4, s5 = [Signal(intbv(0, min=intMin, max=intMax)) for _ in range(10)]

        @always_comb
        def indata():
            d1I.next = i_d1.I
            d1Q.next = i_d1.Q
            d2I.next = i_d2.I
            d2Q.next = i_d2.Q

        @always(i_clk.posedge)
        def regs():
            if i_rst:
                d1I_reg.next = 0
                d1Q_reg.next = 0
                d2I_reg.next = 0
                d2Q_reg.next = 0
                s1.next = 0
                s2.next = 0
                p1.next = 0
                p2.next = 0
                p3.next = 0
                p3_reg.next = 0
                s3.next = 0
                s3_reg.next = 0
                s4.next = 0
                s5.next = 0
            else:
                if i_dv:
                    d1I_reg.next = d1I
                    d1Q_reg.next = d1Q
                    d2I_reg.next = d2I
                    d2Q_reg.next = d2Q
                    s1.next = d1I + d1Q
                    s2.next = d2I + d2Q
                    p1.next = d1I_reg * d2I_reg
                    p2.next = d1Q_reg * d2Q_reg
                    p3.next = s1 * s2
                    p3_reg.next = p3
                    s3.next = p1 - p2
                    s3_reg.next = s3
                    s4.next = p1 + p2
                    s5.next = p3_reg - s4

        @always_comb
        def output():
            o_d.I.next = s3_reg >> shOut
            o_d.Q.next = s5 >> shOut

        if o_dv is not None:
            @always_comb
            def out_opt():
                o_dv.next = i_dv

        return instances()