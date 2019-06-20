""" File with some auxiliary classes and functions
    - IQdata -  Complex data type for myHDL blocks, class
    - RTLblocks - Set of rtl-blocks writing in myHDL, class
"""

import numpy as np
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

    @staticmethod
    @block
    def SFIFO(i_clk, i_rst, i_wren, i_d, i_rden, o_d, o_dv, o_empty, o_full, o_ovflo, DEPTH):
        """ Synchronous (one-clock domain) FIFO (first word fall through type)
        ~~~~ input ports ~~~~
        :param i_clk: clock;
        :param i_rst: synchronous reset for FIFO state;
        :param i_wren: write enable;
        :param i_d: writing data;
        :param i_rden: read enable;
        ~~~~ output ports ~~~~
        :param o_d: read data;
        :param o_dv: output data valid flag;
        :param o_empty: empty flag;
        :param o_full: full flag;
        :param o_ovflo: overflow flag;
        ~~~~ parameters ~~~~
        :param DEPTH: depth of FIFO, real depth will be 2**k where k = np.ceil(np.log2(DEPTH))
        :return:
        """
        tState = enum('EMPTY', 'VALID', 'FULL', 'OVFLO', encoding="binary")

        ADR_WIDTH = int(np.ceil(np.log2(DEPTH)))
        mem = [Signal(i_d.val) for _ in range(DEPTH)]
        wrPnt = Signal(modbv(0)[ADR_WIDTH:])
        rdPnt = Signal(modbv(0)[ADR_WIDTH:])
        state = Signal(tState.EMPTY)
        val = Signal(bool(0))
        empty = Signal(bool(1))
        full = Signal(bool(0))
        ovflo = Signal(bool(0))

        @always(i_clk.posedge)
        def wrToRAM():
            if i_wren:
                mem[wrPnt].next = i_d

        @always_comb
        def rdFromRam():
            o_d.next = mem[rdPnt]

        rdDifWr = Signal(modbv(wrPnt.val))
        wrDifRd = Signal(modbv(wrPnt.val))

        @always_comb
        def pntDif():
            rdDifWr.next = rdPnt - wrPnt
            wrDifRd.next = wrPnt - rdPnt

        @always_seq(i_clk.posedge, reset=i_rst)
        def FSM():
            if state == tState.EMPTY:
                val.next = bool(0)
                empty.next = bool(1)
                full.next = bool(0)
                ovflo.next = bool(0)
                if i_wren:
                    wrPnt.next = wrPnt + 1
                    state.next = tState.VALID
                    val.next = bool(1)
                    empty.next = bool(0)
            elif state == tState.VALID:
                if i_wren:
                    wrPnt.next = wrPnt + 1
                    if i_rden:
                        rdPnt.next = rdPnt + 1
                    else:
                        if rdDifWr == 1:
                            state.next = tState.FULL
                            full.next = bool(1)
                else:
                    if i_rden:
                        rdPnt.next = rdPnt + 1
                        if wrDifRd == 1:
                            state.next = tState.EMPTY
                            val.next = bool(0)
                            empty.next = bool(1)
            elif state == tState.FULL:
                if i_wren:
                    wrPnt.next = wrPnt + 1
                    if not i_rden:
                        state.next = tState.OVFLO
                        val.next = bool(0)
                        ovflo.next = bool(1)
                    else:
                        rdPnt.next = rdPnt + 1
                else:
                    if i_rden:
                        rdPnt.next = rdPnt + 1
                        state.next = tState.VALID
                        full.next = bool(0)
            elif state == tState.OVFLO:
                ovflo.next = bool(1)
                full.next = bool(1)
                val.next = bool(0)

        @always_comb
        def outputs():
            o_dv.next = val
            o_empty.next = empty
            o_full.next = full
            o_ovflo.next = ovflo

        return instances()

if __name__ == "__main__":

    ref = {'d': [], 'empty': True, 'full': False, 'ovflo': False, 'len': 8}

    def readFrom(fifo):
        if fifo['empty'] or fifo['ovflo']:
            return None
        else:
            val = fifo['d'].pop(0)
            if len(fifo['d']) == 0:
                fifo['empty'] = True
            if fifo['full']:
                fifo['full'] = False
            return val


    def wrTo(fifo, val):
        if fifo['full']:
            fifo['ovflo'] = True
        if len(fifo['d']) == fifo['len'] - 1:
            fifo['full'] = True
        if fifo['empty']:
            fifo['empty'] = False
        if not fifo['ovflo']:
            fifo['d'].append(val)


    def rstRef():
        global ref
        ref = {'d': [], 'empty': True, 'full': False, 'ovflo': False, 'len': 8}

    @block
    def testbench(bw=8, autoCheck=True):
        period = 10
        i_clk = Signal(bool(0))
        i_rst = ResetSignal(0, active=bool(1), isasync=False)
        i_wren = Signal(bool(0))
        i_d = Signal(intbv(0, min=-2 ** (bw - 1), max=2 ** (bw - 1)))
        i_rden = Signal(bool(0))
        o_d = Signal(i_d.val)
        o_dv = Signal(bool(0))
        o_empty = Signal(bool(0))
        o_full = Signal(bool(0))
        o_ovflo = Signal(bool(0))
        DEPTH = 8

        uut = RTLblocks.SFIFO(i_clk, i_rst, i_wren, i_d, i_rden, o_d, o_dv, o_empty, o_full, o_ovflo, DEPTH)

        @always(delay(int(period // 2)))
        def clk_gen():
            i_clk.next = not i_clk

        @instance
        def stim():
            rstRef()
            yield i_clk.negedge
            i_rst.next = bool(1)
            yield i_clk.negedge
            i_rst.next = bool(0)
            if autoCheck:
                assert o_empty.val == ref['empty'], "Empty signals mismatch"
                assert o_full.val == ref['full'], "Full signals mismatch"
                assert o_ovflo.val == ref['ovflo'], "Overflow signals mismatch"
            print('write 10 values - overflow check')
            for ii in range(10):
                val = np.random.randint(-2 ** (bw - 1), 2 ** (bw - 1) - 1)
                i_d.next = val
                i_wren.next = bool(1)
                yield i_clk.negedge
                wrTo(ref, val)
                if autoCheck:
                    assert o_empty.val == ref['empty'], f"Empty signals mismatch, ii = {ii}"
                    assert o_full.val == ref['full'], f"Full signals mismatch, ii = {ii}"
                    assert o_ovflo.val == ref['ovflo'], f"Overflow signals mismatch, ii = {ii}"
                    if not (ref['empty'] or ref['ovflo']):
                        assert o_dv.val == True, f"Valid flag must be True, ii = {ii}"
            i_wren.next = bool(0)
            # reset
            rstRef()
            i_rst.next = bool(1)
            yield i_clk.negedge
            i_rst.next = bool(0)
            print('write 8 values - to full state')
            for ii in range(8):
                val = np.random.randint(-2 ** (bw - 1), 2 ** (bw - 1) - 1)
                i_d.next = val
                i_wren.next = bool(1)
                yield i_clk.negedge
                wrTo(ref, val)
                if autoCheck:
                    assert o_empty.val == ref['empty'], f"Empty signals mismatch, ii = {ii}"
                    assert o_full.val == ref['full'], f"Full signals mismatch, ii = {ii}"
                    assert o_ovflo.val == ref['ovflo'], f"Overflow signals mismatch, ii = {ii}"
                    if not (ref['empty'] or ref['ovflo']):
                        assert o_dv.val == True, f"Valid flag must be True, ii = {ii}"
                        assert o_d.val == ref['d'][0]
            i_wren.next = bool(0)
            print('read all values')
            for ii in range(8):
                i_rden.next = bool(1)
                yield i_clk.negedge
                readFrom(ref)
                if autoCheck:
                    assert o_empty.val == ref['empty'], f"Empty signals mismatch, ii = {ii}"
                    assert o_full.val == ref['full'], f"Full signals mismatch, ii = {ii}"
                    assert o_ovflo.val == ref['ovflo'], f"Overflow signals mismatch, ii = {ii}"
                    if not (ref['empty'] or ref['ovflo']):
                        assert o_dv.val == True, f"Valid flag must be True, ii = {ii}"
                        assert o_d.val == ref['d'][0]
            print(' **** Test passed **** ')
            raise StopSimulation

        return instances()

    tb = testbench(bw=8, autoCheck=True)
    tb.config_sim(trace=False)
    tb.run_sim()
